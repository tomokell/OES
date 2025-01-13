%***********************************************************************
% Run the optimal encoding scheme (OES) algorithm to determine optimal encodings 
% for vessel-encoded PCASL, or conventional PCASL, with or without off-resonance 
% correction.
%
% Original code: Eleanor SK Berry, 2018
% Updated by: Thomas W Okell, 2025
%***********************************************************************
%
% [k_coords,frequency,phase,e,hadamard_complete,condition,SNR,ang_locs_tagmode] = ...
%               OES_offres(vesloc,motion,pad_size,had_enc,standard_enc,PCASL_enc, ...
%                          standard_PCASL,vesloc3D,AngTowardsCor,AngTowardsSag,   ...
%                          off_res_flag,off_res,debug)
%
% INPUT ARGUMENTS:
% vesloc         : Matrix containing the coordinates of the vessels in millimeters [x1, x2...; y1 y2 ... <; z1 z2...>]
% motion         : Scalar value indicating the expected motion magnitude in mm, used to 
%                  calculate the encoding frequency mask.
% pad_size       : Integer representing the size of the padded image matrix used within the OES calculation.
% had_enc        : Flag to state whether to use a standard hadamard encoding matrix as the ideal encoding scheme.
% standard_enc   : Flag to indicate whether a standard encoding scheme is being used (4 vessels only),
%                  as described in Okell et al. JCBFM 2013.
% PCASL_enc      : Flag for standard PCASL encoding (i.e. non-selective tag and control only)
% standard_PCASL : Set to true for standard PCASL with no off-resonance correction
% vesloc3D       : Boolean indicating whether the input vessel locations are in 3D.
% AngTowardsCor  : Angle defining the rotation towards the coronal plane for converting 
%                  3D vessel locations to 2D.
% AngTowardsSag  : Angle defining the rotation towards the sagittal plane for converting 
%                  3D vessel locations to 2D.
% off_res_flag   : Boolean indicating whether off-resonance correction is applied.
% off_res        : Vector containing off-resonance values at each vessel location.
% debug          : Boolean flag to enable or disable debug visualizations.
%
% OUTPUT ARGUMENTS:
% k_coords       : Row/column coordinates within the padded matrix of the optimum k-space point
% frequency      : Optimal spatial frequency vectors for each encoding cycle in rads/mm.
% phase          : Phase of each optimal encoding.
% e              : Final achieved encoding matrix, using a realistic modulation function.
% hadamard_complete : Complete ideal Hadamard matrix used for encoding calculations.
% condition      : Condition number of the encoding matrix
% SNR            : Relative signal-to-noise ratio for this encoding, as per
%                  Wong MRM 2007. Note that this includes imperfect
%                  labelling efficiency, so even an ideal encoding will result in SNR
%                  less than one.
% ang_locs_tagmode : Matrix containing tagging angles, shifts from isocenter, and tag modes that can be
%                    saved to a text file and passed to the scanner.

function [k_coords,frequency,phase,e,hadamard_complete,condition,SNR,ang_locs_tagmode] = ...
            OES_offres(vesloc,motion,pad_size,had_enc,standard_enc,PCASL_enc,standard_PCASL,vesloc3D,AngTowardsCor,AngTowardsSag,off_res_flag,off_res,debug)


%*********Set up for encoding calculation***************
matrix_dim = pad_size; %padded pixel matrix size

%convert vessel locations from mm to pixels
[pixves] = convert_vessels_mm2pix(vesloc,pad_size);

%Convert 3D vessels to 2D in a single plane
if vesloc3D == true
    [LRUnitVec, APUnitVec] = CalcObliqueUnitVecs(AngTowardsCor,AngTowardsSag);
    vesloc = Convert3DTo2DVesLocs(vesloc, LRUnitVec, APUnitVec);
end

%Use the ideal encoding scheme set by user
[hadamard_complete, hadamard_for_OES, TagMode] = calc_hadamard_matrices(vesloc,had_enc,standard_enc,PCASL_enc,off_res_flag,off_res);

%Pre allocate variable names to speed up processing time in matlab
r = zeros(1,size(hadamard_for_OES,1));
c = zeros(1,size(hadamard_for_OES,1));
frequency = zeros(size(hadamard_for_OES,1),2);
phase = zeros(1,size(hadamard_for_OES,1));
encoding = zeros(size(hadamard_for_OES,1),size(vesloc,2));
pp = zeros(size(hadamard_for_OES,1),size(vesloc,2));

%define some constants
units = (2*pi/matrix_dim);                                      %units for converting coordinates in to spatial frequencies
k_centre = [(matrix_dim/2)+1 (matrix_dim/2)+1];                 %coordinates of the centre of kspace
delta_k_padded = 1/(matrix_dim);                                %delta k in 1/mm of the zero padded matrix
lamda_limit = motion*4;                                         %minimum possible wavelength [mm]
max_mask_radius = 1/(delta_k_padded*lamda_limit);               %maximum radius of the mask, depends on the minimum possible wavelength [#pixels = 1/((mm^-1)*mm)]
half_FOV_padded = (matrix_dim)/2;
xy = -half_FOV_padded:(2*half_FOV_padded)/matrix_dim:half_FOV_padded; %make a vector (in mm) to fill the axes of later figures
half_kFOV_padded = (matrix_dim)/2;
kxy = -(half_kFOV_padded*delta_k_padded):delta_k_padded:(half_kFOV_padded*delta_k_padded);

%Create the mask and weighting function:
%Circular mask, based on distance of each point from the centre and 
%limited to a range of lower frequencies according to the predicted
%amount of subject motion (or an imposed limit if that range of frequencies is lower)
position_weighting = zeros(matrix_dim,matrix_dim);
mask = false(matrix_dim,matrix_dim);
for g = 1:matrix_dim
    for f = 1:matrix_dim
        position_weighting(g,f) = sqrt((g-k_centre(1,1))^2 + (f-k_centre(1,2))^2);
        if position_weighting(g,f) <= max_mask_radius
            mask(g,f) = true;
        end
    end
end

% Remove the very centre of k-space to avoid selecting zero spatial
% frequency which causes issues in the export of parameters to the sequence
% below
mask(k_centre(1),k_centre(2)) = false;

%normalise the position weighting and define it so that the
%maximum value is 0.5 (at the centre of k space)
p = (-position_weighting/max(max(position_weighting)) + ones(matrix_dim,matrix_dim))/2;

%only remove half of k space if the off resonance correction is not being used
%(k space is not symmetric when calculating it with complex numbers)
if off_res_flag == false   
    mask(1:matrix_dim,1:(matrix_dim/2)) = 0;                 %make all of -kx = 0
    mask(k_centre(1)+1:matrix_dim,k_centre(2)) = 0;  %make -ky below the centre of k space = 0
end



%%
%**********Encoding calculation***************
for m = 1:size(hadamard_for_OES,1)
    
    % Check whether OES calculation is needed for this cycle
    if ((TagMode(m) < 2) && (off_res_flag == false)) || standard_PCASL % Non-selective cycle without off-resonance correction: OES not needed
        frequency(m,1:2) = [0 0];
        phase(m) = (mod(TagMode(m)+1,2)*pi);   % Tag mode = 0 is control: want phase = pi, tag mode 1 is tag: want phase = 0


    else % OES required

        % Step 0: Flag this as now a selective cycle, since we want to use
        % transverse gradients, even if it's a tag all/control all
        % condition
        TagMode(m) = 2;

        %Step 1:
        %turn the vessel locations into +/-1 or exp(i*phi) for off-resonance
        %correction
        Mv = zeros(matrix_dim,matrix_dim);
        for n = 1:size(vesloc,2)
            Mv(pixves(2,n),pixves(1,n)) = hadamard_for_OES(m,n);
        end
    
        %Step 2:
        %Shift the vessels to move the centre from the middle to the edge 
        %in order to do the matlab FT, FT the 'vessels' and then shift the FT back, 
        %so k=0 is at the centre of the image
        k_vessel_centre=fftshift(fft2(ifftshift(Mv)));
    
        %Step 3:
        %normalise k space
        k = abs(k_vessel_centre)/max(max(abs(k_vessel_centre)));
          
        %Step 4:
        %Add the normalised position and kspace values so the maxima closest to the centre will be chosen first 
        weighted_kspace = p + k;
    
        %Step 5:
        %mask normalised k space data
        weighted_kspace = weighted_kspace.*mask;
        if m == 3 % Save out one encoding as an example
            wk = weighted_kspace;
        end
        
        %Step 6:
        %find the maximum in k space closest to the centre and work out the
        %frequency and the phase at this point
        [r(m),c(m)] = find(weighted_kspace == max(max(weighted_kspace)),1);
        frequency(m,1:2) = ([c(m) r(m)] - k_centre)*units; 
        phase(m) = angle(k_vessel_centre(r(m),c(m)));
    
        %Step 6a:
        %add pi to the phi term if using the off resonance correction as the ideal terms are based on desired phase, not +/-1
        if off_res_flag
            phase(m) = phase(m) + pi;   
        end    

    end

    %Step 7:
    %calculate the encodings using this frequency and phase 
    %encoding based on sinusoid = real(e^i(f*x + phi)): encoding(m,d) = real(exp(1i*(frequency(m,1:2)*vesloc(1:2,d) + phase(m))));
    %encoding based on 3 Fourier coeeficients from Wong 2012: 
    %         delta_Mz = 2.092*cos(pp) - 0.322*cos(3*pp) + 0.053*cos(5*pp);
    %         encoding(m,d) = delta_Mz/2;
    %encoding based on the simulated modmat curve:
    encoding(m,:) = encodings_offset(vesloc,frequency(m,1:2),phase(m),off_res);       

    if debug
        if m == 1; figure; end
        prows = ceil(sqrt(size(hadamard_for_OES,1)));
        pcols = ceil(size(hadamard_for_OES,1) / prows);
        subplot(prows,pcols,m)
        VisualiseVesselEncodingOES(matrix_dim,frequency(m,:),phase(m),vesloc,hadamard_complete(m,:),true,off_res);
        title(['Encoding ' num2str(m)])
    end
end


%**************Convert encodings in to phase angles and wavelengths******

k_coords = [c;r];

% Add in static tissue to the encoding matrix
static_tissue = ones(size(hadamard_for_OES,1),1);

if PCASL_enc
    % For standard PCASL, we don't want to separate out the vessels, so take
    % the average encoding of all vessels
    e = [mean(encoding,2),static_tissue];
else % Vessel-specific encoding
    e = [encoding,static_tissue];
end

%calculate condition number and Wong SNR efficiency of the encodings
condition = cond(e);
[~,SNR] = WongSNR(e);

%calculate the tagging angles and the shifts from isocentre (in mm)
k_mm = frequency/(2*pi);
tagang = zeros(1,size(frequency,1));
ab = zeros(1,size(frequency,1));
encs = zeros(2,size(frequency,1));
for n = 1:size(frequency,1)
    % For non-selective cycles, use arbitrary a/b
    if TagMode(n) < 2
        tagang(n) = 0;
        ab(n) = 100;
        half_lambda = 200;
    else % Selective cycle
        tagang(n) = torad2deg( atan2(k_mm(n,1),k_mm(n,2)) );
        ab(n) = -phase(n)/(2*pi*sqrt(k_mm(n,1).^2 + k_mm(n,2).^2));
        half_lambda = 1/(2*sqrt(k_mm(n,1).^2 + k_mm(n,2).^2));
    end

    encs(:,n) = [ab(n)-half_lambda, ab(n)];
end

ang_locs_tagmode = [tagang;encs;TagMode]';


return