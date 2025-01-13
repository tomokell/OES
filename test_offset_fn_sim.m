% Test the OES off-resonance correction by comparing no off-resonance with,
% off-resonance and no correction vs. off-resonance with correction using a
% simulated fieldmap with linear and/or quadratic terms.
% NB. The vessel locations need to be integers to work here

function [meanTotalSNR, stdTotalSNR,h_value, p_value] = test_offset_fn_sim(linear_offset, linear_offset_dirn, quadratic_offset,...
    vesloc,motion,pad_size,had_enc,standard_enc,PCASL_enc,vesloc3D,AngTowardsCor,AngTowardsSag,debug,set_iterations)

if nargin < 14; set_iterations = 20; end

%initialise standard_PCASL factor to zero
standard_PCASL = 0;

% Simulate the fieldmap
fieldmap = simulate_fieldmap(pad_size, linear_offset, linear_offset_dirn, quadratic_offset);

%convert vessels into pixel locations in order to calculate the encoding
%functions with and without field offsets
[pixves] = convert_vessels_mm2pix(vesloc,pad_size);
fieldmap_vesloc = vesloc + 32;

%preassign some variables
meanSNR_enc = zeros(1,set_iterations);
meanSNR_offsetenc = zeros(1,set_iterations);
meanSNR_offsetCorrenc = zeros(1,set_iterations);
no_off_res = zeros(1,size(vesloc,2));

% Encoding calculations
off_res_range = zeros(set_iterations,size(vesloc,2));
for s = 1:set_iterations
    
    disp(['test_offset_fn: Iteration ' num2str(s) ' of ' num2str(set_iterations)]);
    
    %find the off resonance at the vessel locations (needs to be in radians)
    for n = 1:size(vesloc,2)
        off_res(n) = fieldmap(pixves(2,n),pixves(1,n))*s;
    end
    off_res_range(s,:) = off_res;

    %find the encoding functions without any offset present
    if PCASL_enc == true
        standard_PCASL = 1;
    end

    % No off-resonance
    %[k_coords,frequency,phase,e,hadamard_complete,condition,SNR,ang_locs_tagmode]
    off_res_flag = false;
    [~,~,~,~,~,~,meanSNR_enc(s),~] = ...
            OES_offres(vesloc,motion,pad_size,had_enc,...
            standard_enc,PCASL_enc,standard_PCASL,vesloc3D,AngTowardsCor,AngTowardsSag,off_res_flag,no_off_res,debug);
        
    
    % Off-resonance but no correction
    [~,~,~,~,~,~,meanSNR_offsetenc(s),~] = ...
            OES_offres(vesloc,motion,pad_size,had_enc,...
            standard_enc,PCASL_enc,standard_PCASL,vesloc3D,AngTowardsCor,AngTowardsSag,off_res_flag,off_res,debug);

    % With offset correction
    off_res_flag = true;
    if PCASL_enc == true   %re-set the standard PCASL setting to zero for the corrected versions
        standard_PCASL = 0;
    end
    [~,~,~,~,~,~,meanSNR_offsetCorrenc(s),~] = ...
            OES_offres(vesloc,motion,pad_size,had_enc,...
            standard_enc,PCASL_enc,standard_PCASL,vesloc3D,AngTowardsCor,AngTowardsSag,off_res_flag,off_res,debug);    
    
end

[h_value(1), p_value(1)] = ttest(meanSNR_enc,meanSNR_offsetenc);
[h_value(2), p_value(2)] = ttest(meanSNR_enc,meanSNR_offsetCorrenc);
[h_value(3), p_value(3)] = ttest(meanSNR_offsetenc,meanSNR_offsetCorrenc);

meanTotalSNR = [mean(meanSNR_enc) mean(meanSNR_offsetenc) mean(meanSNR_offsetCorrenc)];
stdTotalSNR = [std(meanSNR_enc) std(meanSNR_offsetenc) std(meanSNR_offsetCorrenc)];

return
