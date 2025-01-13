% Example OES calculations using 4 and 9 vessels, with and without off-resonance corrections applied

%% Set parameters that stay the same for all simulations
motion = 4; 
pad_size = 1024;
vesloc3D = 0;
AngTowardsCor = 0;
AngTowardsSag = 0;
off_res_flag = 1; 
debug = true; % This means the encodings in each example will be visualised below
had_enc = 0;

%% Test 4 vessel with no off-resonance
vesloc = [-25 25 -25 25; 25 25 -25 -25]; % Simple four vessel geometry to start
standard_enc = 1; PCASL_enc = 0; standard_PCASL = 0; off_res_flag = 0; off_res = zeros(1,size(vesloc,2)); offset = [];
% When off-res flag is false, the OES calculation is as described in Berry
% MRM 2015.
[k_coords,frequency,phase,e,h,condition,SNR,ang_locs_tagmode] = ...
    OES_offres(vesloc,motion,pad_size,had_enc,standard_enc,PCASL_enc,standard_PCASL, ...
        vesloc3D,AngTowardsCor,AngTowardsSag,off_res_flag,off_res,debug);
    
LargeFigWindow(0.3,0.5); subplot(1,2,1);
imagesc(h); title(['Standard encoding: requested']); caxis([-1 1]); axis equal; axis off; colorbar;
subplot(1,2,2); 
imagesc(e); title(['Standard encoding: achieved. Mean SNR = ' num2str(SNR)]); caxis([-1 1]); axis equal; axis off; colorbar;
colormap gray

% Note that even for this near-perfect encoding, the relative SNR is less than one,
% due to the calculation accounting for imperfect labelling efficiency of
% PCASL.

%% Repeat with off_res_flag on
% Should give near-identical results when there is no off-resonance,
% although the way the encoding is represented and the range of encoding
% directions available is slightly different within the OES algorithm.
off_res_flag = 1;
[k_coords,frequency,phase,e,h,condition,SNR,ang_locs_tagmode] = ...
    OES_offres(vesloc,motion,pad_size,had_enc,standard_enc,PCASL_enc,standard_PCASL, ...
        vesloc3D,AngTowardsCor,AngTowardsSag,off_res_flag,off_res,debug);
    
LargeFigWindow(0.3,0.7); subplot(1,2,1);
imagesc(h); title(['Off-res flag on: requested']); caxis([-1 1]); axis equal; axis off; colorbar;
subplot(1,2,2); 
imagesc(e); title(['Off-res flag on: achieved. Mean SNR = ' num2str(SNR)]); caxis([-1 1]); axis equal; axis off; colorbar;
colormap gray;

%% Try with some fixed off-resonance
off_res = ones(1,size(vesloc,2))*pi; had_enc = 0; standard_enc = 1; PCASL_enc = 0;
[k_coords,frequency,phase,e,h,condition,SNR,ang_locs_tagmode] = ...
    OES_offres(vesloc,motion,pad_size,had_enc,standard_enc,PCASL_enc,standard_PCASL, ...
        vesloc3D,AngTowardsCor,AngTowardsSag,off_res_flag,off_res,debug);

LargeFigWindow(0.3,0.7); subplot(1,2,1);
imagesc(h); title(['Fixed off-res: requested']); caxis([-1 1]); axis equal; axis off; colorbar;
subplot(1,2,2); 
imagesc(e); title(['Fixed off-res: achieved. Mean SNR = ' num2str(SNR)]); caxis([-1 1]); axis equal; axis off; colorbar;
colormap gray;

%% Try with some variable off-resonance
off_res = [0 0 1 1]*pi/4; had_enc = 0; standard_enc = 1; PCASL_enc = 0; standard_PCASL = 0;
[k_coords,frequency,phase,e,h,condition,SNR,ang_locs_tagmode] = ...
    OES_offres(vesloc,motion,pad_size,had_enc,standard_enc,PCASL_enc,standard_PCASL, ...
        vesloc3D,AngTowardsCor,AngTowardsSag,off_res_flag,off_res,debug);


LargeFigWindow(0.3,0.7); subplot(1,2,1);
imagesc(h); title(['Variable off-res: requested']); caxis([-1 1]); axis equal; axis off; colorbar;
subplot(1,2,2); 
imagesc(e); title(['Variable off-res: achieved. Mean SNR = ' num2str(SNR)]); caxis([-1 1]); axis equal; axis off; colorbar;
colormap gray;

%% Test 9 vessel, no off-resonance
vesloc = [-3 -49 -46 -46 37 36 34 -25 9; 39.2075 34.7354 21.3191 8.7994 24.4517 10.1389 2.0912 -7.7494 -9.9855];
off_res = zeros(1,size(vesloc,2));
had_enc = 1; standard_enc = 0; PCASL_enc = 0; standard_PCASL = 0;
[k_coords,frequency,phase,e,h,condition,SNR,ang_locs_tagmode] = ...
    OES_offres(vesloc,motion,pad_size,had_enc,standard_enc,PCASL_enc,standard_PCASL, ...
        vesloc3D,AngTowardsCor,AngTowardsSag,off_res_flag,off_res,debug);

LargeFigWindow(0.3,0.7); subplot(1,2,1);
imagesc(h); title(['9 vessel standard encoding: requested']); caxis([-1 1]); axis equal; axis off; colorbar;
subplot(1,2,2); 
imagesc(e); title(['9 vessel standard encoding: achieved. Mean SNR = ' num2str(SNR)]); caxis([-1 1]); axis equal; axis off; colorbar;
colormap gray;

%% Test 9 vessel, with variable off-resonance
vesloc = [-3 -49 -46 -46 37 36 34 -25 9; 39.2075 34.7354 21.3191 8.7994 24.4517 10.1389 2.0912 -7.7494 -9.9855];
off_res = [0 0 0 0 1 1 1 1 1]*pi/4;
had_enc = 1; standard_enc = 0; PCASL_enc = 0; standard_PCASL = 0;
[k_coords,frequency,phase,e,h,condition,SNR,ang_locs_tagmode] = ...
    OES_offres(vesloc,motion,pad_size,had_enc,standard_enc,PCASL_enc,standard_PCASL, ...
        vesloc3D,AngTowardsCor,AngTowardsSag,off_res_flag,off_res,debug);

LargeFigWindow(0.3,0.7); subplot(1,2,1);
imagesc(h); title(['9 vessel variable off-res: requested']); caxis([-1 1]); axis equal; axis off; colorbar;
subplot(1,2,2); 
imagesc(e); title(['9 vessel variable off-res: achieved. Mean SNR = ' num2str(SNR)]); caxis([-1 1]); axis equal; axis off; colorbar;
colormap gray;
