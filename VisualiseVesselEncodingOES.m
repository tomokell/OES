% Visualise vessel-encodings from OES calculations. Note that purely for
% visualisation purposes, the code below adds a smoothly varying fieldmap
% that matches the off-resonance values at the vessel locations, but does
% not necessarily represent the true underlying fieldmap at other
% locations.

function VisualiseVesselEncodingOES(matrix_dim,frequency,phase,vesloc,ideal_enc_mtx_row,ignore_final_col,off_res)

if nargin < 7; off_res = []; end

if ignore_final_col % Ignore final column (e.g. static tissue)
    ideal_enc_mtx_row = ideal_enc_mtx_row(1:end-1);
end

debug = false;
units = (2*pi/matrix_dim);                                      %units for converting coordinates in to spatial frequencies
k_centre = [(matrix_dim/2)+1 (matrix_dim/2)+1];                 %coordinates of the centre of kspace
half_FOV_padded = (matrix_dim)/2;
xy = -half_FOV_padded:(2*half_FOV_padded)/matrix_dim:(half_FOV_padded-(2*half_FOV_padded)/matrix_dim); %make a vector (in mm) to fill the axes of later figures

% Find the optimal point in the matrix
cr = round(frequency / units + k_centre);

test_freq = zeros(matrix_dim,matrix_dim);
test_freq(cr(2),cr(1)) = exp(1i*phase); % Apply unit magnitude and the relevant phase
spat_freq = fftshift(ifft2(ifftshift(test_freq))); % IFFT back into the image domain

% Find the range of points to be plotted
xminmax = [min(vesloc(1,:)) max(vesloc(1,:))];
xrange = xminmax(2)-xminmax(1);
xIdx = (xy>=xminmax(1)-xrange/2) & (xy<=xminmax(2)+xrange/2);
xred = xy(xIdx);
yminmax = [min(vesloc(2,:)) max(vesloc(2,:))];
yrange = yminmax(2)-yminmax(1);
yIdx = (xy>=yminmax(1)-yrange/2) & (xy<=yminmax(2)+yrange/2);
yred = xy(yIdx);

% Cut down the spatial frequency representation to save space
spat_freq = spat_freq(yIdx,xIdx);

% If off-resonance is present, interpolate the off-resonance values to make
% a map for visualisation
if ~isempty(off_res)
    [x,y] = meshgrid(xred,yred);
    
    % Fit a biharmonic surface (matches perfectly at vessel locations and varies smoothly in between)
    fit_surface = fit([vesloc(1,:)', vesloc(2,:)'], off_res(:), 'biharmonic'); % 2D biharmonic fit

    % Evaluate at all other points
    off_res_map = fit_surface(x, y);

    % Add this phase to the spatial frequency phase contribution
    spat_freq = spat_freq .* exp(1i*off_res_map);
end

% Normalise
maximum = max(max(real(spat_freq)));
spat_freq = real(spat_freq)/maximum;

% Plot
imagesc(xred,yred,spat_freq)
tmp = autumn(512); tmp = tmp(round(end/5):round(4*end/5),:); colormap(tmp)
caxis([-1 1]);
axis xy; axis equal; axis tight
xlabel('x (mm)','fontsize',14);
ylabel('y (mm)','fontsize',14);
set(gca,'fontsize',14);
title('Hadamard encoded vessels + optimised spatial frequency (fft)')
colorbar('fontsize',14);

% Plot vessels
TagIdx = find(ideal_enc_mtx_row == -1);
CntlIdx = find(ideal_enc_mtx_row == 1);
hold on;
plot(vesloc(1,TagIdx), vesloc(2,TagIdx), 'ko','MarkerFaceColor','k','MarkerSize',10)
plot(vesloc(1,CntlIdx),vesloc(2,CntlIdx),'ko','MarkerFaceColor','w','MarkerSize',10)
