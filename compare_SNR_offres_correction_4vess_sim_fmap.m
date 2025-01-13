% Test the OES off-resonance correction algorithm with four vessels using a
% simulated fieldmap

%% Set parameters
n_vessel_geometries = 20; 
motion = 4; 
pad_size = 1024;
vesloc3D = 0;
AngTowardsCor = 0;
AngTowardsSag = 0;
debug = 0;
vesloc = [-7 4 -6 7; 4 4 14 14];
shifts = 2*randn(2,size(vesloc,2),n_vessel_geometries);
shifts(:,:,1) = zeros(size(vesloc)); %set first set of shifts to zero
had_enc = 0;
set_iterations = 20;

%% Run a simulation for each vessel geometry
for v = 1:n_vessel_geometries
    
    disp(['Vessel geometry ' num2str(v) ' of ' num2str(n_vessel_geometries)]);
    vessels = round(vesloc+shifts(:,:,v));
    
    % PCASL encoding
    standard_enc = 0;
    PCASL_enc = 1;

    % simulated fieldmap, 4 vessels, 
    realfmap = 0;
    simfmap = 1;
    linear_offset = pi/3;
    linear_offset_dirn = 1;
    quadratic_offset = 1;

    [meanTotalSNR_simfmap_4PCASL(v,:), ~, ~, ~] = test_offset_fn_sim(linear_offset, linear_offset_dirn, quadratic_offset,...
        vessels,motion,pad_size,had_enc,standard_enc,PCASL_enc,vesloc3D,AngTowardsCor,AngTowardsSag,debug,set_iterations);

    % VEPCASL standard encoding
    standard_enc = 1;
    PCASL_enc = 0;

    % simulated fieldmap, 4 vessels
    [meanTotalSNR_simfmap_4VEPCASL(v,:), ~, ~, ~] = test_offset_fn_sim(linear_offset, linear_offset_dirn, quadratic_offset,...
        vessels,motion,pad_size,had_enc,standard_enc,PCASL_enc,vesloc3D,AngTowardsCor,AngTowardsSag,debug,set_iterations);

end

% Final mean, std, h and p values
SNR_simfmap_4PCASL(1,:) = mean(meanTotalSNR_simfmap_4PCASL);
SNR_simfmap_4PCASL(2,:) = std(meanTotalSNR_simfmap_4PCASL);
[~, p_simfmap_4PCASL(1)] = ttest(meanTotalSNR_simfmap_4PCASL(:,1), meanTotalSNR_simfmap_4PCASL(:,2));
[~, p_simfmap_4PCASL(2)] = ttest(meanTotalSNR_simfmap_4PCASL(:,1), meanTotalSNR_simfmap_4PCASL(:,3));
[~, p_simfmap_4PCASL(4)] = ttest(meanTotalSNR_simfmap_4PCASL(:,2), meanTotalSNR_simfmap_4PCASL(:,3));

SNR_simfmap_4VEPCASL(1,:) = mean(meanTotalSNR_simfmap_4VEPCASL);
SNR_simfmap_4VEPCASL(2,:) = std(meanTotalSNR_simfmap_4VEPCASL);
[~, p_simfmap_4VEPCASL(1)] = ttest(meanTotalSNR_simfmap_4VEPCASL(:,1), meanTotalSNR_simfmap_4VEPCASL(:,2));
[~, p_simfmap_4VEPCASL(2)] = ttest(meanTotalSNR_simfmap_4VEPCASL(:,1), meanTotalSNR_simfmap_4VEPCASL(:,3));
[~, p_simfmap_4VEPCASL(4)] = ttest(meanTotalSNR_simfmap_4VEPCASL(:,2), meanTotalSNR_simfmap_4VEPCASL(:,3));

% Note that the SNR values here account for the imperfect labelling efficiency, so will always be less than one,
% even for a perfect encoding.

%% Plot similar to Berry et al. NeuroIm 2019 Fig2c,e
figure; subplot(1,2,1);
% Squash a little
Pos = get(gcf,'Position');
set(gcf,'Position',[Pos(1:2) Pos(3)*1.2 Pos(4)]);
set(gca, 'color','none');
ylabel('SNR +/- S.D.');
set(gca,'xticklabel',{'','',''},'xtick',1:3);
set(gca,'fontsize',30,'FontName','Times New Roman');
axis([0 4 0 1])
hold on
bar(1,SNR_simfmap_4PCASL(1,1),0.9,'facecolor',[0.3 0.2 0.4], 'linestyle','--','linewidth',4);
bar(2,SNR_simfmap_4PCASL(1,2),0.9,'facecolor',[0.5 0.4 0.6], 'linestyle','-.','linewidth',4);
bar(3,SNR_simfmap_4PCASL(1,3),0.9,'facecolor',[0.7 0.6 0.8], 'linestyle',':','linewidth',4);
errorbar([SNR_simfmap_4PCASL(1,:)], [SNR_simfmap_4PCASL(2,:)],'.k','linewidth',5);
title 'PCASL'
legend({'No offset','Offset','Corrected'},'location','southoutside','orientation','horizontal');
hold off

subplot(1,2,2);
set(gca, 'color','none');
set(gca,'xticklabel',{'','',''},'xtick',1:3);
set(gca,'fontsize',30,'FontName','Times New Roman');
axis([0 4 0 1])
hold on
bar(1,SNR_simfmap_4VEPCASL(1,1),0.9,'facecolor',[0.3 0.2 0.4], 'linestyle','--','linewidth',4);
bar(2,SNR_simfmap_4VEPCASL(1,2),0.9,'facecolor',[0.5 0.4 0.6], 'linestyle','-.','linewidth',4);
bar(3,SNR_simfmap_4VEPCASL(1,3),0.9,'facecolor',[0.7 0.6 0.8], 'linestyle',':','linewidth',4);
errorbar([SNR_simfmap_4VEPCASL(1,:)], [SNR_simfmap_4VEPCASL(2,:)],'.k','linewidth',5);
title 'VEPCASL'
h = gca;
h.YAxis.Visible = 'off'; % remove y-axis

% Reposition legend
subplot(1,2,1); AxPos = get(gca,'Position');
h = legend; h.Position = h.Position + [0.4 -0.1 0 0];
