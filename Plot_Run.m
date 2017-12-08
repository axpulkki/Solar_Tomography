%PLOT_RUN
%

%   Antti Pulkkinen, July 2017.

cc;

% Run to plot.
run_directory = './DataCube_swmf_CCMC_lowres_3_sc';

% Load the run.
load(sprintf('%s/run_data.mat',run_directory));

[n,m,s] = size(data);

eval(sprintf('mkdir %s/slices_frames',run_directory'));

for ssSlices = 1:s,
    
    figure;
    subplot(1,2,1), pcolor(squeeze(x_data(:,:,ssSlices))/Rs,squeeze(y_data(:,:,ssSlices))/Rs,(squeeze(data(:,:,ssSlices))));
    colorbar(gca,'SouthOutside'); title(sprintf('True z: %01.1f Rs',z_data(1,1,ssSlices)/Rs)); axis equal;
    caxis([min(min(squeeze(data(:,:,ssSlices)))) max(max(squeeze(data(:,:,ssSlices))))]); shading('interp');
    xlabel('x [Rs]'); ylabel('y [Rs]');
    
    subplot(1,2,2), pcolor(squeeze(x_data(:,:,ssSlices))/Rs,squeeze(y_data(:,:,ssSlices))/Rs,(squeeze(Ne_inverted(:,:,ssSlices))));
    colorbar(gca,'SouthOutside'); title(sprintf('Reconstructed z: %01.1f Rs',z_data(1,1,ssSlices)/Rs)); axis equal;
    caxis([min(min(squeeze(data(:,:,ssSlices)))) max(max(squeeze(data(:,:,ssSlices))))]); shading('interp');
    xlabel('x [Rs]'); ylabel('y [Rs]');
    
    eval(sprintf('print -dpng %s/slices_frames/frame_%01.0f.png',run_directory,ssSlices)); close;
    
end;
