%PLOT_RMS_SPACECRAFT_NO

cc;

spacecraft_no = 2:8;

for ii = spacecraft_no,
    
    load(sprintf('./DataCube_Spherical_Gaussian_run_%01.0f/run_data.mat',ii));
    
    rms_spacecraft(ii) = rms_difference(end);
    
end;

figure; plot(spacecraft_no,rms_spacecraft(spacecraft_no)); hold on; plot(spacecraft_no,rms_spacecraft(spacecraft_no),'.','color','k','markersize',10);
grid on; xlabel('Number of spacecraft'); ylabel('RMS [#/m^3]');
title('RMS difference as a function of number of spacecraft');
