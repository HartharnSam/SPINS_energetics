clearvars; close all; clc; 
%cd('C:\Users\samha\OneDrive - Newcastle University\02_PhD_Project\06_Fission_Bolus\02_Raw_data\Thin_20L_33_220721')
t1 = 0; t2 = 150;
Nt = length(t1:t2);

rhob = NaN(Nt, 256);
mass = NaN(Nt, 1); BPE = mass;
filtname = 'nofilt';
filtname = 'remfilt';
filtname = 'filt';

switch filtname
    case 'nofilt'
        isFilt = false; isRemFilt = false;
    case 'remfilt'
        isFilt = true; isRemFilt = true;
    case 'filt'
        isFilt = true; isRemFilt = false;
end

for ii = t1:t2
    [energy, rhob_temp, mass(ii+1)] = sort_energetics(ii, [0 14.5], isFilt, isRemFilt);
    BPE(ii+1) = energy.BPE_Total;
    rhob(ii+1, :) = rhob_temp(1, :);
end

%%
close all; 
z = zgrid_reader;
cd('C:\Users\samha\OneDrive - Newcastle University\02_PhD_Project\07_CanadaMixing\02_Raw_data\HV_fission_order_08_strength_025')
figure; 
tiledlayout(4, 1);
nexttile;
%subplot(4, 1, 1);
pcolor(t1:t2, z(1, :), rhob'); shading flat; c = colorbar; ylabel(c, 'rhob');
axis tight; box on;
%print(['../../04_Output/01_Sort_Hill/rhob_hovmoller_', filtname, '.png'], '-dpng');

%
nexttile;
%subplot(4, 1, 2);
pcolor(t1:t2, z(1, :), (rhob-rhob(1, :))'); shading flat; c = colorbar; ylabel(c, 'rhob-rhob(1)');
axis tight; newbluewhitered; box on;
%print(['../../04_Output/01_Sort_Hill/rhob_anom_hovmoller_', filtname, '.png'], '-dpng');

%
nexttile;
%subplot(4, 1, 3)
plot(t1:t2, mass/mass(1)*100);
ylabel('mass(t)/mass(0) (%)')
nexttile;%subplot(4, 1, 4)
plot(t1:t2, BPE/BPE(1)*100);
ylabel('BPE(t)/BPE(1) (%)');

print(['../../04_Output/01_Sort_Hill/mass_dev_', filtname, '.png'], '-dpng');
