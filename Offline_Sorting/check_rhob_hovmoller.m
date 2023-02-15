%CHECK_RHOB_HOVMOLLER - Visual inspector for weird things happening in the sorting algorithm. 
% allows the user to check a hovmoller plot of the background density plot, change in total BPE
% change in mass

clearvars; close all; clc; 
t1 = first_output(); t2 = last_output();
Nt = length(t1:t2);
params = spins_params;

% Pre-allocate arrays 
rhob = NaN(Nt, params.Nz);
mass = NaN(Nt, 1); BPE = mass;

% Choice of filter to apply (uncomment the one you want)
%filtname = 'nofilt';
%filtname = 'remfilt';
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
    rhob(ii+1, :) = rhob_temp(1, :); % just take the far tank value (it is the same everywhere, and then truncated by the slope)
end

%% Now plot
close all;
z = zgrid_reader;
% Plot the BPE hovmoller
figure; 
tiledlayout(4, 1);
nexttile;
pcolor(t1:t2, z(1, :), rhob'); shading flat; c = colorbar; ylabel(c, 'rhob');
axis tight; box on;
xlabel('t'); ylabel('z (m)');
  
% Plot BPE-BPE(1) hovmoller
nexttile;
pcolor(t1:t2, z(1, :), (rhob-rhob(1, :))'); shading flat; c = colorbar; ylabel(c, 'rhob-rhob(1)');
axis tight; newbluewhitered; box on;
xlabel('t'); ylabel('z (m)');

% Plot line plot of mass deviation
nexttile;
plot(t1:t2, mass/mass(1)*100);
ylabel('mass(t)/mass(0) (%)')
xlabel('t');

% Plot line plot of BPE deviation
nexttile;%subplot(4, 1, 4)
plot(t1:t2, BPE/BPE(1)*100);
ylabel('BPE(t)/BPE(1) (%)');
xlabel('t');

print('rhob_plot.png', '-dpng');
