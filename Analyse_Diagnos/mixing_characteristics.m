%MIXING_CHARACTERISTICS - Calculates some characteristics which are
%relevant to mixing 
% Sam Hartharn-Evans, 2022

clc; clearvars; close all; 
spinsstartup;

[x, z] = spinsgrid2d;
params = spins_params;

if isfield('params', 'hill_height')
    % When using a mapped version with a hill, we want something away from
    % the hill and from the wave generation zone
    undisturbed_tank = ((params.Lx - (params.hill_height/params.hill_slope)) + params.L_adj)/2;
    tank_ind = nearest_index(x(:, 1), undisturbed_tank);
else
    tank_ind = params.Nx/2;
end

%% Locate the pycnocline & vertical grid spacing
pyc_top = params.pyc_loc + params.h_halfwidth;
pyc_bot = params.pyc_loc - params.h_halfwidth;

pyc_top_ind = nearest_index(z(tank_ind, :), pyc_top);
pyc_bot_ind = nearest_index(z(tank_ind, :), pyc_bot);
n_pyc_zcoords = pyc_top_ind - pyc_bot_ind;
d_pyc_zcoords = (params.h_halfwidth*2)/n_pyc_zcoords;

dx = params.Lx / params.Nx;
n_pyc_xcoords = params.h_halfwidth*2 / dx;

%% Print the above
fprintf('Pycnocline Thickness: %2.2f m \n', 2*params.h_halfwidth);
fprintf('Nz in pycnocline: %i2 \n', n_pyc_zcoords);
fprintf('Average Dz in pycnocline: %4.4f m \n', d_pyc_zcoords);
fprintf('L_pyc / Dx: %4.4f \n', n_pyc_xcoords);

%% Schmidt Number
Sc = params.visco / params.kappa_rho;

fprintf('------------------------\n');
fprintf('Schmidt Number: %4.4f m \n', Sc);

%% Characterise Filter
nyquist_freq = 2*pi/params.Lx;
ks = [0:params.Nx/2-1 0 -params.Nx/2+1:-1]*nyquist_freq;

aliasing_wavelength = 2/params.Nx;
fprintf('------------------------\n');
fprintf('Aliasing Wavelength: %4.4f m \n', aliasing_wavelength);

% Transisions with wavelengths smaller thant he Nyquist wavelength cannot
% be resolved
filtalpha = params.f_strength; 
filtbeta = params.f_order;
f_cutoff = params.f_cutoff; 
knyq = max(abs(ks));
k_ny = (pi/params.Lx)*(params.Nx-1);

dummy = ones(size(ks));
if filtalpha>0
    kcut = f_cutoff*knyq;
    myfiltu = dummy.*(abs(ks)<kcut)+exp(-filtalpha*(((abs(ks)-kcut)/(knyq-kcut)).^filtbeta)).*(abs(ks)>=kcut);
    title('Exponential Filter')
    xlabel('Wavenumber')
    ylabel('Filter Strength')
    xlim([-knyq knyq])
else
    f_cutoff = 0;
    kcut = f_cutoff*knyq;
    myfiltu = dummy.*(abs(ks)<kcut)+exp(-filtalpha*(((abs(ks)-kcut)/(knyq-kcut)).^filtbeta)).*(abs(ks)>=kcut);    
    title('Indicative Hyperviscosity')
end

plot(ks, myfiltu, 'b.')
print('filter_representation.png', '-dpng');