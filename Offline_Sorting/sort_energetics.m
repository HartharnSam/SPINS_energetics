function [energy_structure, rhob, tot_mass] = sort_energetics(ii, xlimits, filter_rho, filt_remove)
% Script which employs the Winters (1995)-like sorting algorithm for a
% mapped case to calculate APE & BPE, and then also calculates KE and
% Packages nicely into a structure
%
% Syntax:  [energy_structure, rhob, tot_mass] = sort_energetics(ii, xlimits, filter_rho, filt_remove)
%
% Inputs:
%    ii - Time to calculate energetics for
%    xlimits - x region over which to calculate the sorting algorithm for
%    filter_rho - boolean switch to replace anomalously high/low densities
%    with the max/min expected densities respectively
%    filt_remove - boolean switch to instead remove those values and
%    replace with NaN (not advised)
%
% Outputs:
%    energy_structure - Output structure containing each energy component
%    rhob - background density profile rhob(x, z)
%    tot_mass - current total mass of model after filtering
%
% Other m-files required: cheb, clencurt, nearest_index, spins_params,
% spins_reader_new, spinsgrid2d
%
% Author: Marek Stastna / Sam Hartharn-Evans
% GitHub: https://github.com/HartharnSam
% Dec-2021; Last revision: 12-Jul-2022
% MATLAB Version: 9.12.0.1975300 (R2022a) Update 3
%
%--------------------------
%------START SCRIPT -------
%--------------------------
% Read in and cut down the grid
params = spins_params;
[x, z] = spinsgrid2d;
if nargin > 1
    xind1 = nearest_index(x(:, 1), xlimits(1));
    xind2 = nearest_index(x(:, 1), xlimits(2));
else
    xind1 = 1;
    xind2 = length(x);
end

xinds = xind1:xind2;
x = x(xinds, :); z = z(xinds, :);

if nargin < 3
    filter_rho = false;
end

% Read in spins.conf - and calculate the related infos
sz = size(x);
Nx = sz(1);
Nz = sz(2);
Nzc = Nz - 1; % Nz for chebyshev
dx = x(2) - x(1);
z = z - params.min_z;

% Read in rho, u and w data
rho_0 = params.rho_0;
rho = spins_reader_new('rho', ii, xinds, []);
rho = (rho*rho_0 + rho_0); % Confirmed - this is what compute_BPE does

u = spins_reader_new('u',ii, xinds, []);
w = spins_reader_new('w',ii, xinds, []);

% Filter out erroneously high density (that results from model
% mixing/filter)
if nargin < 4
    filt_remove = false;
end
if filter_rho
    % Read in a reference density profile
    rho0 = (spins_reader_new('rho',0))*rho_0 + rho_0;
    rhomax = max(rho0(:));
    rhomin = min(rho0(:));

    if ~filt_remove
        % Tries a version where we replace any overshoots with the closest
        % relevant density
        %  Produces a "myind" - matrix of where density overshoots
        myind = 1.0*(rho>rhomax);%+(myfact*delrho));
        rho = rhomax*myind + rho.*(1-myind);

        myind = 1.0*(rho<rhomin);%-(myfact*delrho));
        rho = rhomin*myind + rho.*(1-myind);

    else
        % Tries a version where we remove any overshoots
        myind = 1.0*(rho>rhomax | rho<rhomin);
        rho = rho.*(1-myind);
        u = u.*(1-myind); % The purpose of doing this is just to treat KE, PE etc in the same way as BPE
        w = w.*(1-myind);
    end
end

ke = 0.5*rho_0.*(u.^2+w.^2); % calculate Kinetic Energy

%% Calculate Chebyschev volumes
% Compute the area associated with each Chebyshev point using the values
% halfway between the point below and above
[~,z1dc] = cheb(Nzc);
[~,wci] = clencurt(Nzc);

% A normalised height applied to each Chebyshev point
arc(1) = 0.5*(z1dc(1)-z1dc(2)); % .5 comes from clencurt returning a total weight of 2
arc(Nzc+1) = arc(1);
for jj=2:(Nzc) % Central differencing
    arc(jj) = 0.5*(z1dc(jj-1)-z1dc(jj))+0.5*(z1dc(jj)-z1dc(jj+1)); % d_z1dc
end

% Then normalise this by the change in depth at each x position
arcphys = NaN(Nx, Nz);
for jj = 1:Nx
    Lznow = max(z(jj, :))-min(z(jj, :));
    arcphys(jj, :) = arc*Lznow; % Thickness of each grid point
end

if filt_remove
    arcphys(logical(myind)) = NaN; % Remove 
end
arcphysv = arcphys(:);

%% Now do the sorting algorithm
% For chebyshev it takes a bit of work

[rhosortedc, rhosortedci] = sort(rho(:), 'descend');

% Now create the zsortedc
zsortedc = NaN(1, length(arcphysv));
zsortedc(1) = arcphysv(rhosortedci(1));
for jj = 2:(length(arcphysv))  % For each grid point
    zsortedc(jj) = zsortedc(jj-1)+arcphysv(rhosortedci(jj)); % point weights sorted by density
end
zsortedc = zsortedc-min(zsortedc);
zsortedc = zsortedc/max(zsortedc);

% This assumes the maximum depth is at the left end point
% Create a physical version of the sorted grid
zmin = min(z(1,:));
zmax = max(z(1,:));
zsortedc = zmin+(zmax-zmin)*zsortedc;

if filt_remove
    zsortedc = zsortedc(~logical(myind(rhosortedci)));
    rhosortedc = rhosortedc(~logical(myind(rhosortedci)));
end
% now create a working grid so the repeated interpolations don't take
% forever - decrease the resolutions of zsortedc/rhosortedc
zsortedworking = linspace(zmin,zmax,256);
rhosortedworking = interp1(zsortedc, rhosortedc, zsortedworking,'linear');

% Finally create the background density profile at each x value
% And get the total APE as well as the tot APE at each x value
% I compute the KE info as well
energy_structure = 0;
ke_tot = 0;
bpe_tot = 0;
rhob = NaN(Nx, Nz);
ape_hor = NaN(1, Nx);
pe_tot = 0;
pe_hor = NaN(Nx, Nz);
ke_hor = NaN(1, Nx);
bpe_hor = NaN(1, Nx);
wi_tot = NaN(Nz, Nx);
zi = z;

for jj = 1:Nx
    zminnow = min(z(jj,:));
    zmaxnow = max(z(jj,:));
    rhob(jj,:) = interp1(zsortedworking, rhosortedworking, zi(jj,:), 'spline'); % background (sorted) stratification

    % get the local chain rule expression - Weight of each grid point
    winow = wci*(zmaxnow-zminnow)*0.5*dx;

    % Get APE
    % integrate vertically
    dummy = 9.81*sum(winow.*(rho(jj,:) - rhob(jj,:)).*zi(jj,:));
    % notice that if the surface is at z=0 z will have negative values
    ape_hor(jj) = dummy;
    energy_structure = energy_structure+dummy;

    % Get BPE
    % integrate vertically
    dummy = 9.81*sum(winow.*(rhob(jj,:)).*z(jj,:));
    % notice that if the surface is at z=0 z will have negative values
    bpe_hor(jj) = dummy;
    bpe_tot = bpe_tot+dummy;

    % Get PE
    dummy = 9.81*sum(winow.*rho(jj, :).*z(jj,:));
    pe_hor(jj) = dummy;
    pe_tot = pe_tot+dummy;

    % Get KE
    % integrate vertically
    dummy = sum(winow.*ke(jj,:));
    ke_hor(jj) = dummy;
    ke_tot = ke_tot+dummy;
    wi_tot(:, jj) = winow;
end

energy_structure = struct('APE_Total', energy_structure, 'PE_Total', pe_tot, 'KE_Total', ke_tot, ...
    'BPE_Total', bpe_tot);

tot_mass = sum(wi_tot(:).*rho(:));

end
