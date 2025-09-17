%% demo_PL_CUNEC_1st.m
% CUNEC First-Order Path Loss with Jointly Correlated Log-Normal Shadowing
% This script computes and visualizes the FIRST-ORDER (corner) path-loss component,
% adds it to the zeroth-order baseline, and plots the total path loss versus
% AP distance from the corner. Monte Carlo variability is shown as a shaded band.

clear; clc;
rng(123);  % Reproducibility

%% ---------------- Geometry (3D, meters) ----------------
% One UE near origin, APs along the negative x-axis at y=70 m
N_AP    = 100;
p_UE_1  = [0, 0, 1.5; 2 0 1.5];                     % [2 x 3]
p_AP_1  = zeros(N_AP,3);
p_AP_1(:,1) = - (1:N_AP)';                          % x = -1:-N_AP
p_AP_1(:,2) = 70;                                   % y = 70
p_AP_1(:,3) = 1.5;                                  % z = 1.5
p_corner_1  = [0, 70, 1.5];                         % corner aligned with AP row

% For consistency, use the same geometry for zeroth-order
p_UE_0 = p_UE_1;
p_AP_0 = p_AP_1;

%% ---------------- Environment (meters) ----------------
building_len = 20;      % b
street_w     = 20;      % w (latent for first-order; used by zeroth-order)
building_h   = 15;      % h
near_far     = "far";   % enables offset adjustment in the first-order term

%% ---------------- Model parameters & correlations ----------------
% Expected to define:
%   FSPL_1m_3pt5GHz (scalar, dB)
%   mu_0,  sigma_0,  C_0  (zeroth-order: 9x1, 9x1, 9x9)
%   mu_1,  sigma_1,  C_1  (first-order: 10x1, 10x1, 10x10)
run('load_model_parameters.m');
run('load_correlations.m');

%% ---------------- Monte Carlo settings ----------------
R = 10;  % Use >1 to show variability; keep modest due to correlated-field generation cost.

%% ---------------- Compute zeroth-order PL ----------------
% Returns:
%   D0: [N_UE x N_AP] UE–AP distances (m)
%   PL0: [R x N_UE x N_AP] path-loss realizations (dB)
[D0, PL0] = calc_PL_CUNEC_0th( ...
    p_corner_1, p_UE_0, ...
    building_len, street_w, building_h, ...
    FSPL_1m_3pt5GHz, mu_0, sigma_0, C_0, R);

%% ---------------- Compute first-order PL (corner term) ----------------
% Returns:
%   PL1: [R x N_UE x N_AP] first-order realizations (dB)
PL1 = calc_PL_CUNEC_1st( ...
    p_AP_1, p_UE_1, p_corner_1, near_far, ...
    building_len, building_h, ...
    mu_1, sigma_1, C_1, R);

assert(all(size(PL1) == [R, size(p_UE_1,1), size(p_AP_1,1)]), 'PL1 size mismatch.');

%% ---------------- Total path loss (zeroth + first) ----------------
PL_total = PL0 + PL1;  % [R x N_UE x N_AP]

% For this demo, we have a single UE (index 1). Extract the AP profile.
pl_total_AP = squeeze(PL_total(:,1,:));     % [R x N_AP]

% Monte Carlo statistics
pl_mean = mean(pl_total_AP, 1);             % [1 x N_AP]
pl_std  = std(pl_total_AP, 0, 1);           % [1 x N_AP]

% Distance from the corner for each AP (used as x-axis)
d_AP2corner = sqrt(sum((p_AP_1 - p_corner_1).^2, 2));  % [N_AP x 1]

%% ---------------- Plot (publication-ready) ----------------
% Style
set(0,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',11);
set(0,'DefaultLineLineWidth',1.6, 'DefaultFigureColor','w');

figure; hold on; grid on; box on;

% Shaded MC band: mean ± 1 SD
x = d_AP2corner(:)';                     % 1 x N_AP
[xx, idx] = sort(x, 'ascend');           % sort for a clean band
m  = pl_mean(idx);
s  = pl_std(idx);
upper = m + s;
lower = m - s;
fill([xx, fliplr(xx)], [upper, fliplr(lower)], [0.85 0.85 0.95], ...
     'EdgeColor','none', 'FaceAlpha',0.7, 'DisplayName','MC band (±1 SD)');

% Total mean curve
plot(xx, m, 'k-', 'DisplayName','Total PL (mean)');

xlabel('AP Distance from Corner (m)');
ylabel('Path Loss (dB)');
title('CUNEC Total Path Loss = Zeroth + First (Jointly Correlated Shadowing)');
xlim([min(xx) max(xx)]);
legend('Location','northwest');

% Save figures
out_png = 'cunec_firstorder_totalPL.png';
out_fig = 'cunec_firstorder_totalPL.fig';
exportgraphics(gcf, out_png, 'Resolution',300);
