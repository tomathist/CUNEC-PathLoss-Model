%% demo_PL_CUNEC_0th.m
% CUNEC Zeroth-Order Path Loss with Spatially Correlated Log-Normal Shadowing
% This script computes and visualizes path loss using the CUNEC zeroth-order model
% with shadowing correlated jointly across UEs and APs.

clear; clc;
rng(123);                             % Reproducibility (optional)

%% Geometry (3D positions in meters)
M_UE  = 2;                          % # UEs
N_AP  = 200;                          % # APs

% UE at (x,y,z) = (0,-10,15)
p_UE_0 = [0, -10, 15;
    0, -10, 15];

% APs along y-axis at z=1.5 m, spaced 2 m in y
p_AP_0 = [zeros(N_AP,1), (2*(1:N_AP)).', 1.5*ones(N_AP,1)];

%% Environment (meters)
building_len = 50;                    % building length b
street_w     = 15;                    % street width w
building_h   = 15;                    % building height h

%% Model params & correlations (user-provided loaders)
% Expects: FSPL_1m_3pt5GHz  (scalar, dB)
%          mu_0, sigma_0    (9x1 vectors)
%          C_0              (9x9 correlation matrix)
run('load_model_parameters.m');
run('load_correlations.m');

%% Monte Carlo realizations
R = 2;                                % keep small; each draw adds a full correlated field

%% Compute path loss(s)
[D, PL_CUNEC_zeroth] = calc_PL_CUNEC_0th(p_AP_0, p_UE_0, ...
    building_len, street_w, building_h, ...
    FSPL_1m_3pt5GHz, mu_0, sigma_0, C_0, R);

%% Plot
figure('Color','w'); hold on; grid on;
plot(D(1,:), squeeze(PL_CUNEC_zeroth(1,1,:)), 'o-', 'LineWidth', 1.8);
xlabel('Distance from Transmitter (m)'); ylabel('Path Loss (dB)');
title('CUNEC Zeroth-Order Path Loss with Correlated Shadowing');
legend('Zeroth-Order (1 realization)','Location','northwest'); box on;

% Heatmap of the shadowing-added PL for the single UE over all APs (optional)
if R==1
    figure('Color','w');
    imagesc(squeeze(PL_CUNEC_zeroth(1,1,:)).');
    axis tight; colorbar; box on;
    xlabel('AP index'); ylabel('UE index (fixed at 1)'); title('PL (dB)');
end

