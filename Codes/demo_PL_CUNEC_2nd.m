%% demo_PL_CUNEC_2nd.m
% CUNEC Zeroth + First + Second-Order Path Loss with Jointly Correlated Shadowing
% - Uses a single AP grid for all orders so sums align in [R x N_UE x N_AP]
% - Two UEs; APs extend south from the 2nd corner along x = 70 m
% - Plots total PL vs. AP distance from the second corner with MC mean ± 1 SD

clear; clc;
rng(123);  % reproducibility

%% ---------------- Geometry (3D, meters) ----------------
N_AP      = 100;

% UEs (two examples)
p_UE      = [0,  0, 1.5;   ...   % UE #1
             0, 50, 1.5];        % UE #2
N_UE      = size(p_UE,1);

% Corners
p_corner_1 = [ 0, 70, 1.5];      % first corner
p_corner_2 = [70, 70, 1.5];      % second corner

% APs: along negative y from the 2nd corner at x = 70 m (southbound “leg”)
p_AP       = zeros(N_AP,3);
p_AP(:,1)  = 70;                          % x
p_AP(:,2)  = 70 - (1:N_AP)';              % y: 69, 68, ..., (south)
p_AP(:,3)  = 1.5;                         % z (AP height)

%% ---------------- Environment (meters) ----------------
building_len = 20;      % b
street_w     = 20;      % w (used by 0th order)
building_h   = 15;      % h
near_far     = "far";   % first-order offset gate enabled

%% ---------------- Model parameters & correlations ----------------
% Expected from the loaders:
%   FSPL_1m_3pt5GHz (scalar, dB)
%   mu_0, sigma_0, C_0   (0th: 9x1, 9x1,  9x9)
%   mu_1, sigma_1, C_1   (1st: 10x1,10x1,10x10)
%   mu_2, sigma_2, C_2   (2nd: 6x1,  6x1,  6x6)
run('load_model_parameters.m');
run('load_correlations.m');

%% ---------------- Monte Carlo settings ----------------
R = 20;  % keep modest: each draw generates a correlated UE–AP field

%% ---------------- Compute each order on the SAME AP/UE grid ----------------
% Zeroth order
[D0, PL0] = calc_PL_CUNEC_0th( ...
    p_corner_1, p_UE, ...
    building_len, street_w, building_h, ...
    FSPL_1m_3pt5GHz, mu_0, sigma_0, C_0, R);

% First order (corner = p_corner_1)
PL1 = calc_PL_CUNEC_1st( ...
    p_corner_2, p_UE, p_corner_1, near_far, ...
    building_len, building_h, ...
    mu_1, sigma_1, C_1, R);

% Second order (corner = p_corner_2)
PL2 = calc_PL_CUNEC_2nd( ...
    p_AP, p_UE, p_corner_2, ...
    building_len, mu_2, sigma_2, C_2, R);

%% ---------------- Total path loss ----------------
PL_total = PL0 + PL1 + PL2;   % [R x N_UE x N_AP]

%% ---------------- Plot: total PL vs. distance from the 2nd corner ----------------
% x-axis: AP distance to the second corner
d_AP2corner2 = sqrt(sum((p_AP - p_corner_2).^2, 2));   % [N_AP x 1]
[x_sorted, idx] = sort(d_AP2corner2, 'ascend');

% Style
set(0,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',11);
set(0,'DefaultLineLineWidth',1.6, 'DefaultFigureColor','w');

figure; hold on; grid on; box on;

% For each UE, show mean ± 1 SD band and mean curve
for u = 1:N_UE
    pl_u = squeeze(PL_total(:,u,:));     % [R x N_AP]
    m_u  = mean(pl_u, 1);                % [1 x N_AP]
    s_u  = std(pl_u, 0, 1);              % [1 x N_AP]

    m_u  = m_u(idx);
    s_u  = s_u(idx);

    upper = m_u + s_u;
    lower = m_u - s_u;

    % Shaded band
    h_fill = fill([x_sorted', fliplr(x_sorted')], [upper, fliplr(lower)], ...
                  [0.85 0.85 0.95], 'EdgeColor','none', 'FaceAlpha',0.55);
    % Mean line
    h_line = plot(x_sorted, m_u, 'k-', 'DisplayName', sprintf('UE %d (mean)', u));

    % Put the band behind the line
    uistack(h_fill, 'bottom');
end

xlabel('AP distance from second corner (m)');
ylabel('Total path loss (dB)');
title('CUNEC Total Path Loss = Zeroth + First + Second (Jointly Correlated Shadowing)');
legend('Location','northwest');

% Save figures
out_png = 'cunec_totalPL_2corners.png';
out_fig = 'cunec_totalPL_2corners.fig';
exportgraphics(gcf, out_png, 'Resolution', 300);

