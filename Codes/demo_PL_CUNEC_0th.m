%% CUNEC Zeroth-Order Path Loss Model with Correlated Shadowing
% This script computes and visualizes path loss using the CUNEC zeroth-order model
% with spatially correlated log-normal shadowing across multiple Monte Carlo realizations.

%% Input Parameters
p_UE_0          = [0, -10, 15;
                   0, -40, 15];
p_AP_0          = zeros(200,3);
for i = 1:300
    p_AP_0(i,3) = 1.5;
    p_AP_0(i,2) = -10+i/2;
end

% Uncomment if you want to look at UE Trajectory
temp            = p_AP_0;
p_AP_0          = p_UE_0;
p_UE_0          = temp;

building        = 50;
width           = 15;           % Fixed street width (latent variable)
height          = 15;           % Fixed building height (latent variable)

% building        = 50;
% width           = 25;           % Fixed street width (latent variable)
% height          = 25;           % Fixed building height (latent variable)

N_realizations  = 100;         % Number of Monte Carlo realizations

%% Load Model Parameters
run('load_model_parameters.m'); % Loads FSPL_1m_3pt5GHz, mu_0, sigma_0, etc.
run('load_correlations.m');     % Loads C_0 (correlation matrix)

%% Compute Zeroth-Order Path Loss with Correlated Shadowing
[D, PL_CUNEC_zeroth] = calc_PL_CUNEC_0th( ...
    p_AP_0, ...
    p_UE_0, ...
    building, ...
    width, ...
    height, ...
    FSPL_1m_3pt5GHz, ...
    mu_0, ...
    sigma_0, ...
    C_0, ...
    N_realizations);

% Uncomment if you want to look at UE trajectory
PL_CUNEC_zeroth = permute(PL_CUNEC_zeroth, [1, 3, 2]);
D = D';

PL_CUNEC_zeroth_1 = squeeze(PL_CUNEC_zeroth(:,1,:));

% Compute average path loss across all realizations
mean_PL_zeroth  = squeeze(mean(PL_CUNEC_zeroth, 1));

%% Interpolate to Nearest Integer Distances
target_distances = 30:1:120; % integers from 30 to 120 meters

mean_PL_interp = zeros(size(target_distances));
for idx = 1:length(target_distances)
    [~, nearest_idx] = min(abs(D(1,:) - target_distances(idx))); % find closest
    mean_PL_interp(idx) = mean_PL_zeroth(1, nearest_idx);
    pathloss_final(:,idx) = PL_CUNEC_zeroth_1(:,nearest_idx);
end

mean_PL_CUNEC = mean_PL_interp;

%% Plot Results
figure; hold on; grid on;

% Plot mean path loss curve
plot(target_distances, mean_PL_interp, 'b-o', ...
    'LineWidth', 2, 'DisplayName', 'Zeroth-Order (Mean)');

% Plot formatting
xlabel('Distance from Transmitter (m)', 'FontSize', 12);
ylabel('Path Loss (dB)', 'FontSize', 12);
title('CUNEC Zeroth-Order: Path Loss', 'FontSize', 14);
legend('Location', 'northwest');
set(gca, 'FontSize', 12);

