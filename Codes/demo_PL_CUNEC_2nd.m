%% CUNEC First-Order Path Loss Model with Correlated Shadowing
% This script computes and visualizes path loss using the CUNEC first-order model
% with spatially correlated log-normal shadowing across multiple Monte Carlo realizations.

%% Input Parameters
p_UE_1          = [0 0 1.5];
for i = 1:100
    p_UE_1(i,3) = 1.5;
    p_UE_1(i,2) = -i;
end

p_AP_1          = zeros(100,3);
for i = 1:100
    p_AP_1(i,2) = 70;
    p_AP_1(i,1) = -i;
end

p_corner_1 = [0 70 1.5];
p_corner_2 = [210 20 1.5];
near_far = "far";

building        = 20;
width           = 10;
height          = 15;            % Fixed street width (latent variable)
N_realizations = 1000;          % Number of Monte Carlo realizations

%% Load Model Parameters
run('load_model_parameters.m'); % Loads mu_1, sigma_1, etc.
run('load_correlations.m');     % Loads C_1 (correlation matrix)

%% Compute First-Order Path Loss with Correlated Shadowing
[D, PL_CUNEC_zeroth] = calc_PL_CUNEC_0th( ...
    p_corner_1, ...
    p_UE_1, ...
    building, ...
    width, ...
    height, ...
    FSPL_1m_3pt5GHz, ...
    mu_0, ...
    sigma_0, ...
    C_0, ...
    N_realizations);

PL_CUNEC_first = calc_PL_CUNEC_1st( ...
    p_AP_1, ...
    p_UE_1, ...
    p_corner_1, ...
    near_far, ...
    building, ...
    height, ...
    mu_1, ...
    sigma_1, ...
    C_1, ...
    N_realizations);

mean_0 = mean(PL_CUNEC_zeroth);
mean_1 = squeeze(mean(PL_CUNEC_first,1));
mean_PL = mean_0+mean_1;

% Compute average path loss across all realizations
mean_PL_first  = squeeze(mean(PL_CUNEC_first, 1))+80;

mean_PL_CUNEC = mean(mean_PL_first(:,5:70));
pathloss_final = mean_PL_first(:,5:70);

%% Plot Results
figure; hold on; grid on;

% Plot mean path loss curve
plot(mean_PL_first(1,:), 'b-o', ...
    'LineWidth', 2, 'DisplayName', 'First-Order (Mean)');

xlim([5 70]);

% Plot formatting
xlabel('Distance from First Corner (m)', 'FontSize', 12);
ylabel('Added Path Loss From Corner (dB)', 'FontSize', 12);
title('CUNEC First-Order: Added Path Loss From Corner', 'FontSize', 14);
legend('Location', 'northwest');
set(gca, 'FontSize', 12);
