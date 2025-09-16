function [D, PL_CUNEC_zeroth] = calc_PL_CUNEC_0th(p_AP_0, p_UE_0, b, w, h, FSPL_1m, mu, sigma, C, N)
% calc_PL_CUNEC_0th
% Computes stochastic path loss values using the zeroth-order CUNEC model.
% This model accounts for environmental variability and parameter correlation
% using conditional multivariate normal sampling conditioned on fixed building and street parameters.
%
% Inputs:
%   p_AP_0          - [N_AP × 3] matrix of AP 3D positions (meters)
%   p_UE_0          - [N_UE × 3] matrix of UE 3D positions (meters)
%   b               - Building length (meters, scalar)
%   w               - Street width (meters, scalar)
%   h               - Building height (meters, scalar)
%   FSPL_1m         - Free-space path loss at 1 meter, in dB (scalar)
%   mu              - Mean vector of 9 parameters
%   sigma           - Standard deviation vector (length 9)
%   C               - 9×9 correlation matrix
%   N               - Number of random path loss realizations to generate
%
% Output:
%   PL_CUNEC_zeroth - [N × N_UE × N_AP] matrix of simulated path loss values in dB

    %% Input Validation
    if size(p_AP_0,2) ~= 3
        error('"p_AP_0" must be N_AP x 3 matrix of 3D positions.');
    end
    if size(p_UE_0,2) ~= 3
        error('"p_UE_0" must be N_UE x 3 matrix of 3D positions.');
    end
    if ~isscalar(b) || b <= 0
        error('"b" (building length) must be a positive scalar.');
    end
    if ~isscalar(w) || w <= 0
        error('"w" (street width) must be a positive scalar.');
    end
    if ~isscalar(h) || h <= 0
        error('"h" (height) must be a positive scalar.');
    end
    if ~isscalar(FSPL_1m) || FSPL_1m <= 0
        error('"FSPL_1m" must be a positive scalar.');
    end
    if ~isvector(mu) || length(mu) ~= 9
        error('"mu" must be a vector of length 9.');
    end
    if ~isvector(sigma) || length(sigma) ~= 9
        error('"sigma" must be a vector of length 9.');
    end
    if any(sigma <= 0)
        error('All elements of "sigma" must be strictly positive.');
    end
    if ~ismatrix(C) || ~isequal(size(C), [9, 9]) || ~issymmetric(C)
        error('"C" must be a symmetric 9×9 matrix.');
    end
    if ~isscalar(N) || N <= 0 || mod(N, 1) ~= 0
        error('"N" must be a positive integer.');
    end

    %% Step 1: Full Covariance Matrix
    Sigma                        = diag(sigma) * C * diag(sigma);

    %% Step 2: Condition on Fixed Width and Height
    x_given                      = [b, w, h];                          % Fixed values
    idx_fixed                    = [1, 2, 3];                          % Indices for width, height
    idx_free                     = setdiff(1:9, idx_fixed);            % Indices to sample

    mu_fixed                     = mu(idx_fixed);
    mu_free                      = mu(idx_free);

    Sigma_ff                     = Sigma(idx_fixed, idx_fixed);
    Sigma_uu                     = Sigma(idx_free, idx_free);
    Sigma_uf                     = Sigma(idx_free, idx_fixed);
    Sigma_fu                     = Sigma(idx_fixed, idx_free);

    mu_cond                      = mu_free + Sigma_uf / Sigma_ff * (x_given(:) - mu_fixed(:));
    Sigma_cond                   = Sigma_uu - Sigma_uf / Sigma_ff * Sigma_fu;
    Sigma_cond                   = nearestSPD(Sigma_cond);             % Ensure SPD

    %% Step 3: Sample from Conditional Distribution
    % Sampling: [a_0, Delta_0, sigma_S_0_AP, sigma_S_0_UE, corr_dist_0_AP, corr_dist_0_UE]
    X                            = mvnrnd(mu_cond', Sigma_cond, N);    % N x 6

    % Extract common parameters
    a_0_realizations             = max(X(:, 1), 0.01);
    Delta_0_realizations         = X(:, 2);
    sigma_0_AP_realizations      = max(X(:, 3), 0);
    sigma_0_UE_realizations      = max(X(:, 4), 0);
    d_corr_0_AP_realizations     = max(X(:, 5), 1);
    d_corr_0_UE_realizations     = max(X(:, 6), 1);

    %% Step 4: Compute Path Loss
    % p_UE_0: M x 3
    % p_AP_0: N x 3
    
    % Expand and compute pairwise distance
    D                            = sqrt( ...
                                       sum(p_AP_0.^2, 2).' + ...       % (1 x N) AP norms
                                       sum(p_UE_0.^2, 2) - ...         % (M x 1) UE norms
                                       2 * (p_UE_0 * p_AP_0.') ...     % (M x N) inner products
                                   );

    log_D0                       = 10 * log10(D);                      % 1 × D
    PL_mean                      = zeros(N, size(D,1), size(D,2));     % [N × N_UE × N_AP]
    for i = 1:N
        PL_mean(i, :, :)         = FSPL_1m + Delta_0_realizations(i) + a_0_realizations(i) * log_D0;
    end
    %% Step 5: Add Correlated Shadowing
    if ndims(PL_mean) == 3
        % Multiple realizations
        [N, U, A] = size(PL_mean);
        PL_CUNEC_zeroth = zeros(N, U, A);
        for i = 1:N
%             shadowing = generate_correlated_shadowing(p_AP_0, p_UE_0, ...
%                 sigma_0_UE_realizations(i), sigma_0_AP_realizations(i), ...
%                 d_corr_0_UE_realizations(i), d_corr_0_AP_realizations(i));
            shadowing = generate_correlated_shadowing(p_AP_0, p_UE_0, ...
               sigma_0_AP_realizations(i), sigma_0_UE_realizations(i), ...
               d_corr_0_AP_realizations(i), d_corr_0_UE_realizations(i));
            PL_CUNEC_zeroth(i, :, :) = squeeze(PL_mean(i, :, :)) + shadowing;
        end
    else
        % Single realization
        [U, A] = size(PL_mean);
        PL_CUNEC_zeroth = zeros(U, A);
        for i = 1:N
%             shadowing = generate_correlated_shadowing(p_AP_0, p_UE_0, ...
%                 sigma_0_UE_realizations(i), sigma_0_AP_realizations(i), ...
%                 d_corr_0_UE_realizations(i), d_corr_0_AP_realizations(i));
            shadowing = generate_correlated_shadowing(p_AP_0, p_UE_0, ...
                sigma_0_AP_realizations(i), sigma_0_UE_realizations(i), ...
                d_corr_0_AP_realizations(i), d_corr_0_UE_realizations(i));
            PL_CUNEC_zeroth(i,:) = PL_mean(i, :) + shadowing;
        end
    end
end
