function PL_CUNEC_first = calc_PL_CUNEC_1st(p_AP_1, p_UE_1, p_corner_1, near_far, b, h, mu, sigma, C, N)
% calc_PL_CUNEC_1st
% Computes stochastic path loss values using the first-order CUNEC model.
% This model incorporates geometry and environment-specific parameters and samples
% from a multivariate normal distribution to simulate path loss behavior.
%
% Inputs:
%   p_AP_1       - [N_AP × 3] matrix of AP 3D positions (meters)
%   p_UE_1       - [N_UE × 3] matrix of UE 3D positions (meters)
%   p_corner_1   - [1 × 3] 3D position of the relevant corner (meters)
%   b            - Building length (meters, scalar)
%   h            - Building height (meters, scalar)
%   mu           - Mean vector of 9 parameters
%   sigma        - Standard deviation vector (length 9)
%   C            - 9×9 correlation matrix
%   N            - Number of random path loss realizations to generate
%
% Outputs:
%   D               - [N_UE × N_AP] distance matrix (placeholder)
%   PL_CUNEC_first  - [N × N_UE × N_AP] matrix of path loss values (dB)

    %% Input Validation
    if size(p_AP_1,2) ~= 3
        error('"p_AP_1" must be N_AP × 3 matrix of 3D positions.');
    end
    if size(p_UE_1,2) ~= 3
        error('"p_UE_1" must be N_UE × 3 matrix of 3D positions.');
    end
    if ~isequal(size(p_corner_1), [1, 3])
        error('"p_corner_1" must be a 1 × 3 vector of corner position.');
    end
    if ~isscalar(b) || b <= 0
        error('"b" (building length) must be a positive scalar.');
    end
    if ~isscalar(h) || h <= 0
        error('"h" (building height) must be a positive scalar.');
    end
    if ~isvector(mu) || length(mu) ~= 10
        error('"mu" must be a vector of length 10.');
    end
    if ~isvector(sigma) || length(sigma) ~= 10 || any(sigma <= 0)
        error('"sigma" must be a vector of 10 strictly positive values.');
    end
    if ~isequal(size(C), [10, 10]) || ~issymmetric(C)
        error('"C" must be a symmetric 10 × 10 matrix.');
    end
    if ~isscalar(N) || N <= 0 || mod(N,1) ~= 0
        error('"N" must be a positive integer.');
    end

    %% Step 1: Full Covariance Matrix
    Sigma                        = diag(sigma) * C * diag(sigma);

    %% Step 2: Condition on Fixed Width and Height
    x_given                      = [b, h];                          % Fixed values
    idx_fixed                    = [1, 2];                          % Indices for building, height
    idx_free                     = setdiff(1:10, idx_fixed);         % Indices to sample

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
    % Sampling: [a_1, Delta_1, kappa_1, offset_1, sigma_S_1_AP, sigma_S_1_UE, corr_dist_1_AP, corr_dist_1_UE]
    X                            = mvnrnd(mu_cond', Sigma_cond, N);    % N x 8

    % Extract common parameters
    a_1_realizations             = max(X(:, 1), 0.01);
    Delta_1_realizations         = max(X(:, 2), 0);
    kappa_1_realizations         = max(X(:, 3), 0.0001);
    offset_1_realizations        = max(X(:, 4), 0);
    sigma_1_AP_realizations      = max(X(:, 5), 0.01);
    sigma_1_UE_realizations      = max(X(:, 6), 0.01);
    d_corr_1_AP_realizations     = max(X(:, 7), 1);
    d_corr_1_UE_realizations     = max(X(:, 8), 1);

    if near_far == "far"
        adjust_factor = 1;
    else
        adjust_factor = 0;
    end

    %% Output placeholders
    N_UE = size(p_UE_1, 1);
    N_AP = size(p_AP_1, 1);
    D = zeros(N_UE, N_AP);
    PL_mean = zeros(N, N_UE, N_AP);

    d_corner_1 = sqrt(sum((p_UE_1 - p_corner_1).^2, 2));
    d_1 = sqrt(sum((p_AP_1 - p_corner_1).^2, 2));

    for i = 1:N
        adjust_term = (1/pi * atan(1 * (d_corner_1 - 100)) + 0.5).*offset_1_realizations(i) .* adjust_factor;  % Apply adjust * side only when condition is true
        term1 = (Delta_1_realizations(i) + adjust_term) .* (1 - exp(-kappa_1_realizations(i) * d_1))';                            % [N x 1]
        term2 = (1/pi * atan(0.01 * (d_corner_1 - 70)) + 0.5) ...
                  .* a_1_realizations(i) .* 10 .* log10(d_1');  
        PL_mean(i, :, :) = term1 + term2;
    end

    %% Step 5: Add Correlated Shadowing
    if ndims(PL_mean) == 3
        % Multiple realizations
        [N, U, A] = size(PL_mean);
        PL_CUNEC_first = zeros(N, U, A);
        for i = 1:N
%             shadowing = generate_correlated_shadowing(p_AP_1, p_UE_1, ...
%                 sigma_1_UE_realizations(i), sigma_1_AP_realizations(i), ...
%                 d_corr_1_UE_realizations(i), d_corr_1_AP_realizations(i));
            shadowing = generate_correlated_shadowing(p_AP_1, p_UE_1, ...
                sigma_1_AP_realizations(i), sigma_1_UE_realizations(i), ...
                d_corr_1_AP_realizations(i), d_corr_1_UE_realizations(i));
            PL_CUNEC_first(i, :, :) = squeeze(PL_mean(i, :, :)) + shadowing;
            std_S(i) = std(shadowing(:));
            std_S_AP(i) = mean(std(shadowing));
            std_S_UE(i) = mean(std(shadowing'));
        end
    else
        % Single realization
        [U, A] = size(PL_mean);
        PL_CUNEC_first = zeros(U, A);
        for i = 1:N
%             shadowing = generate_correlated_shadowing(p_AP_1, p_UE_1, ...
%                 sigma_1_UE_realizations(i), sigma_1_AP_realizations(i), ...
%                 d_corr_1_UE_realizations(i), d_corr_1_AP_realizations(i));
            shadowing = generate_correlated_shadowing(p_AP_1, p_UE_1, ...
                sigma_1_AP_realizations(i), sigma_1_UE_realizations(i), ...
                d_corr_1_AP_realizations(i), d_corr_1_UE_realizations(i));
            PL_CUNEC_first(i,:) = PL_mean(i, :) + shadowing;
        end
    end
end