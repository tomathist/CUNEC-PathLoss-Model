function PL_CUNEC_second = calc_PL_CUNEC_2nd(p_AP_2, p_UE_2, p_corner_2, b, mu, sigma, C, N)
    %% Step 1: Full Covariance Matrix
    Sigma                        = diag(sigma) * C * diag(sigma);

    %% Step 2: Condition on Fixed Width and Height
    x_given                      = [b];                          % Fixed values
    idx_fixed                    = [1];                          % Indices for building, height
    idx_free                     = setdiff(1:6, idx_fixed);         % Indices to sample

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
    X                            = mvnrnd(mu_cond', Sigma_cond, N);

    % Extract common parameters
    a_2_realizations             = max(X(:, 1), 0.01);
    sigma_2_AP_realizations      = max(X(:, 2), 0.01);
    sigma_2_UE_realizations      = max(X(:, 3), 0.01);
    d_corr_2_AP_realizations     = max(X(:, 4), 1);
    d_corr_2_UE_realizations     = max(X(:, 5), 1);

    %% Output placeholders
    N_UE = size(p_UE_2, 1);
    N_AP = size(p_AP_2, 1);
    PL_mean = zeros(N, N_UE, N_AP);

    d_2 = sqrt(sum((p_AP_2 - p_corner_2).^2, 2));
    for i = 1:N
        PL_mean(i,:,:) = a_2_realizations(i) .* 10 .* log10(d_2');  
    end

    PL_mean = squeeze(PL_mean);
    %% Step 5: Add Correlated Shadowing
    if ndims(PL_mean) == 3
        % Multiple realizations
        [N, U, A] = size(PL_mean);
        PL_CUNEC_second = zeros(N, U, A);
        for i = 1:N
            shadowing = generate_correlated_shadowing(p_AP_2, p_UE_2, ...
                sigma_2_UE_realizations(i), sigma_2_AP_realizations(i), ...
                d_corr_2_UE_realizations(i), d_corr_2_AP_realizations(i));
%             shadowing = generate_correlated_shadowing(p_AP_2, p_UE_2, ...
%                 sigma_2_AP_realizations(i), sigma_2_UE_realizations(i), ...
%                 d_corr_2_AP_realizations(i), d_corr_2_UE_realizations(i));
            PL_CUNEC_second(i, :, :) = squeeze(PL_mean(i, :, :)) + shadowing;
        end
    else
        % Single realization
        [U, A] = size(PL_mean);
        PL_CUNEC_second = zeros(U, A);
        for i = 1:N
            shadowing = generate_correlated_shadowing(p_AP_2, p_UE_2, ...
                sigma_2_UE_realizations(i), sigma_2_AP_realizations(i), ...
                d_corr_2_UE_realizations(i), d_corr_2_AP_realizations(i));
%             shadowing = generate_correlated_shadowing(p_AP_2, p_UE_2, ...
%                 sigma_2_AP_realizations(i), sigma_2_UE_realizations(i), ...
%                 d_corr_2_AP_realizations(i), d_corr_2_UE_realizations(i));
            PL_CUNEC_second(i,:) = PL_mean(i,:) + shadowing;
        end
    end
end