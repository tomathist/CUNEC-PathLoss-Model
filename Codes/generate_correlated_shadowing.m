function shadowing = generate_correlated_shadowing(p_AP_0, p_UE_0, sigma_0_AP, sigma_0_UE, d_corr_0_AP, d_corr_0_UE)
% generate_correlated_shadowing
% Generates correlated shadowing field with specified standard deviations
% across UE rows and AP columns.
%
% Inputs:
%   p_AP_0         - [N_AP × 3] matrix of AP positions
%   p_UE_0         - [N_UE × 3] matrix of UE positions
%   sigma_0_AP     - Target standard deviation across APs
%   sigma_0_UE     - Target standard deviation across UEs
%   d_corr_0_AP    - Correlation distance across APs
%   d_corr_0_UE    - Correlation distance across UEs
%
% Output:
%   shadowing      - [N_UE × N_AP] correlated shadowing matrix (dB)

    N_AP = size(p_AP_0, 1);
    N_UE = size(p_UE_0, 1);

    if N_UE == 1  % Only one UE: generate 1D correlated shadowing across APs
        d_AP = vecnorm(p_AP_0 - p_AP_0(1,:), 2, 2);
        R_AP = exp(-abs(d_AP - d_AP') ./ d_corr_0_AP);
        [V_AP, D_AP] = eig(R_AP);
        z = randn(N_AP, 1);
        shadowing_1D = sqrt(sigma_0_AP) * V_AP * sqrt(D_AP) * z;
        shadowing = shadowing_1D';  % 1 × N_AP

    elseif N_AP == 1  % Only one AP: generate 1D correlated shadowing across UEs
        d_UE = vecnorm(p_UE_0 - p_UE_0(1,:), 2, 2);
        R_UE = exp(-abs(d_UE - d_UE') ./ d_corr_0_UE);
        [V_UE, D_UE] = eig(R_UE);
        z = randn(N_UE, 1);
        shadowing_1D = sqrt(sigma_0_UE) * V_UE * sqrt(D_UE) * z;
        shadowing = shadowing_1D';

    else  % Full 2D shadowing
        d_AP = vecnorm(p_AP_0 - p_AP_0(1,:), 2, 2);
        d_UE = vecnorm(p_UE_0 - p_UE_0(1,:), 2, 2);
        R_AP = exp(-abs(d_AP - d_AP') / d_corr_0_AP);
        R_UE = exp(-abs(d_UE - d_UE') / d_corr_0_UE);

        [V_AP, D_AP] = eig(R_AP);
        [V_UE, D_UE] = eig(R_UE);

        Z = randn(N_UE, N_AP);
        shadowing_base = V_UE * sqrt(D_UE) * Z * sqrt(D_AP) * V_AP';

        row_std = mean(std(shadowing_base, 0, 2));
        col_std = mean(std(shadowing_base, 0, 1));

        scale_factor = sqrt((sigma_0_UE / row_std) * (sigma_0_AP / col_std));
        shadowing = scale_factor * shadowing_base;
    end
end
