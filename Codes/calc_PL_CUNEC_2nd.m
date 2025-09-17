function PL_CUNEC_second = calc_PL_CUNEC_2nd( ...
    p_AP_2, p_UE_2, p_corner_2, ...
    b, mu, sigma, C, R)
% calc_PL_CUNEC_2nd
% Second-order CUNEC path loss with *jointly* correlated log-normal shadowing.
% Correlated fields are generated via gen_shadowing_joint_aniso (external).
%
% Inputs
%   p_AP_2        [N_AP x 3]  AP 3D positions (m)
%   p_UE_2        [N_UE x 3]  UE 3D positions (m)
%   p_corner_2    [1 x 3]     3D position of the second-order corner (m)
%   b             scalar      building length (m)  [fixed/conditioned]
%   mu,sigma      6x1         parameter means and stds
%   C             6x6         parameter correlation matrix (symmetric, PSD)
%   R             integer     # Monte Carlo realizations
%
% Outputs
%   D                [N_UE x N_AP] Euclidean UE–AP distances (m)
%   PL_CUNEC_second  [R x N_UE x N_AP] second-order path loss realizations (dB)
%
% Notes
% - Assumed parameter order:
%       [ b, a2, σ_S_AP, σ_S_UE, d_corr_AP, d_corr_UE ].
%   We condition on fixed b (index 1) and sample the remaining 5 parameters.
% - Second-order mean PL (per realization r):
%       PL_mean_r = a2(r) * 10*log10( max(d_AP2corner, 1) )'
%   replicated across the UE dimension and then summed with a jointly
%   correlated dB shadowing field S (UE–AP).
%
% Requirements
% - A function gen_shadowing_joint_aniso(rUE, rAP, dcorrUE, dcorrAP, sigmaUE, sigmaAP)
%   must be on path.

    %% --- Validation (consistent with 0th/1st order) ---
    mustBeSize(p_AP_2, [NaN,3], 'p_AP_2 must be N_AP x 3');
    mustBeSize(p_UE_2, [NaN,3], 'p_UE_2 must be N_UE x 3');
    if ~isequal(size(p_corner_2),[1,3]), error('p_corner_2 must be 1x3.'); end

    mustBePosScalar(b,'b');

    mustBeSize(mu,[6,1],'mu must be 6x1');
    mustBeSize(sigma,[6,1],'sigma must be 6x1');
    if any(sigma<=0), error('All sigma elements must be > 0.'); end
    if ~isequal(size(C),[6,6]) || ~issymmetric(C)
        error('C must be symmetric 6x6.');
    end
    if ~(isscalar(R) && R>=1 && R==round(R)), error('R must be a positive integer.'); end

    N_AP = size(p_AP_2,1);
    N_UE = size(p_UE_2,1);

    %% --- Full covariance of parameters and conditioning on fixed b ---
    Sigma = diag(sigma) * C * diag(sigma);
    Sigma = (Sigma+Sigma')/2;

    idx_fixed = 1;            % b
    idx_free  = 2:6;          % [a2, σ_AP, σ_UE, d_corr_AP, d_corr_UE]

    mu_f = mu(idx_fixed);
    mu_u = mu(idx_free);
    S_ff = Sigma(idx_fixed, idx_fixed);
    S_uu = Sigma(idx_free,  idx_free);
    S_uf = Sigma(idx_free,  idx_fixed);
    S_fu = Sigma(idx_fixed, idx_free);

    x_fixed = b;
    S_ff = makeSPD(S_ff);
    S_uu = makeSPD(S_uu);

    mu_cond    = mu_u + S_uf / S_ff * (x_fixed - mu_f);
    Sigma_cond = S_uu - S_uf / S_ff * S_fu;
    Sigma_cond = makeSPD(Sigma_cond);

    %% --- Sample the 5 free parameters R times ---
    % Order (idx_free): [a2, σ_AP, σ_UE, d_corr_AP, d_corr_UE]
    X = mvnrnd(mu_cond.', Sigma_cond, R);     % [R x 5]

    a2  = max(X(:,1), 0.01);                  % slope-like coefficient
    sAP = max(X(:,2), 0.01);                  % AP shadowing std (dB)
    sUE = max(X(:,3), 0.01);                  % UE shadowing std (dB)
    dAP = max(X(:,4), 1);                     % AP corr distance (m)
    dUE = max(X(:,5), 1);                     % UE corr distance (m)

    %% --- UE–AP distances (to return, like other orders) ---
    AP2 = sum(p_AP_2.^2,2);                         % [N_AP x 1]
    UE2 = sum(p_UE_2.^2,2);                         % [N_UE x 1]
    D2  = UE2 + AP2.' - 2*(p_UE_2 * p_AP_2.');      % [N_UE x N_AP]
    D2  = max(D2, 0);
    D   = sqrt(D2);

    %% --- AP distance to the 2nd-order corner (drives the mean term) ---
    d_AP2corner        = sqrt(sum((p_AP_2 - p_corner_2).^2, 2));  % [N_AP x 1]
    d_AP2corner_clamp  = max(d_AP2corner, 1);                      % log10 safety
    log_dAP2corner     = 10*log10(d_AP2corner_clamp);              % [N_AP x 1]

    %% --- Mean PL and correlated shadowing (per realization) ---
    PL_CUNEC_second = zeros(R, N_UE, N_AP);

    for r = 1:R
        % Mean (replicate across UEs): [N_UE x N_AP]
        PL_mean_r = a2(r) .* (ones(N_UE,1) * log_dAP2corner.');    % broadcast across UEs

        % Jointly correlated shadowing (anisotropic UE–AP field), dB
        S = gen_shadowing_joint_aniso( ...
                p_UE_2, p_AP_2, ...
                dUE(r), dAP(r), ...
                sUE(r), sAP(r));                                   % [N_UE x N_AP]

        PL_CUNEC_second(r,:,:) = PL_mean_r + S;
    end
end

% ------------------- helpers -------------------

function A = makeSPD(A)
    % Make symmetric positive definite (Higham-style fix)
    A = (A+A')/2;
    [V,D] = eig(A);
    d = diag(D);
    d = max(d, eps(max(d)));
    A = V*diag(d)*V';
    A = (A+A')/2;
end

function mustBePosScalar(x,name)
    if ~(isscalar(x) && isfinite(x) && x>0)
        error('%s must be a positive scalar.', name);
    end
end

function mustBeSize(X,sz,msg)
    if ~(size(X,2)==sz(2) && (isnan(sz(1)) || size(X,1)==sz(1)))
        error(msg);
    end
end
