function [D, PL_CUNEC_zeroth] = calc_PL_CUNEC_0th( ...
    p_AP_0, p_UE_0, b, w, h, FSPL_1m, mu, sigma, C, ...
    R)
% calc_PL_CUNEC_0th
% Zeroth-order CUNEC path loss with *jointly* correlated log-normal shadowing.
% Correlated fields are generated via gen_shadowing_joint_aniso (external).
%
% Inputs
%   p_AP_0        [N_AP x 3]  AP 3D positions (m)
%   p_UE_0        [N_UE x 3]  UE 3D positions (m)
%   b,w,h         scalars     building length, street width, building height (m)
%   FSPL_1m       scalar      free-space path loss @1 m (dB)
%   mu,sigma      9x1         parameter means and stds (CUNEC parameterization)
%   C             9x9         parameter correlation matrix (symmetric, PSD)
%   R             integer     # Monte Carlo realizations
%
% Outputs
%   D                [N_UE x N_AP] Euclidean distances (m)
%   PL_CUNEC_zeroth  [R x N_UE x N_AP] path loss realizations (dB)
%
% Notes
% - CUNEC "zeroth-order" mean PL:  FSPL_1m + Δ0 + b0 * 10*log10(D_clamped)
% - Shadowing is added as a *dB* field with joint UE-AP correlation.
% - The conditional MVN step treats [b,w,h] as fixed and samples the remaining 6 params:
%     [b0, Δ0, σ_S_AP, σ_S_UE, d_corr_AP, d_corr_UE].
%
% Requirements
% - A function gen_shadowing_joint_aniso(rUE, rAP, dcorrUE, dcorrAP, sigmaUE_mean, sigmaAP_mean)
%   must be on path (your simple one-shot scaling version is fine).

    %% --- Validation ---
    mustBeSize(p_AP_0,  [NaN,3], 'p_AP_0 must be N_AP x 3');
    mustBeSize(p_UE_0,  [NaN,3], 'p_UE_0 must be N_UE x 3');
    mustBePosScalar(b,'b');  mustBePosScalar(w,'w');  mustBePosScalar(h,'h');
    mustBePosScalar(FSPL_1m,'FSPL_1m');
    mustBeSize(mu,[9,1],'mu must be 9x1'); mustBeSize(sigma,[9,1],'sigma must be 9x1');
    if any(sigma<=0), error('All sigma elements must be > 0.'); end
    if ~isequal(size(C),[9,9]) || ~issymmetric(C)
        error('C must be symmetric 9x9.');
    end
    if ~(isscalar(R) && R>=1 && R==round(R)), error('R must be a positive integer.'); end
    
    N_AP = size(p_AP_0,1);
    N_UE = size(p_UE_0,1);

    %% --- Full covariance of parameters ---
    Sigma = diag(sigma) * C * diag(sigma);
    Sigma = (Sigma+Sigma')/2;

    % Indices: fix x_fixed = [b,w,h] (1:3), sample the rest (4:9)
    idx_fixed = [1,2,3];
    idx_free  = 4:9;

    mu_f = mu(idx_fixed);
    mu_u = mu(idx_free);
    S_ff = Sigma(idx_fixed, idx_fixed);
    S_uu = Sigma(idx_free,  idx_free);
    S_uf = Sigma(idx_free,  idx_fixed);
    S_fu = Sigma(idx_fixed, idx_free);

    % Condition on fixed [b,w,h]
    x_fixed = [b; w; h];
    % Regularize S_ff if needed (invertibility)
    S_ff = makeSPD(S_ff);
    S_uu = makeSPD(S_uu);

    mu_cond    = mu_u + S_uf / S_ff * (x_fixed - mu_f);
    Sigma_cond = S_uu - S_uf / S_ff * S_fu;
    Sigma_cond = makeSPD(Sigma_cond);

    %% --- Sample the 6 free parameters R times ---
    % Order: [b0, Δ0, σ_S_AP, σ_S_UE, d_corr_AP, d_corr_UE]
    X = mvnrnd(mu_cond.', Sigma_cond, R);       % R x 6
    b0      = max(X(:,1), 0.01);                % path loss slope (lower bounded)
    Delta0  = X(:,2);                            % intercept offset
    sAP     = max(X(:,3), 0);                    % shadowing std across APs (dB)
    sUE     = max(X(:,4), 0);                    % shadowing std across UEs (dB)
    dAP     = max(X(:,5), 1);                    % corr distance across APs (m)
    dUE     = max(X(:,6), 1);                    % corr distance across UEs (m)

    %% --- Pairwise UE–AP distances (m) ---
    % Use base MATLAB (no pdist2): ||u-a|| via norm identities, clamped at >=1 m for FSPL_1m
    AP2 = sum(p_AP_0.^2,2);          % [N_AP x 1]
    UE2 = sum(p_UE_0.^2,2);          % [N_UE x 1]
    D2  = UE2 + AP2.' - 2*(p_UE_0 * p_AP_0.');   % [N_UE x N_AP]
    D2  = max(D2, 0);                                % numerical safety
    D   = sqrt(D2);
    Dclamp = max(D, 1);                               % avoid log10(0); matches FSPL @ >= 1 m

    %% --- Zeroth-order mean PL (vectorized over R) ---
    % PL_mean(r,i,j) = FSPL_1m + Delta0(r) + b0(r) * 10*log10(Dclamp(i,j))
    L10D = 10*log10(Dclamp);                        % [N_UE x N_AP]
    PL_CUNEC_zeroth = zeros(R, N_UE, N_AP);

    for r = 1:R
        PL_mean_r = FSPL_1m + Delta0(r) + b0(r)*L10D;     % [N_UE x N_AP]

        % Shadowing per realization (can be heavy for large grids)
        S = gen_shadowing_joint_aniso( ...
                p_UE_0, p_AP_0, ...
                dUE(r), dAP(r), ...
                sUE(r), sAP(r));

        PL_CUNEC_zeroth(r,:,:) = PL_mean_r + S;
    end
end

% ------------------- helpers -------------------

function A = makeSPD(A)
    % Make symmetric positive definite (Higham-style fix)
    A = (A+A')/2;
    [V,D] = eig(A);
    D = diag(D);
    D = max(D, eps(max(D)));   % push up tiny/negative eigenvalues
    A = V*diag(D)*V';
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

