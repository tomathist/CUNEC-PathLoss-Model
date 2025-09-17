function PL_CUNEC_first = calc_PL_CUNEC_1st( ...
    p_AP_1, p_UE_1, p_corner_1, near_far, ...
    b, h, mu, sigma, C, R)
% calc_PL_CUNEC_1st
% First-order CUNEC path loss with *jointly* correlated log-normal shadowing.
% Correlated fields are generated via gen_shadowing_joint_aniso (external).
%
% Inputs
%   p_AP_1        [N_AP x 3]  AP 3D positions (m)
%   p_UE_1        [N_UE x 3]  UE 3D positions (m)
%   p_corner_1    [1 x 3]     3D position of the relevant corner (m)
%   b, h          scalars     building length, building height (m)
%   mu,sigma      10x1        parameter means and stds (CUNEC parameterization)
%   C             10x10       parameter correlation matrix (symmetric, PSD)
%   R             integer     # Monte Carlo realizations
%
% Outputs
%   PL_CUNEC_first   [R x N_UE x N_AP] path loss realizations (dB)
%
% Notes
% - Parameters are assumed ordered as:
%     [b, h, a1, Δ1, κ1, offset1, σ_S_AP, σ_S_UE, d_corr_AP, d_corr_UE].
% - We condition on fixed [b,h] (indices 1:2) and sample the remaining 8 params.
% - Mean PL form used here (per-realization i):
%     PL_mean = (Δ1_i + adjust_i(d_UE2corner)) .* (1 - exp(-κ1_i * d_AP2corner))' ...
%               + gate(d_UE2corner) .* a1_i .* 10*log10(max(d_AP2corner,1))'
%   where d_AP2corner is AP–corner distance and d_UE2corner is UE–corner distance.
% - Shadowing is added as a *dB* field with joint UE–AP correlation (anisotropic).
%
% Requirements
% - A function gen_shadowing_joint_aniso(rUE, rAP, dcorrUE, dcorrAP, sigmaUE, sigmaAP)
%   must be on path.

    %% --- Validation (mirror style of 0th-order) ---
    mustBeSize(p_AP_1,  [NaN,3], 'p_AP_1 must be N_AP x 3');
    mustBeSize(p_UE_1,  [NaN,3], 'p_UE_1 must be N_UE x 3');
    if ~isequal(size(p_corner_1),[1,3]), error('p_corner_1 must be 1x3.'); end

    mustBePosScalar(b,'b');  mustBePosScalar(h,'h');

    mustBeSize(mu,[10,1],'mu must be 10x1');
    mustBeSize(sigma,[10,1],'sigma must be 10x1');
    if any(sigma<=0), error('All sigma elements must be > 0.'); end
    if ~isequal(size(C),[10,10]) || ~issymmetric(C)
        error('C must be symmetric 10x10.');
    end
    if ~(isscalar(R) && R>=1 && R==round(R)), error('R must be a positive integer.'); end

    % near_far flag (accept string/char/logical); "far" enables the offset adjust term
    if isstring(near_far) || ischar(near_far)
        useAdjust = strcmpi(string(near_far),'far');
    elseif islogical(near_far)
        useAdjust = near_far;
    else
        error('near_far must be "near"|"far" (string/char) or logical.');
    end
    adjust_factor = double(useAdjust);

    N_AP = size(p_AP_1,1);
    N_UE = size(p_UE_1,1);

    %% --- Full covariance of parameters ---
    Sigma = diag(sigma) * C * diag(sigma);
    Sigma = (Sigma+Sigma')/2;

    % Indices: fix x_fixed = [b,h] (1:2), sample the rest (3:10)
    idx_fixed = [1,2];
    idx_free  = 3:10;

    mu_f = mu(idx_fixed);
    mu_u = mu(idx_free);
    S_ff = Sigma(idx_fixed, idx_fixed);
    S_uu = Sigma(idx_free,  idx_free);
    S_uf = Sigma(idx_free,  idx_fixed);
    S_fu = Sigma(idx_fixed, idx_free);

    % Condition on fixed [b,h]
    x_fixed = [b; h];
    S_ff = makeSPD(S_ff);
    S_uu = makeSPD(S_uu);

    mu_cond    = mu_u + S_uf / S_ff * (x_fixed - mu_f);
    Sigma_cond = S_uu - S_uf / S_ff * S_fu;
    Sigma_cond = makeSPD(Sigma_cond);

    %% --- Sample the 8 free parameters R times ---
    % Order (idx_free): [a1, Δ1, κ1, offset1, σ_AP, σ_UE, d_corr_AP, d_corr_UE]
    X = mvnrnd(mu_cond.', Sigma_cond, R);      % [R x 8]

    a1      = max(X(:,1), 0.01);               % slope-like coefficient
    Delta1  = X(:,2);                           % intercept-ish term (can be negative)
    kappa1  = max(X(:,3), 1e-4);                % exponential rate
    offset1 = max(X(:,4), 0);                   % offset (≥0)
    sAP     = max(X(:,5), 0.01);                % AP shadowing std (dB)
    sUE     = max(X(:,6), 0.01);                % UE shadowing std (dB)
    dAP     = max(X(:,7), 1);                   % AP corr distance (m)
    dUE     = max(X(:,8), 1);                   % UE corr distance (m)

    %% --- Pairwise UE–AP distances (m), to return like 0th-order ---
    AP2 = sum(p_AP_1.^2,2);               % [N_AP x 1]
    UE2 = sum(p_UE_1.^2,2);               % [N_UE x 1]
    D2  = UE2 + AP2.' - 2*(p_UE_1 * p_AP_1.');   % [N_UE x N_AP]
    D2  = max(D2, 0);

    %% --- Distances to the corner (used by 1st-order mean model) ---
    % UE -> corner: [N_UE x 1], AP -> corner: [N_AP x 1]
    d_UE2corner = sqrt(sum((p_UE_1 - p_corner_1).^2, 2));  % [N_UE x 1]
    d_AP2corner = sqrt(sum((p_AP_1 - p_corner_1).^2, 2));  % [N_AP x 1]
    d_AP2corner_clamp = max(d_AP2corner, 1);               % for log10 safety

    %% --- First-order mean PL (vectorized over R) ---
    % Gate/adjust functions (match your existing shapes):
    % - gate(d_UE2corner): smooth gate in [0,1]
    % - adjust(d_UE2corner): used only when near_far=="far"
    gate_UE  = (1/pi * atan(0.01 * (d_UE2corner - 70)) + 0.5);   % [N_UE x 1]
    adjustUE = (1/pi * atan(1.00 * (d_UE2corner - 100)) + 0.5);  % [N_UE x 1]

    % Precompute AP terms
    log_dAP  = 10*log10(d_AP2corner_clamp);                      % [N_AP x 1]
    exp_AP   = (1 - exp(-kappa1.*d_AP2corner'));                % [R x N_AP], broadcast κ per r

    % Allocate
    PL_CUNEC_first = zeros(R, N_UE, N_AP);

    % Build mean + shadowing per realization
    for r = 1:R
        % term1: (Δ1 + adjust*offset) * (1 - exp(-κ * d_AP2corner))'
        adjust_term_r = adjust_factor .* (adjustUE .* offset1(r));             % [N_UE x 1]
        term1_r = (Delta1(r) + adjust_term_r) .* exp_AP(r,:);                  % [N_UE x N_AP] via implicit expansion

        % term2: gate(d_UE2corner) * a1 * 10*log10(d_AP2corner)'
        term2_r = gate_UE .* (a1(r) .* log_dAP.');                              % [N_UE x N_AP]

        PL_mean_r = term1_r + term2_r;                                          % [N_UE x N_AP]

        % Shadowing per realization (UE-AP jointly correlated, anisotropic)
        S = gen_shadowing_joint_aniso( ...
                p_UE_1, p_AP_1, ...
                dUE(r), dAP(r), ...
                sUE(r), sAP(r));                                                % [N_UE x N_AP]

        PL_CUNEC_first(r,:,:) = PL_mean_r + S;
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
