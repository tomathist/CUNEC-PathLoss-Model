function S_gen = gen_shadowing_joint_aniso(rUE, rAP, dcorrUE, dcorrAP, sigmaUE_mean, sigmaAP_mean)
%GEN_SHADOWING_JOINT_ANISO
% Generate correlated shadowing (dB) with anisotropic exponential joint covariance,
% then scale with S_gen = D_UE * S_tilde * D_AP (one-shot).
%
% Works for 1D, 2D, or 3D positions. Robust for M==1 or N==1.

    % ---- checks & setup ----
    if nargin < 6, error('Need 6 inputs: rUE, rAP, dcorrUE, dcorrAP, sigmaUE_mean, sigmaAP_mean'); end
    if dcorrUE <= 0 || dcorrAP <= 0, error('dcorrUE/dcorrAP must be > 0'); end
    if sigmaUE_mean < 0 || sigmaAP_mean < 0, error('sigmaUE_mean/sigmaAP_mean must be >= 0'); end
    eps0 = 1e-12;
    rUE = pad3(rUE);  [M,~] = size(rUE);
    rAP = pad3(rAP);  [N,~] = size(rAP);
    MN = M*N;

    % ---- pairwise distances ----
    DUE = pdist2(rUE, rUE, 'euclidean');   % [M x M]
    DAP = pdist2(rAP, rAP, 'euclidean');   % [N x N]

    % ---- joint covariance (MN x MN) ----
    termUE = kron((DUE./dcorrUE).^2, ones(N));
    termAP = kron(ones(M), (DAP./dcorrAP).^2);
    Tau = sqrt(termUE + termAP);
    Sigma_joint = exp(-Tau);
    Sigma_joint = (Sigma_joint + Sigma_joint.')/2;     % symmetry hygiene
    Sigma_joint(1:MN+1:end) = 1;                       % exact unit diag

    % ---- sample vec(S) ~ N(0, Sigma_joint) ----
    L = chol_psd(Sigma_joint);
    S_tilde = reshape(L * randn(MN,1), [M,N]);
    if ~(M==1 && N==1)  % keep single-sample case as-is
        S_tilde = S_tilde - mean(S_tilde,'all');
    end

    % ---- scale (handle degenerate axes!) ----
    if M==1 && N==1
        % Nothing to match statistically; return the scalar draw
        S_gen = S_tilde*sqrt(sigmaUE_mean)*sqrt(sigmaAP_mean);

    elseif M==1
        % Only one UE: match row-wise std (across APs); don't touch columns
        rstd = std(S_tilde, 0, 2);                    % scalar
        s = sigmaUE_mean / max(rstd, eps0);
        S_gen = s * S_tilde;

    elseif N==1
        % Only one AP: match column-wise std (across UEs); don't touch rows
        cstd = std(S_tilde, 0, 1);                    % scalar
        s = sigmaAP_mean / max(cstd, eps0);
        S_gen = S_tilde * s;

    else
        % Full 2D: one-shot row/col scaling
        row_std = std(S_tilde, 0, 2);                 % [M x 1]
        col_std = std(S_tilde, 0, 1);                 % [1 x N]
        D_UE = diag( sigmaUE_mean ./ max(row_std, eps0) );
        D_AP = diag( sigmaAP_mean ./ max(col_std, eps0) );
        S_gen = sqrt(D_UE) * S_tilde * sqrt(D_AP);
    end
end

% ---------- helpers ----------
function P3 = pad3(P)
    % Accept [Kx1], [Kx2], [Kx3], or row vector -> treat as [Kx1]
    if isempty(P), error('Position matrix cannot be empty'); end
    if size(P,1)==1 && size(P,2) > 3, P = P(:); end
    switch size(P,2)
        case 3, P3 = P;
        case 2, P3 = [P, zeros(size(P,1),1)];
        case 1, P3 = [P, zeros(size(P,1),2)];
        otherwise, error('Positions must have 1, 2, or 3 columns (or a row vector).');
    end
end

function L = chol_psd(A)
    % Cholesky with diagonal jitter if needed; eig fallback.
    try
        L = chol(A,'lower'); return;
    catch
    end
    n = size(A,1);
    jitter = 1e-12 * max(trace(A)/max(n,1), 1);
    for k=1:8
        try
            L = chol(A + jitter*eye(n),'lower'); return;
        catch
            jitter = jitter*10;
        end
    end
    [V,D]=eig((A+A')/2); D = max(D,0); L = V*sqrt(D);
end
