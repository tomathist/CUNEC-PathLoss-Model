function Ahat = nearestSPD(A)
    % From Higham, "Computing a nearest symmetric positive semidefinite matrix"
    B      = (A + A') / 2;
    [V, D] = eig(B);
    D      = diag(max(diag(D), 0));  % Remove negative eigenvalues
    Ahat   = V * D * V';
    Ahat   = (Ahat + Ahat') / 2;  % Ensure symmetry
end