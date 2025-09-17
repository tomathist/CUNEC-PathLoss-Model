%% demo_shadowing_joint_aniso.m
% Demo: generate and visualize a correlated shadowing field
% using gen_shadowing_joint_aniso (anisotropic exponential joint covariance).
%
% NOTE ON SIZE: M=N=100 => MN=10,000 and Sigma_joint is 10,000x10,000
% (~800 MB in doubles), plus Cholesky. This is heavy but often feasible on
% a modern machine. If you see memory/time issues, try M=N=40 first.

clear; clc;

%% 1) Reproducibility (optional)
rng(123);  % comment out for a fresh draw every run

%% 2) Define UE/AP positions (1-D lines along x, meters)
M = 100;                        % number of UEs
N = 100;                        % number of APs

% UEs: x = 1..M (y=z=0); APs: x = 1..N (y=z=0)
rUE = [(1:M).', zeros(M,1)];    % [M x 2] -> internally padded to 3D
rAP = [(1:N).', zeros(N,1)];    % [N x 2] -> internally padded to 3D

%% 3) Correlation distances and target stds (dB)
dcorrUE = 15;                   % UE-axis corr distance (meters)
dcorrAP = 5;                    % AP-axis corr distance (meters)
sigmaUE_mean = 6.5;             % target avg row-wise std (across APs)
sigmaAP_mean = 5.3;             % target avg col-wise std (across UEs)

%% 4) Generate shadowing field (dB)
tic;
S = gen_shadowing_joint_aniso(rUE, rAP, dcorrUE, dcorrAP, sigmaUE_mean, sigmaAP_mean);
t_gen = toc;

%% 5) Quick sanity checks
rowStdMean = mean(std(S,0,2));          % average std across rows (UE-wise)
colStdMean = mean(std(S,0,1));          % average std across cols (AP-wise)
Smean      = mean(S,'all');

fprintf('Generation time: %.2f s\n', t_gen);
fprintf('Mean(S): %.3f dB (should be near 0 in expectation)\n', Smean);
fprintf('Avg row std: %.3f dB (target %.3f)\n', rowStdMean, sigmaUE_mean);
fprintf('Avg col std: %.3f dB (target %.3f)\n', colStdMean, sigmaAP_mean);

%% 6) Visualize the field
figure('Name','Shadowing Field','Color','w');
imagesc(S);
axis image; box on;
xlabel('AP index'); ylabel('UE index');
title('Correlated Shadowing S (dB)');
cb = colorbar; ylabel(cb, 'Shadowing (dB)');

%% 7) Inspect one UE row and one AP column
iUE = round(M/2);   % middle UE
jAP = round(N/2);   % middle AP

figure('Name','Row slice','Color','w');
plot(1:N, S(iUE,:), 'LineWidth', 1.2);
grid on; xlabel('AP index'); ylabel('Shadowing (dB)');
title(sprintf('Row slice: UE %d across all APs', iUE));

figure('Name','Column slice','Color','w');
plot(1:M, S(:,jAP), 'LineWidth', 1.2);
grid on; xlabel('UE index'); ylabel('Shadowing (dB)');
title(sprintf('Column slice: AP %d across all UEs', jAP));

%% 8) Tips if this is slow or memory-heavy:
% - Reduce M and N (e.g., M=N=40) to shrink Sigma_joint.
% - If positions are on a grid and you can accept separability,
%   switch to a Kronecker sampler (UE factor * Z * AP factor') to avoid
%   building the full (MN)x(MN) covariance.
