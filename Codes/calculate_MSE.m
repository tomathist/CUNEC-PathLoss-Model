%all_pathloss = path_loss_avg';
%all_pathloss = all_pathloss(:,5:60);
%mean_pathloss = mean(all_pathloss);

%mean_PL_AlphaBeta = mean_PL_AlphaBeta(:,5:60);
% Compute MSE between A and B
MSE_AB = mean((mean_PL_AlphaBeta - mean_pathloss).^2);

% Compute MSE between A and C
MSE_CUNEC = mean((mean_PL_CUNEC - mean_pathloss).^2);

% Display results
disp(['MSE between ray tracing and alpha-beta: ', num2str(MSE_AB)]);
disp(['MSE between ray tracing and CUNEC: ', num2str(MSE_CUNEC)]);

all_pathloss = mean_pathloss;

% Compute squared Euclidean distances (sum of squared errors)
D = pdist2(all_pathloss, pathloss_final, 'squaredeuclidean');

% Convert to mean squared error
MSE = D / size(all_pathloss, 2);

% For each row of all_pathloss, find the best match in pathloss_final
[bestMSE, bestIdx] = min(MSE, [], 2);

% Report lowest MSE overall
fprintf('Lowest MSE of CUNEC instances: %.6f\n', min(bestMSE));

% (Optional) also show which row of pathloss_final was the best match
%fprintf('Best match is row %d of pathloss_final\n', bestIdx(bestMSE == min(bestMSE)));