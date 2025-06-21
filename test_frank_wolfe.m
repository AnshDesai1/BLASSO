clear all;
close all;
clc;

%% Generate Signal
N = 4;
d = 0.5;
K = 300;
theta = deg2rad([0, 15, 20]);
az_angs=-40:.1:40;
SNR = [10, 15, 20];
SssDim = 3;
x = signal_gen(N, d, rad2deg(theta), SNR, K);

% Use average signal as input
R = (x * x') / K;            % Sample covariance
[U, S, ~] = svd(R);          % Principal eigenvector
y_obs = U(:,1) * sqrt(S(1)); % Dominant source direction
y_grid = (0:N-1)';

%% Frank-Wolfe
lambda = 0.01 * norm(y_obs)^2 / length(y_obs);
max_iter = 100;
tol = 1e-8;
[theta_est, a_est] = sliding_frank_wolfe(y_obs, y_grid, d, lambda, max_iter, tol);

%% MUSIC
Vs=U(:,1:SssDim);            % Signal Subspace
Vn=U(:,SssDim+1:end);        % Noise Subspace
eigValSpec=10*log10(diag(S)); % Eigenvalue Spectrum

A=linear_dir_vec(N,d,az_angs); %Array Manifold of Steering Vectors
I=eye(N);
for ii=1:length(az_angs)
    a=A(:,ii); % Assess spectrum at the angle value corresponding to this particular steering vector
    Pmus(ii)=1/(a'*Vn*Vn'*a); % Calculate MUSIC Spectrum Using Noise Subspace
    Pmus_signalSubSpace(ii)=1/(a'*(I-Vs*Vs')*a); % Calculate Spectrum Using Signal Subspace
end
[~, locs] = findpeaks(10*log10(abs(Pmus_signalSubSpace)), 'SortStr', 'descend', 'NPeaks', SssDim);
theta_est_music = deg2rad(az_angs(locs));

%% Compare
fprintf('\nTrue angle locations (radians):\n');
disp(theta);
fprintf('\nEstimated angle locations (radians):\n');
disp(theta_est');
fprintf('MUSIC estimated angles (deg):\n');
disp(theta_est_music);

figure;
hold on;
xlim([-pi, pi]);
ylim([0, 1.1]);
title('True vs Estimated Spike Locations');
xlabel('\theta (radians)');
ylabel('Normalized Spike Height');

% Plot one line of each type, store handles
h1 = plot([theta(1), theta(1)], [0, 1], 'b--', 'LineWidth', 2); % True
for i = 2:length(theta)
    plot([theta(i), theta(i)], [0, 1], 'b--', 'LineWidth', 2);
end

h2 = plot([theta_est(1), theta_est(1)], [0, 0.8], 'r-', 'LineWidth', 2); % SFW
for i = 2:length(theta_est)
    plot([theta_est(i), theta_est(i)], [0, 0.8], 'r-', 'LineWidth', 2);
end

h3 = plot([theta_est_music(1), theta_est_music(1)], [0, 0.6], 'm-.', 'LineWidth', 2); % MUSIC
for i = 2:length(theta_est_music)
    plot([theta_est_music(i), theta_est_music(i)], [0, 0.6], 'm-.', 'LineWidth', 2);
end

% Proper legend using representative handles
legend([h1, h2, h3], {'Ground Truth', 'SFW Estimated', 'MUSIC'}, 'Location', 'northeast');

grid on;