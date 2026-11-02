%MAIN  Engineering Mathematics - MATLAB homework 2: the heat equation with a
%      point source (pdepe) and the 2-D wave equation (finite differences).

clear; close all; clc;
addpath(fileparts(mfilename('fullpath')));

%% Q1  u_t = alpha u_xx + f,  u_x(0,t) = u_x(1,t) = 0,  u(x,0) = 0,  0 < t < 10
% f = 100 at x = 0.5: modelled as a band of width h = 0.01 (strength 100 h = 1),
% which keeps the source independent of the grid.
h = 0.01;
x = unique([linspace(0, 1, 401), 0.5 - h/2, 0.5 + h/2]);   % mesh points at the band edges
t = linspace(0, 10, 201);
alphas = [0.001 0.01 0.1 1];
fig = figure('Visible', 'off', 'Position', [100 100 1100 420]);
for k = 1:numel(alphas)
    a = alphas(k);
    pde = @(x, t, u, dudx) deal(1, a*dudx, 100*(abs(x - 0.5) < h/2));
    sol = pdepe(0, pde, @(x) 0, @(xl, ul, xr, ur, t) deal(0, 1, 0, 1), x, t);
    total = trapz(x, sol(end, :));                          % insulated rod: int u dx = 1 * t
    fprintf('Q1  alpha = %-6g  u(0.5, 10) = %7.2f   u(0, 10) = %6.2f   int u dx = %.3f (expected 10)\n', ...
            a, interp1(x, sol(end, :), 0.5), sol(end, 1), total);
    subplot(1, 2, 1); plot(x, sol(end, :), 'LineWidth', 1.3); hold on;
    subplot(1, 2, 2); plot(t, sol(:, x == 0.5), 'LineWidth', 1.3); hold on;
    if a == 0.01
        f3 = figure('Visible', 'off', 'Position', [100 100 700 500]);
        surf(x(1:4:end), t, sol(:, 1:4:end), 'EdgeColor', 'none'); view(135, 30);
        xlabel('x'); ylabel('t'); zlabel('u'); title('u(x,t) for \alpha = 0.01'); colormap hot;
        save_figure(f3, 'q1_surface');
        figure(fig);
    end
end
subplot(1, 2, 1); grid on; xlabel('x'); ylabel('u(x, 10)'); set(gca, 'YScale', 'log'); ylim([1e-2 1e3]);
title('Profile at t = 10 s'); legend(arrayfun(@(a) sprintf('\\alpha = %g', a), alphas, 'UniformOutput', false));
subplot(1, 2, 2); grid on; xlabel('t (s)'); ylabel('u(0.5, t)'); title('Temperature at the source');
save_figure(fig, 'q1_alpha_sweep');

%% Q2  u_tt = u_xx + u_yy on [-5,5]^2, u = 0 on the boundary, u(x,y,0) = exp(-100(x^2+y^2)), u_t = 0
L = 5;  n = 401;
[X, Y] = meshgrid(linspace(-L, L, n));                    % a true 2-D grid
dx = X(1, 2) - X(1, 1);
dt = dx/2;                                                % CFL in 2-D: c dt / dx <= 1/sqrt(2)
lap = @(U) [zeros(1, n); zeros(n-2, 1), U(3:end, 2:end-1) + U(1:end-2, 2:end-1) + ...
            U(2:end-1, 3:end) + U(2:end-1, 1:end-2) - 4*U(2:end-1, 2:end-1), zeros(n-2, 1); zeros(1, n)] / dx^2;
U0 = exp(-100*(X.^2 + Y.^2));
U1 = U0 + 0.5*dt^2*lap(U0);                              % second-order start with u_t = 0
% Discrete energy that the leapfrog scheme conserves exactly: kinetic term from the
% time difference, potential term from the product of gradients at the two levels.
energy = @(Un, Up) 0.5*(sum(((Un - Up)/dt).^2, 'all') + ...
         (sum(diff(Un, 1, 2).*diff(Up, 1, 2), 'all') + sum(diff(Un, 1, 1).*diff(Up, 1, 1), 'all'))/dx^2)*dx^2;
T_end = 8;  steps = round(T_end/dt);
snap_t = [0 1 3 5 6.5 8];  snaps = {};  E = zeros(1, steps);
Up = U0;  Uc = U1;
for s = 1:steps
    Un = 2*Uc - Up + dt^2*lap(Uc);
    Un([1 end], :) = 0;  Un(:, [1 end]) = 0;
    E(s) = energy(Un, Uc);
    Up = Uc;  Uc = Un;
    if any(abs(s*dt - snap_t) < dt/2), snaps{end+1} = Uc; end %#ok<SAGROW>
end
fprintf('Q2  %d x %d grid, dt = %.4f, energy drift over %g s: %.2e (relative)\n', n, n, dt, T_end, ...
        (max(E) - min(E))/mean(E));
fig = figure('Visible', 'off', 'Position', [100 100 1150 650]);
snaps = [{U0}, snaps];
for k = 1:min(6, numel(snaps))
    subplot(2, 3, k); imagesc(X(1, :), Y(:, 1), snaps{k}); axis image; colormap(gca, turbo);
    clim([-0.02 0.02]); title(sprintf('t = %.1f', snap_t(k))); xlabel('x'); ylabel('y');
end
sgtitle('2-D wave equation: an expanding ring that reflects (with sign flip) at the fixed walls');
save_figure(fig, 'q2_wave_snapshots');
