
clear; clc; close all

%% Ring-shaped domain: 0.2^2 < x^2 + y^2 <= 1^2

theta_grid = linspace(0, 2*pi, 80);
rho = linspace(0.2, 1.00, 15);   % avoid exactly r = 1 if desired

[Theta, Rho] = meshgrid(theta_grid, rho);

x_ = Rho .* cos(Theta);
y_ = Rho .* sin(Theta);

xx = x_(:);
yy = y_(:);
S = [xx, yy];
m = zeros(size(xx));

%% Inner and outer boundaries for plotting

theta_bd = 0:0.01:2*pi;

bd_x_inner = 0.2 * cos(theta_bd);
bd_y_inner = 0.2 * sin(theta_bd);

bd_x_outer = 1 * cos(theta_bd);
bd_y_outer = 1 * sin(theta_bd);

%% Squared exponential covariance

sqexp = @(A,B,lambda) exp(-pdist2(A,B).^2 ./ (2*lambda^2));

%% cGRF covariance function for annulus

outer_ring_cGRF_cov = @(S, lambda, h) local_outer_ring_cov(S, lambda, h, sqexp);
annulus_cGRF_cov = @(S, lambda, h) local_annulus_cov(S, lambda, h, sqexp);

%% Visualization of cGRF with constraints on inner and outer circles (Paper)

rng(2)

lambda = 0.5; h = 0.5;
K = annulus_cGRF_cov(S, lambda, h);
K(isnan(K)) = 0;
f = mvnrnd(m, K, 2);

% figure;
t = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for j = 1:2
   nexttile(j)
    surf(x_, y_, reshape(f(j,:), size(x_)), 'FaceColor', 'interp', 'FaceAlpha', 0.8)
    zlim([min(f(:)) max(f(:))]);
    hold on
    plot3(bd_x_inner, bd_y_inner, zeros(size(bd_x_inner)), '-', 'Color', 'red', 'LineWidth', 3)
    plot3(bd_x_outer, bd_y_outer, zeros(size(bd_x_outer)), '-', 'Color', 'red', 'LineWidth', 3)
    % view(2)
    axis square
    if j == 2
        xlabel("x_1")
        ylabel("x_2")
    end
    zlabel("u^A")
    % if j == 1
    %     zlabel("u^A")
    % end
    % xlabel("x_1")
    % ylabel("x_2")
    xticks([-1 0 1])
    yticks([-1 0 1])

    % colorbar;

    fontsize(20, "points")
    hold off
end

%% Sensitivity Analysis (Paper)

rng(10)

lambda_fixed = 0.5;
strengths = [0.1, 0.2, 0.5, 1, 2];

nS = numel(strengths);
f_strength = zeros(2, length(xx), nS);
v_strength = zeros(1, length(xx), nS);

for i = 1:nS
    h = strengths(i);

    K = outer_ring_cGRF_cov(S, lambda_fixed, h);

    K(isnan(K)) = 0;
    K = (K + K') / 2;

    f_strength(:,:,i) = mvnrnd(m, K, 2);
    v_strength(1,:,i) = diag(K);
end

clims_f = [min(f_strength(:)), max(f_strength(:))];
clims_v = [min(v_strength(:)), max(v_strength(:))];

figure;
t = tiledlayout(2, nS, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:nS
    for j = 1:2
        nexttile((j-1)*nS + i)

        surf(x_, y_, reshape(f_strength(j,:,i), size(x_)), ...
            'FaceColor', 'interp', 'EdgeColor', 'none')
        hold on
        
        % plot3(bd_x_inner, bd_y_inner, zeros(size(bd_x_inner)), ...
        %     '-', 'Color', 'red', 'LineWidth', 3)
        plot3(bd_x_outer, bd_y_outer, zeros(size(bd_x_outer)), ...
            '-', 'Color', 'red', 'LineWidth', 3)

        view(2)
        axis square
        clim(clims_f)

        if j == 1
            title(sprintf("h = %.2g", strengths(i)))
            xticks([])
        else
            xlabel("x_1")
        end

        if i == nS
            colorbar;
        end

        if i == 1
            ylabel("x_2")
        else
            yticks([])
        end

        fontsize(15, "points")
        hold off
    end
end

sgtitle(sprintf("Samples from cGRFs on annulus with fixed $\\lambda = %.2g$ and varying h", lambda_fixed), ...
    "Interpreter", "latex", 'FontSize', 20)

figure;
t = tiledlayout(1, nS, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:nS
    nexttile(i)

    surf(x_, y_, reshape(v_strength(1,:,i), size(x_)), ...
        'FaceColor', 'interp', 'EdgeColor', 'none')
    hold on

    % plot3(bd_x_inner, bd_y_inner, zeros(size(bd_x_inner)), ...
    %     '-', 'Color', 'red', 'LineWidth', 3)
    plot3(bd_x_outer, bd_y_outer, zeros(size(bd_x_outer)), ...
        '-', 'Color', 'red', 'LineWidth', 3)

    view(2)
    axis square
    clim(clims_v)

    title(sprintf("h = %.2g", strengths(i)))
    xlabel("x_1")

    if i == nS
        colorbar;
    end

    if i == 1
        ylabel("x_2")
    else
        yticks([])
    end

    fontsize(15, "points")
    hold off
end

sgtitle(sprintf("Pointwise variance on annulus with fixed $\\lambda = %.2g$ and varying h", lambda_fixed), ...
    "Interpreter", "latex", 'FontSize', 20)


%% Functions

function Kc = local_annulus_cov(S, lambda, h, sqexp)

    N = size(S,1);

    r = sqrt(sum(S.^2, 2));

    F1 = 0.2 * S ./ r;
    F2 = 1.0 * S ./ r;


    KSS = sqexp(S, S, lambda);

    KSF1 = sqexp(S, F1, lambda);
    KSF2 = sqexp(S, F2, lambda);

    KF1S = KSF1';
    KF2S = KSF2';

    KF1F1 = sqexp(F1, F1, lambda);
    KF1F2 = sqexp(F1, F2, lambda);
    KF2F1 = sqexp(F2, F1, lambda);
    KF2F2 = sqexp(F2, F2, lambda);

    a = ones(N,1);
    d = ones(N,1);
    b = diag(KF1F2); 

    detK = a .* d - b.^2;

    inv11 = d ./ detK; % M
    inv12 = -b ./ detK;
    inv22 = a ./ detK;

    ksf1 = diag(sqexp(S, F1, h)); % v
    ksf2 = diag(sqexp(S, F2, h));

    A1 = ksf1 .* inv11 + ksf2 .* inv12; % w
    A2 = ksf1 .* inv12 + ksf2 .* inv22;

    Kc = KSS ...
        - A1 .* KF1S ...
        - A2 .* KF2S ...
        - A1' .* KSF1 ...
        - A2' .* KSF2 ...
        + A1 .* KF1F1 .* A1' ...
        + A1 .* KF1F2 .* A2' ...
        + A2 .* KF2F1 .* A1' ...
        + A2 .* KF2F2 .* A2';

end

function Kc = local_outer_ring_cov(S, lambda, h, sqexp)

    r = sqrt(sum(S.^2, 2));

    F = S ./ r;

    KSS = sqexp(S, S, lambda);
    KSF = sqexp(S, F, lambda);
    KFS = KSF';
    KFF = sqexp(F, F, lambda);

    kff_diag = ones(size(r));

    A = diag(sqexp(S, F, h)) ./ kff_diag;

    Kc = KSS ...
        - A .* KFS ...
        - A' .* KSF ...
        + A .* KFF .* A';

end
