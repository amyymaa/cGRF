
clear; clc; close all

%% Load boundary vertices

load("/Users/yuema/Desktop/cGRF/Code/kernel/boundary_dense.mat")
load("/Users/yuema/Desktop/cGRF/Code/kernel/boundary_sparse.mat")

%% Interior grid over the non-convex polygon

dx = 0.005;
lon = sparse(:,1);
lat = sparse(:,2);

minx = min(lon); maxx = max(lon);
miny = min(lat); maxy = max(lat);

[xg, yg] = meshgrid(minx:dx:maxx, miny:dx:maxy);
inside = inpolygon(xg, yg, lon, lat);

X_in = xg(inside);
Y_in = yg(inside);
S = [X_in, Y_in];
n = size(S,1);
m = zeros(n,1);


%% Visualization of linear approximation

plot(dense(:,1), dense(:,2), LineWidth=2);
hold on
plot(sparse(:,1), sparse(:,2), LineWidth=1.5);

for i = 1:length(sparse)
    label = sprintf('(%.2f, %.2f)', sparse(i,1), sparse(i,2));
    text(sparse(i,1), sparse(i,2), label, ...
        'VerticalAlignment','bottom', ...
        'HorizontalAlignment','left');
end
hold off
axis square
xlabel("x_1") 
ylabel("x_2") 
yticks([52.8 52.9 53.0 53.1 53.2])
fontsize(15, "points")
set(gcf, 'GraphicsSmoothing', 'on')

%% Boundary functions

bd_top    = f_top();
bd_bottom = f_bottom();
bd_left   = f_left();
bd_right  = f_right();

%% Base covariance and cGRF covariance

lambda = 0.2;   % base squared-exponential length-scale
h = 0.2;           % decay length-scale for weights

% cov = @(A,B,ell) exp(-pdist2(A,B).^2 ./ (2*ell^2));
% cov = @(A,B,ell) (1 + (sqrt(3) * pdist2(A,B) ./ ell)) .* exp(-(sqrt(3) * pdist2(A,B) ./ ell));
cov = @(A,B,ell) (1 + (sqrt(5) * pdist2(A,B) ./ ell) + (sqrt(5) * pdist2(A,B) ./ ell).^2/3) .* exp(-(sqrt(5) * pdist2(A,B) ./ ell));

K = constr_cov(S, lambda, h, cov, bd_top, bd_bottom, bd_left, bd_right);

K(isnan(K)) = 0;
K = (K + K') / 2;


%% Draw three samples

rng(9)
f = mvnrnd(m, K, 3);

% Triangulation restricted to polygon
T = delaunay(X_in, Y_in);
cx = mean(X_in(T), 2);
cy = mean(Y_in(T), 2);
insideTri = inpolygon(cx, cy, sparse(:,1), sparse(:,2));
T2 = T(insideTri,:);

clim_min = min(f(:));
clim_max = max(f(:));


figure;
t = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
ax_all = gobjects(1,3);

for i = 1:3
    ax = nexttile;
    ax_all(i) = ax;

    trisurf(triangulation(T2, X_in, Y_in, f(i,:)'), ...
        'EdgeColor', 'none', ...
        'FaceColor', 'interp', ...
        'FaceAlpha', 0.8);

    view(2);
    shading interp;
    axis square tight;
    clim([clim_min, clim_max]);
    xlim([-7.8, -6.8]);
    ylim([52.75, 53.25]);

    hold on;
    plot(sparse(:,1), sparse(:,2), 'r-', 'LineWidth', 2.5);
    hold off;

    % title(sprintf('Sample %d', i));
    % xlabel('$x_1$', 'Interpreter', 'latex');
    xlabel("x_1")
    if i == 1
        % ylabel('$x_2$', 'Interpreter', 'latex');
        ylabel("x_2")
    else
        yticks([]);
    end

    if i == 3
        colorbar;
    end
end

fontsize(15, 'points')

%% Functions

function Kc = constr_cov(S, lambda, h, sqexp, bd_top, bd_bottom, bd_left, bd_right)

    n = size(S,1);
    x = S(:,1);
    y = S(:,2);

    % Four projections for each interior point s = (x,y)
    F_T = [x, bd_top(x)];
    F_B = [x, bd_bottom(x)];
    F_L = [bd_left(y), y];
    F_R = [bd_right(y), y];
    F = cat(3, F_T, F_B, F_L, F_R);
    Fcell = {F_T, F_B, F_L, F_R};

    % Precompute weights
    A = zeros(n,4);
    for i = 1:n
        Fi = squeeze(F(i,:,:))';
        KF_i = sqexp(Fi, Fi, h);
        k_i  = sqexp(S(i,:), Fi, h)';

        a_i = KF_i \ k_i;
        A(i,:) = a_i';
    end

    % Precompute covariance blocks for constrained formulation
    KSS = sqexp(S, S, lambda);
    KSF = cell(1,4);
    KFS = cell(1,4);
    KFF = cell(4,4);
    for j = 1:4
        KSF{j} = sqexp(S, Fcell{j}, lambda);
        KFS{j} = KSF{j}';                     
        for k = 1:4
            KFF{j,k} = sqexp(Fcell{j}, Fcell{k}, lambda);
        end
    end
    
    % Constrained formulation
    Kc = KSS;
    for j = 1:4
        aj = A(:,j);
        Kc = Kc - aj .* KFS{j} - aj' .* KSF{j};   
    end
    for j = 1:4
        aj = A(:,j);
        for k = 1:4
            ak = A(:,k);
            Kc = Kc + (aj * ak') .* KFF{j,k};
        end
    end

    Kc = (Kc + Kc') / 2;
end

function bd_top = f_top()
    p1 = [-7.63, 53.17];
    p2 = [-7.07, 53.17];
    f12 = linear(p1, p2);
    bd_top = @(x) f12(x);
end

function bd_bottom = f_bottom()
    p1 = [-7.67, 52.78];
    p2 = [-7.20, 52.89];
    p3 = [-6.97, 52.81];
    f12 = linear(p1, p2);
    f23 = linear(p2, p3);
    bd_bottom = @(x) f12(x) .* (x <= -7.20) + f23(x) .* (x > -7.20);
end

function bd_left = f_left()
    p1 = [53.17, -7.63];
    p2 = [53.10, -7.59];
    p3 = [53.01, -7.71];
    p4 = [52.93, -7.66];
    p5 = [52.86, -7.73];
    p6 = [52.78, -7.67];
    f12 = linear(p1, p2);
    f23 = linear(p2, p3);
    f34 = linear(p3, p4);
    f45 = linear(p4, p5);
    f56 = linear(p5, p6);
    bd_left = @(y) f12(y) .* (y >= 53.10) ...
        + f23(y) .* ((y < 53.10) & (y >= 53.01)) ...
        + f34(y) .* ((y < 53.01) & (y >= 52.93)) ...
        + f45(y) .* ((y < 52.93) & (y >= 52.86)) ...
        + f56(y) .* (y < 52.86);
end

function bd_right = f_right()
    p1 = [53.17, -7.07];
    p2 = [53.00, -7.09];
    p3 = [52.95, -6.95];
    p4 = [52.81, -6.97];
    f12 = linear(p1, p2);
    f23 = linear(p2, p3);
    f34 = linear(p3, p4);
    bd_right = @(y) f12(y) .* (y >= 53.00) ...
        + f23(y) .* ((y < 53.00) & (y >= 52.95)) ...
        + f34(y) .* (y < 52.95);
end

function f = linear(p1, p2)
    x1 = p1(1); y1 = p1(2);
    x2 = p2(1); y2 = p2(2);
    m = (y2 - y1) / (x2 - x1);
    b = y1 - m * x1;
    f = @(x) m * x + b;
end
