
clc; clear; close all;
addpath("/Users/yuema/Desktop/cGRF/Code/kernel/")

%% (Paper) Figure for A, f

figure;
tiledlayout(1, 5, 'TileSpacing', 'compact', 'Padding', 'compact');

% 1. One-sided boundary: {x1 = 0}
nexttile;
hold on; axis equal; box on;
fill([0 1 1 0], [0 0 1 1], [0.9 0.9 0.9], 'EdgeColor', 'none');
rectangle('Position', [0 0 1 1], 'EdgeColor', 'k', 'LineWidth', 1.5); 
plot([0 0], [0 1], 'r-', 'LineWidth', 2); % x1 = 0
title('One-sided');
xlim([-0.1 1.1]); ylim([-0.1 1.1]);
% Add arrows pointing left with label f^A
y_arrow = [0.2 0.5 0.8];
x_arrow = 0.2 * ones(size(y_arrow));
quiver(x_arrow, y_arrow, -0.2*ones(size(y_arrow)), zeros(size(y_arrow)), ...
       0, 'k', 'LineWidth', 2.5, 'MaxHeadSize', 4); % arrows left
text(0.25, 0.5, '$f$', 'Interpreter', 'latex', 'FontSize', 15);
xlabel('x_1')
ylabel('x_2')

xticks([0 0.5 1])
yticks([0 0.5 1])

% 2. Two-sided parallel: {x1 = 0 and x1 = 1}
nexttile;
hold on; axis equal; box on;
fill([0 1 1 0], [0 0 1 1], [0.9 0.9 0.9], 'EdgeColor', 'none');
rectangle('Position', [0 0 1 1], 'EdgeColor', 'k', 'LineWidth', 1.5);
plot([0 0], [0 1], 'r-', 'LineWidth', 2); % x1 = 0
plot([1 1], [0 1], 'r-', 'LineWidth', 2); % x1 = 1
title('Parallel');
xlim([-0.1 1.1]); ylim([-0.1 1.1]);
xlabel('x_1')
% Add arrows pointing left with label f^A
y_arrow = [0.2 0.5 0.8];
x_arrow = 0.2 * ones(size(y_arrow));
quiver(x_arrow, y_arrow, -0.2*ones(size(y_arrow)), zeros(size(y_arrow)), ...
       0, 'k', 'LineWidth', 2.5, 'MaxHeadSize', 4); % arrows left
text(0.25, 0.5, '$f_{1}$', 'Interpreter', 'latex', 'FontSize', 15);
% Add arrows pointing left with label f^A
y_arrow = [0.2 0.5 0.8];
x_arrow = 0.8 * ones(size(y_arrow));
quiver(x_arrow, y_arrow, 0.2*ones(size(y_arrow)), zeros(size(y_arrow)), ...
       0, 'k', 'LineWidth', 2.5, 'MaxHeadSize', 4); % arrows left
text(0.6, 0.5, '$f_{2}$', 'Interpreter', 'latex', 'FontSize', 15);
xlabel('x_1')

xticks([0 0.5 1])
yticks([0 0.5 1])

% % 3. Two-sided perpendicular: {x1 = 0 and x2 = 1}
% nexttile;
% hold on; axis equal; box on;
% fill([0 1 1 0], [0 0 1 1], [0.9 0.9 0.9], 'EdgeColor', 'none');
% rectangle('Position', [0 0 1 1], 'EdgeColor', 'k', 'LineWidth', 1.5);
% plot([0 0], [0 1], 'r-', 'LineWidth', 2); % x1 = 0
% plot([0 1], [0 0], 'r-', 'LineWidth', 2); % x2 = 0
% title('Perpendicular');
% xlim([-0.1 1.1]); ylim([-0.1 1.1]);
% xlabel('x_1')
% 
% xticks([0 0.5 1])
% yticks([0 0.5 1])

% 4. One-sided boundary of non-convex domain
nexttile;
hold on; axis equal; box on;
fill([0 1 1 0 0.5], [0 0 1 1 0.5], [0.9 0.9 0.9], 'EdgeColor', 'none');
plot([0 1 1 0], [0 0 1 1], 'black', 'LineWidth', 1.5);
plot([0 0.5 0], [0 0.5 1], 'r-', 'LineWidth', 2);
title('Angled');
xlim([-0.1 1.1]); ylim([-0.1 1.1]);
% Add arrows pointing left with label f^A
y_arrow = [0.2 0.5 0.8];
x_arrow = [0.4 0.7 0.4];
quiver(x_arrow, y_arrow, -0.2*ones(size(y_arrow)), zeros(size(y_arrow)), ...
       0, 'k', 'LineWidth', 2.5, 'MaxHeadSize', 4); % arrows left
text(0.75, 0.5, '$f$', 'Interpreter', 'latex', 'FontSize', 15);
xlabel('x_1')

xticks([0 0.5 1])
yticks([0 0.5 1])


% 5. Partial boundary: {x1 = -sqrt(1 - x2^2)} on unit disk
nexttile;
hold on; axis equal; box on;
theta = linspace(0, 2*pi, 300);
x = cos(theta); y = sin(theta);
fill(x, y, [0.9 0.9 0.9], 'EdgeColor', 'none');
plot(x, y, 'k', 'LineWidth', 1.5); % unit circle
x2_vals = linspace(-1, 1, 300);
x1_vals = -sqrt(1 - x2_vals.^2); % partial boundary
plot(x1_vals, x2_vals, 'r-', 'LineWidth', 2);
title('Partial');
xlim([-1.1 1.1]); ylim([-1.1 1.1]);
y_arrow = [-0.5 -0 0.5];
x_arrow = -sqrt(1 - y_arrow.^2) + 0.3;
quiver(x_arrow, y_arrow, -0.3*ones(size(y_arrow)), zeros(size(y_arrow)), ...
       0, 'k', 'LineWidth', 2.5, 'MaxHeadSize', 4); % larger arrows
text(-0.6, 0, '$f$', 'Interpreter', 'latex', 'FontSize', 15);
xlabel('x_1')

xticks([-1 0 1])
yticks([-1 0 1])

% % 6. Full boundary of unit disk
% nexttile;
% hold on; axis equal; box on;
% fill(x, y, [0.9 0.9 0.9], 'EdgeColor', 'none');
% plot(x, y, 'r-', 'LineWidth', 2); % full boundary
% title('Full');
% xlim([-1.1 1.1]); ylim([-1.1 1.1]);
% xlabel('x_1')
% 
% xticks([-1 0 1])
% yticks([-1 0 1])

% 7. Full boundary of ring-shaped domain
nexttile;
hold on; axis equal; box on;

% Draw annulus
theta = linspace(0, 2*pi, 400);

r_inner = 0.2;
r_outer = 1;

x_outer = r_outer * cos(theta);
y_outer = r_outer * sin(theta);

x_inner = r_inner * cos(theta);
y_inner = r_inner * sin(theta);

% Fill annulus
fill([x_outer fliplr(x_inner)], ...
     [y_outer fliplr(y_inner)], ...
     [0.9 0.9 0.9], 'EdgeColor', 'none');

% Draw boundaries
plot(x_outer, y_outer, 'r-', 'LineWidth', 2);
plot(x_inner, y_inner, 'r-', 'LineWidth', 2);

title('Concentric');
xlabel('x_1')

xlim([-1.1 1.1]);
ylim([-1.1 1.1]);

% Sample interior points for arrows
% theta_arrow = [pi/6, pi/3, pi/2];
theta_arrow = [pi*5/6];

x01 = 0.5 * cos(theta_arrow);
y01 = 0.5 * sin(theta_arrow);
x02 = 0.7 * cos(theta_arrow);
y02 = 0.7 * sin(theta_arrow);

% f_1 : projection to inner circle
x1 = r_inner * cos(theta_arrow);
y1 = r_inner * sin(theta_arrow);

quiver(x01, y01, ...
       x1 - x01, y1 - y01, ...
       0, 'k', 'LineWidth', 2.5, 'MaxHeadSize', 4);

text(-0.4, 0.4, '$f_{1}$', ...
    'Interpreter', 'latex', ...
    'FontSize', 15);

% f_2 : projection to outer circle
x2 = r_outer * cos(theta_arrow);
y2 = r_outer * sin(theta_arrow);

quiver(x02, y02, ...
       x2 - x02, y2 - y02, ...
       0, 'k', 'LineWidth', 2.5, 'MaxHeadSize', 4);

text(-0.8, 0.2, '$f_{2}$', ...
    'Interpreter', 'latex', ...
    'FontSize', 15);

hold off

xticks([-1 0 1])
yticks([-1 0 1])

% fontsize(15, "points")
fontsize(30, "points")

%% (Paper) 1d state and its linear transformation (SE)

num_samples = 100;
num_functions = 3;
a1 = 0; a2 = 1;
l = 0.5;
x = linspace(a1, a2, num_samples);

% m = zeros(num_samples*2, 1);
m = [cm(x, l, a1, a2, [1, 0, 0, 0], 0, 1); ...
    cdm(x, l, a1, a2, [1, 0, 0, 0], 0, 1)];
K = 1 * [cQQ(x, x, l, a1, a2, [1, 0, 0, 0]), ...
     cQR(x, x, l, a1, a2, [1, 0, 0, 0]); ...
     cQR(x, x, l, a1, a2, [1, 0, 0, 0])', ...
     cRR(x, x, l, a1, a2, [1, 0, 0, 0])];

f = mvnrnd(m, K, num_functions);

subplot(1,2,1)
% subplot(1,2,1)
for i=1:num_functions
    plot(x, f(i,1:num_samples), "LineWidth", 2)
    hold on;
end
hold off;
xlabel("x")
ylabel("u^A")
ylim([min(min(f)) max(max(f))])
axis square

subplot(1,2,2)
% subplot(1,2,2)
for i=1:num_functions
    plot(x, f(i,(num_samples+1):end), "LineWidth", 2)
    hold on;
end
% scatter(1, 1, 50, "filled", "red")
hold off;
xlabel("x")
ylabel("u_x^A")
ylim([min(min(f)) max(max(f))])
axis square

% subplot(1,3,3)
% for i=1:num_functions
%     plot(x, f(i,1:num_samples)+f(i,(num_samples+1):end), "LineWidth", 2)
%     hold on;
% end
% % scatter(0, 0, 50, "filled", "red")
% hold off;
% xlabel("x")
% ylabel("u^A+u_x^A")
% ylim([min(min(f)) max(max(f))])
% fontsize(15, "points")
% axis square

%% (Paper) GB & 2d Continuous versus Pointwise Constraining

%%% 1d

num_samples = 100;
num_functions = 3;
a1 = 0; a2 = 1;
nu = 2.5; alpha = 20; l = 2.5;
x = linspace(a1, a2, num_samples);

% m = zeros(num_samples, 1);
m = cm(x, l, a1, a2, [0, 1, 0, 1], 0, 0);
K = (1/alpha) * cQQ_matern(x, x, l, a1, a2, [0, 1, 0, 1]);
f_1d_cGRF = mvnrnd(m, K, num_functions);

x1 = [a1, a2];
f1 = zeros(2, 1);
m1 = zeros(2, 1);
K11 = (1/alpha) * QQ_matern(x1, x1, nu);
K12 = (1/alpha) * QQ_matern(x1, x, nu);
K22 = (1/alpha) * QQ_matern(x, x, nu);
[m_1d_cond, K_1d_cond, ~] = GP(f1, m1, m, K11, K12, K22, eps);
f_1d_cond = mvnrnd(m_1d_cond, K_1d_cond, num_functions);

%%% Gaussian bridge

figure;
tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile
for i=1:num_functions
    plot(x, f_1d_cond(i,:), "LineWidth", 1.5)
    hold on;
end
xlabel("x")
ylabel("u^{\{0,1\}}_{GB}")
% yticks([-0.5 0 0.5])
ylim([min(min(f_1d_cond)) max(max(f_1d_cond))])
plot([0,1], [0,0], "r.", MarkerSize=20)
hold off;
fontsize(20,"points")
axis square

%%% continuous constraining

figure;
tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile
for i=1:num_functions
    plot(x, f_1d_cGRF(i,:), "LineWidth", 1.5)
    hold on;
end
xlabel("x")
ylabel("u^{\{0,1\}}")
% yticks([-0.5 0 0.5])
ylim([min(min(f_1d_cGRF)) max(max(f_1d_cGRF))])
% plot([0,1], [0,0], "r.", MarkerSize=20)
hold off;
fontsize(20,"points")
axis square

%%% 2d 

num_samples = 30;
a1 = 0; a2 = 1;
nu = 2.5;
x = linspace(0, 1, num_samples);
y = linspace(0, 1, num_samples);
xx = repmat(x, 1, length(y))';
yy = repelem(y, 1, length(x))';

num_samples1 = 9;
a = linspace(a1, a2, num_samples1);
b = zeros(1, length(a)) + a1;
c = zeros(1, length(a)) + a2;
xx1 = cat(2, a, a, b, c)';
yy1 = cat(2, b, c, a, a)';
f1 = zeros(length(xx1), 1);
m1 = zeros(length(xx1), 1);


%%% continuous constraining

m = zeros(length(xx), 1);
% sequential cGRF construction over two coordinates; product covariance
% structure preserved for rectangular domains.
K = cQQ_matern(xx, xx, nu, a1, a2, [0, 1, 0, 1]) .* ...
    cQQ_matern(yy, yy, nu, a1, a2, [0, 1, 0, 1]);
f_cGRF = mvnrnd(m, K, num_functions);

%%% pointwise constraining (conditioning)


K11 = QQ_matern(xx1, xx1, nu) .* QQ_matern(yy1, yy1, nu);
K12 = QQ_matern(xx1, xx, nu) .* QQ_matern(yy1, yy, nu);
K22 = QQ_matern(xx, xx, nu) .* QQ_matern(yy, yy, nu);
[m_cond, K_cond, ~] = GP(f1, m1, m, K11, K12, K22, eps);

f_cond = mvnrnd(m_cond, K_cond, num_functions);


figure;
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile
surf(y, x, reshape(f_cGRF(1,:), num_samples, num_samples), "FaceAlpha", 0.8)
hold on
k = boundary(xx, yy);
plot3(yy(k), xx(k), f_cGRF(1,k), 'r-', "LineWidth", 3)
hold off
xlabel("x_1")
ylabel("x_2")
zlabel("u^A")
% zticks([-0.1 0 0.1])
zlim([min(min(f_cGRF)) max(max(f_cGRF))])
fontsize(20,"points")
axis square

nexttile
surf(y, x, reshape(f_cGRF(2,:), num_samples, num_samples), "FaceAlpha", 0.8)
hold on
k = boundary(xx, yy);
plot3(yy(k), xx(k), f_cGRF(2,k), 'r-', "LineWidth", 3)
hold off
xlabel("x_1")
ylabel("x_2")
% zlabel("u(x_1,x_2)")
% zticks([-0.2 0 0.2])
zlim([min(min(f_cGRF)) max(max(f_cGRF))])
fontsize(20,"points")
axis square

nexttile
surf(y, x, reshape(f_cGRF(3,:), num_samples, num_samples), "FaceAlpha", 0.8)
hold on
k = boundary(xx, yy);
plot3(yy(k), xx(k), f_cGRF(3,k), 'r-', "LineWidth", 3)
hold off
xlabel("x_1")
ylabel("x_2")
% zlabel("u(x_1,x_2)")
% zticks([-0.2 0 0.2])
zlim([min(min(f_cGRF)) max(max(f_cGRF))])
fontsize(20,"points")
axis square


figure;
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile
surf(x, y, reshape(f_cond(1,:), num_samples, num_samples), "FaceAlpha", 0.8)
hold on
plot3(xx1, yy1, f1, "r.", MarkerSize=18)
k = boundary(xx, yy);
plot3(yy(k), xx(k), f_cond(1,k), 'r-', "LineWidth", 3)
hold off
xlabel("x_1")
ylabel("x_2")
zlabel("u")
% zticks([-0.1 0 0.1])
zlim([min(min(f_cond)) max(max(f_cond))])
fontsize(20,"points")
axis square

nexttile
surf(x, y, reshape(f_cond(2,:), num_samples, num_samples), "FaceAlpha", 0.8)
hold on
plot3(xx1, yy1, f1, "r.", MarkerSize=18)
k = boundary(xx, yy);
plot3(yy(k), xx(k), f_cond(2,k), 'r-', "LineWidth", 3)
hold off
xlabel("x_1")
ylabel("x_2")
% zlabel("u(x_1,x_2)")
% zticks([-0.2 0 0.2])
zlim([min(min(f_cond)) max(max(f_cond))])
fontsize(20,"points")
axis square

nexttile
surf(x, y, reshape(f_cond(3,:), num_samples, num_samples), "FaceAlpha", 0.8)
hold on
plot3(xx1, yy1, f1, "r.", MarkerSize=18)
k = boundary(xx, yy);
plot3(yy(k), xx(k), f_cond(3,k), 'r-', "LineWidth", 3)
hold off
xlabel("x_1")
ylabel("x_2")
% zlabel("u(x_1,x_2)")
% zticks([-0.2 0 0.2])
zlim([min(min(f_cond)) max(max(f_cond))])
fontsize(20,"points")
axis square


%% (Paper) Introduction (Robin)

num_samples = 20;
num_functions = 3;
a1 = 0; a2 = 1;
l = 0.5;

x = linspace(0, 1, num_samples);
t = linspace(0, 1, num_samples);
xx = repmat(x, 1, length(t))';
tt = repelem(t, 1, length(x))';
coeffs = [1, 1, 1, 0];
m = zeros(length(xx)*2, 1);
K = [cQQ(tt, tt, l, a1, a2, coeffs).*QQ(xx, xx, l), ...
     cQR(tt, tt, l, a1, a2, coeffs).*QQ(xx, xx, l); ...
     cQR(tt, tt, l, a1, a2, coeffs)'.*QQ(xx, xx, l), ...
     cRR(tt, tt, l, a1, a2, coeffs).*QQ(xx, xx, l)];
f_robin = mvnrnd(m, K, num_functions);

figure;
tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

plot_simu_model(t, x, f_robin(1,:), true)
plot_simu_model(t, x, f_robin(2,:), true) 

%% Disk-shaped (Not included in paper)

theta = linspace(0, 2*pi, 60);
seq = 3*cumsum(1:20); seq = seq - seq(1);
rho = 1-seq / seq(end); 
[Theta, Rho] = meshgrid(theta, rho);
l = 0.5; 
x_ = Rho .* cos(Theta);
y_ = Rho .* sin(Theta);
xx = x_(:);
yy = y_(:);
m = cm(xx, l, -sqrt(1-yy.^2), sqrt(1-yy.^2), [0, 1, 0, 1], (-sqrt(1-yy.^2)+yy).^2, (sqrt(1-yy.^2)+yy).^2);
K = cQQ(xx, xx, l, -sqrt(1-yy.^2), sqrt(1-yy.^2), [0, 1, 0, 1]);
f_disk = mvnrnd(m, K, num_functions);

figure;
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile
surf(x_, y_, reshape(f_disk(1,:), size(x_)), 'FaceColor', 'interp', 'FaceAlpha', 0.8)
hold on
theta=0:0.01:2*pi;
v=null([0,0,0]);
points=repmat([0,0,0]',1,size(theta,2))+1*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
plot3(points(1,:),points(2,:),(points(1,:)+points(2,:)).^2,'-', 'Color', 'red', "LineWidth", 3);
% idx = points(1,:) < 0;
% plot3(points(1,idx),points(2,idx),(points(1,idx)+points(2,idx)).^2,'-',"Color",'red', "LineWidth", 2.5);
xlabel("x_1") 
ylabel("x_2") 
zlabel("u")
fontsize(15, "points")
hold off
axis square

nexttile
surf(x_, y_, reshape(f_disk(2,:), size(x_)), 'FaceColor', 'interp', 'FaceAlpha', 0.8)
hold on
theta=0:0.01:2*pi;
v=null([0,0,0]);
points=repmat([0,0,0]',1,size(theta,2))+1*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
plot3(points(1,:),points(2,:),(points(1,:)+points(2,:)).^2,'-', 'Color', 'red', "LineWidth", 3);
xlabel("x_1") 
ylabel("x_2") 
zlabel("u")
fontsize(15, "points")
hold off
axis square


%% Unit triangular domain (Paper)

num_samples = 20;
nu = 2.5;

x = linspace(0, 1, num_samples);
y = flip(1-x);
xx = repmat(x, 1, length(y))';
yy = repelem(y, 1, length(x))';
idx = -xx+1 >= yy;

m = zeros(length(xx(idx)), 1);
K = cQQ_matern(xx(idx), xx(idx), nu, zeros(size(yy(idx))), 1-yy(idx), [0, 1, 0, 1]);
K(isnan(K)) = 0; % manually specify K when a1=a2, which has to be mutually exclusive

figure;
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for i=1:2

nexttile
f = mvnrnd(m, K, 1);
zz = zeros(length(xx), 1);
zz(idx) = f;
zz(idx == 0) = NaN;
T = delaunay(xx(idx), yy(idx));
trisurf(triangulation(T, xx(idx), yy(idx), zz(idx)), 'EdgeColor', 'black', 'FaceColor', 'interp', 'FaceAlpha', 0.8)
hold on
% line(x, zeros(size(x)), zeros(size(x)), 'Color', 'red', 'LineWidth', 2.5);
line(zeros(size(x)), x, zeros(size(x)), 'Color', 'red', 'LineWidth', 3);
line(x, 1-x, zeros(size(x)), 'Color', 'red', 'LineWidth', 3);
xlabel("x_1") 
ylabel("x_2") 
if i == 1
zlabel("u^A")
end
fontsize(15, "points")
hold off
axis square
end

%% 2d product covariance with periodic (Paper)

num_samples = 30;
a1 = 0; a2 = 1;
nu = 2.5;
x = linspace(0, 1, num_samples);
y = linspace(0, 1, num_samples);
xx = repmat(x, 1, length(y))';
yy = repelem(y, 1, length(x))';

m = zeros(length(xx), 1);
K = cQQ_matern(xx, xx, nu, a1, a2, [0, 1, 0, 1]) .* ...
    (0.1 * QQ_periodic(yy, yy, 1, 1));

figure;
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for i=1:2

nexttile
f = mvnrnd(m, K, 1);
surf(y, x, reshape(f,num_samples, num_samples), 'FaceAlpha', 0.8)
hold on
k = boundary(xx, yy);
plot3(yy(k), xx(k), f(k), 'r-', "LineWidth", 3)
hold off
xlabel("x_1")
ylabel("x_2")
if i == 1
zlabel("u^A")
end
fontsize(15,"points")
axis square

end


%% Functions

function plot_simu_model(t,x,f, deriv)

n = length(t);

nexttile
for i=1:size(f,1)
u = reshape(f(i, 1:n*n), n, n);
if i~=1
    surf(t, x, u, 'FaceColor', 'interp', 'FaceAlpha', 0.8);
    continue
end
surf(t, x, u, 'FaceColor', 'interp', 'FaceAlpha', 0.8);
hold on
end
hold off
xlabel("x_1")
ylabel("x_2")
zlabel("u^A")
axis square

nexttile
for i=1:size(f,1)
ut = reshape(f(i, (1+n*n):end), n, n);
if i~=1
    surf(t, x, ut, 'FaceColor', 'interp', 'FaceAlpha', 0.8);
    continue
end
if deriv == true
    surf(t, x, ut, 'FaceColor', 'interp', 'FaceAlpha', 0.8);
    hold on
else
    surf(t, x, zeros(size(ut)), 'FaceColor', 'interp', 'FaceAlpha', 0.5);
    hold on
end
plot3(zeros(size(t))+t(end), x, ut(:,end), 'red', 'LineWidth', 3)
end
hold off
xlabel("x_1")
ylabel("x_2")
zlabel("u_{x_1}^A")
axis square

nexttile
for i=1:size(f,1)
u = reshape(f(i, 1:n*n), n, n);
ut = reshape(f(i, (1+n*n):end), n, n);
if i~=1
    surf(t, x, u+ut, 'FaceColor', 'interp', 'FaceAlpha', 0.8);
    continue
end
surf(t, x, u+ut, 'FaceColor', 'interp', 'FaceAlpha', 0.8);
hold on
plot3(zeros(size(t))+t(1), x, u(:,1)+ut(:,1), 'red', 'LineWidth', 3)
end
hold off
xlabel("x_1")
ylabel("x_2")
zlabel("u^A+u_{x_1}^A")
axis square

fontsize(15, "points")

end

