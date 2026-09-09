%% 
clear; clc;
addpath("/Users/yuema/Desktop/GP/Code/kernel/")
addpath("/Users/yuema/Desktop/GP/Code/Tensile-smoothing/")

%%

A = readmatrix("tensileDog.csv");

idx_train = [1,4,7,10];
idx_val = [2,5,8,11];
idx_test = [3,6,9];

%% Visualization of mean

n = 100;
y = linspace(-10, 10, n);
t = linspace(0, 1, n);
yy = repmat(y, 1, length(t))';
tt = repelem(t, 1, length(y))';
m = cm(yy, 1, -10, 10, [0,1,0,1], zeros(size(tt)), tt*0.05);
surf(t,y,reshape(m,n,n))

%% Optimizer set-up

% options = optimoptions('fmincon', 'StepTolerance', 1e-6);
options = [];

%% Setup

path = "result_fixed_design/";
ratio = "sparse";
noise = 0.01;
params_initial =  [1, 1, 1, 1, 0.001];
metric = "mse";
rep = 100;

MSE_uncon = [];
MSE_con = [];
params_uncon = [];
params_con = [];

%%
count = 0;
for i = 1:rep
    
    data_train = data_loader(A, idx_train, false, ratio, noise);
    data_val = data_loader(A, idx_val, false, ratio, noise);
    data_test = data_loader(A, idx_test, true, ratio, noise); % ratio and noise unrelated if testing==true
    
    % uncontrained
    [params_best_uncon, ~] = fmincon(@(x) krigging(@cov_uncon, @mean_uncon, @GP, data_train, data_val, x, metric, true), ...
        params_initial, [], [], [], [], [0.05, 1, 0.1, 1, 0], [10, 1, 20, 1, 0.1], [], options);
    mse_uncon = krigging(@cov_uncon, @mean_uncon, @GP, data_train, data_test, params_best_uncon, metric, true);
    
    % contrained
    [params_best_con, ~] = fmincon(@(x) krigging(@cov_con, @mean_con, @GP, data_train, data_val, x, metric, true), ...
        params_initial, [], [], [], [], [0.05, 1, 0.1, 1, 0], [10, 1, 20, 1, 0.1], [], options);
    mse_con = krigging(@cov_con, @mean_con, @GP, data_train, data_test, params_best_con, metric, true);

    MSE_uncon = [MSE_uncon; mse_uncon];
    MSE_con = [MSE_con; mse_con];
    params_uncon = [params_uncon; params_best_uncon];
    params_con = [params_con; params_best_con];

    display(i)
    count = count + (mse_uncon > mse_con);
    display(count)
    
    % save(path+"MSE_uncon_ratio"+num2str(ratio)+"noise"+num2str(noise)+".mat", "MSE_uncon")
    % save(path+"MSE_con_ratio"+num2str(ratio)+"noise"+num2str(noise)+".mat", "MSE_con")
    % save(path+"param_uncon_ratio"+num2str(ratio)+"noise"+num2str(noise)+".mat", "params_uncon")
    % save(path+"param_con_ratio"+num2str(ratio)+"noise"+num2str(noise)+".mat", "params_con")

end

% save(path+"MSE_uncon_ratio"+num2str(ratio)+"noise"+num2str(noise)+".mat", "MSE_uncon")
% save(path+"MSE_con_ratio"+num2str(ratio)+"noise"+num2str(noise)+".mat", "MSE_con")
% save(path+"param_uncon_ratio"+num2str(ratio)+"noise"+num2str(noise)+".mat", "params_uncon")
% save(path+"param_con_ratio"+num2str(ratio)+"noise"+num2str(noise)+".mat", "params_con")

%% Visualization of one example setting: result

path = "results/";
noise = 0.0001;
ratio = "medium";
load(path+"MSE_uncon_ratio"+num2str(ratio)+"noise"+num2str(noise)+".mat", "MSE_uncon")
load(path+"MSE_con_ratio"+num2str(ratio)+"noise"+num2str(noise)+".mat", "MSE_con")

tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax1 = nexttile;
bx = boxplot([MSE_con, MSE_uncon], 'OutlierSize', 2);
set(bx, 'LineWidth', 1.5)
blue   = [0 0.4470 0.7410];
orange = [0.8500 0.3250 0.0980];
set(findobj(bx,'Tag','Median'), 'Color', orange, 'LineWidth', 1.5);
set(findobj(bx,'Tag','Outliers'), 'MarkerEdgeColor', orange);
set(findobj(bx,'Tag','Box'), 'Color', blue, 'LineWidth', 1.5);
set(findobj(bx,'Tag','Whisker'), 'Color', blue);
set(findobj(bx,'Tag','Upper Whisker'), 'Color', blue);
set(findobj(bx,'Tag','Lower Whisker'), 'Color', blue);
set(findobj(bx,'Tag','Upper Adjacent Value'), 'Color', blue);
set(findobj(bx,'Tag','Lower Adjacent Value'), 'Color', blue);
xticklabels(["cGRF", "GRF"])
title("\sigma^2="+num2str(noise)+"^2")
ylabel("log(MSE)");
ylim([-20 -8])
fontsize(15,"points")

%% Visualization of one example setting: design

tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax2 = nexttile;
hold on
plot([-1, 1, 1, 0.5, 0.5, 1, 1, -1, -1, -0.5, -0.5, -1, -1], ...
     [-10, -10, -7, -5, 5, 7, 10, 10, 7, 5, -5, -7, -10], "LineWidth", 1.5)
load("design_"+ratio+".mat", "points_new")
scatter(points_new(:,1), points_new(:,2), 15, "filled");
hold off
axis equal; 
xlim([-1.3 1.3]); ylim([-10.3 10.3]);
title("Design")
yyaxis right
set(gca, 'YTick', [], 'YColor', 'k')
h = ylabel(ratio);
set(h,'Rotation',-90,'VerticalAlignment','bottom','HorizontalAlignment','center');
fontsize(15,"points")

%% Visualization of all settings: result and design

path = "results/";

noises = [0.0001, 0.001, 0.01];
ratios = ["dense", "medium", "low", "sparse"];
t = tiledlayout(4, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

for i=1:length(ratios)
    ratio = ratios(i);

    for j=1:length(noises)
        noise = noises(j);
        load(path+"MSE_uncon_ratio"+num2str(ratio)+"noise"+num2str(noise)+".mat", "MSE_uncon")
        load(path+"MSE_con_ratio"+num2str(ratio)+"noise"+num2str(noise)+".mat", "MSE_con")
        
        nexttile
        bx = boxplot([MSE_con, MSE_uncon], 'OutlierSize', 2);
        set(bx, 'LineWidth', 1.5)
        % --- Adjust components ---
        % Default MATLAB colors
        blue   = [0 0.4470 0.7410];
        orange = [0.8500 0.3250 0.0980];

        % Median line
        set(findobj(bx,'Tag','Median'), 'Color', orange, 'LineWidth', 1.5);

        % Outliers
        set(findobj(bx,'Tag','Outliers'), 'MarkerEdgeColor', orange);

        % Box
        set(findobj(bx,'Tag','Box'), 'Color', blue, 'LineWidth', 1.5);

        % Whiskers and caps (keep same as box, optional)
        set(findobj(bx,'Tag','Whisker'), 'Color', blue);
        set(findobj(bx,'Tag','Upper Whisker'), 'Color', blue);
        set(findobj(bx,'Tag','Lower Whisker'), 'Color', blue);
        set(findobj(bx,'Tag','Upper Adjacent Value'), 'Color', blue);
        set(findobj(bx,'Tag','Lower Adjacent Value'), 'Color', blue);


        if i ~= 4
            set(gca, 'XTick', [])
        else
            xticklabels(["cGRF", "GRF"])
        end
        if i == 1
            title("\sigma^2="+num2str(noise)+"^2")
        end
        if j == 1
            ylabel("log(MSE)");
        else
            set(gca, 'YTick', [])
        end
        ylim([-20 -8])
    end
    nexttile
    plot([-1, 1, 1, 0.5, 0.5, 1, 1, -1, -1, -0.5, -0.5, -1, -1], ...
         [-10, -10, -7, -5, 5, 7, 10, 10, 7, 5, -5, -7, -10], "LineWidth", 1.5)
    hold on
    load("design_"+ratio+".mat", "points_new")
    scatter(points_new(:,1), points_new(:,2), 15, "filled");
    hold off
    pbaspect([0.5 1 1])
    xlim([-1.1 1.1]); ylim([-10.5 10.5]);
    if ratio ~= "sparse"
        set(gca, 'XTick', [])
    end
    if ratio == "dense"
        title("Design")
    end
    yyaxis right
    set(gca, 'YTick', [], 'YColor', 'k')
    h = ylabel(ratio);
    set(h,'Rotation',-90,'VerticalAlignment','bottom','HorizontalAlignment','center');
   
end
fontsize(15,"points")


%% Function

function data = data_loader(metadata, times, testing, ratio, noise)
    n = size(metadata, 1);
    data = zeros(n*length(times), 4);
    for t = 1:length(times)
        data(((t-1)*n+1):(t*n), :) = metadata(:, [2+(times(t)-1)*4+4,1,2,2+(times(t)-1)*4+2]);
    end

    if testing == false
        load("design_"+ratio+".mat", "points_new")
        idx = ismember(data(:,2:3), points_new, "rows");
        data = data(idx, :);
        data(:, 4) = data(:, 4) + noise*randn(size(data(:, 4)));
    end
    
end

function K = cov_uncon(data_j, data_k, params)
    lt = params(1); ly = params(3); alpha = params(4);
    K = (1/alpha) * QQ(data_j(:,1), data_k(:,1), lt) .* QQ(data_j(:,3), data_k(:,3), ly);
end

function mu = mean_uncon(data_i, params)
    mu = zeros(size(data_i,1), 1);
end

function K = cov_con(data_j, data_k, params)
    lt = params(1); ly = params(3); alpha = params(4);
    K = (1/alpha) * QQ(data_j(:,1), data_k(:,1), lt) .* cQQ(data_j(:,3), data_k(:,3), ly, -10, 10, [0, 1, 0, 1]);
end

function mu = mean_con(data_i, params)
    lt = params(1); ly = params(3); alpha = params(4);
    mu = cm(data_i(:,3), ly, -10, 10, [0,1,0,1], zeros(size(data_i(:,1))), data_i(:,1)*0.05);
end

function result = krigging(fcov, fmu, GP, data, data_test_i, params_best, metric, training)

    mu = fmu(data, params_best);
    K = fcov(data, data, params_best);

    K_ = fcov(data_test_i, data_test_i, params_best);
    mu_ = fmu(data_test_i, params_best);

    C = fcov(data, data_test_i, params_best);

    [f_mu, f_K, ~] = GP(data(:,4), mu, mu_, K, C, K_, params_best(end));
    f_se = sqrt(diag(f_K));
    
    if metric == "mse"
        se_best = (data_test_i(:,4) - f_mu).^2;
        mse_best = log(mean(se_best));
    elseif metric == "nlpd"
        se_best = -(log(normpdf(data_test_i(:,4), f_mu, f_se)));
        mse_best = mean(se_best);
    end
    
    if training
        result = mse_best;
    else
        result = {mse_best, se_best, f_mu, f_se};
    end
end

