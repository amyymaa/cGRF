%%

clear; clc;
addpath("/Users/yuema/Desktop/cGRF/Code/kernel/")
addpath("/Users/yuema/Desktop/cGRF/Code/Tensile-smoothing-real/data/")
addpath("/Users/yuema/Desktop/cGRF/Code/Tensile-smoothing-real/data1/")
dataFolder = "/Users/yuema/Desktop/cGRF/Code/Tensile-smoothing-real/data1/";


%% Design

preproc.use_preprocessing = false;

ref_fname = fullfile(dataFolder, sprintf('Poly-%04d_1.csv', 100));
[Xref, Yref, Vref, info_ref] = read_and_preprocess_file(ref_fname, preproc);
Nref = numel(Xref);

design_grid_sizes = [9, 20];

designs = struct([]);
threshold = 1; % to use knn with k = 1 and distance <= threshold
nx = 9; ny = 27;
idx_design = equal_spaced_design_indices(Xref, Yref, nx, ny, threshold);
designs(1).grid_size = [nx, ny];
designs(1).idx = idx_design(:);
designs(1).n_points = numel(idx_design);

xB = [-9, 9, 9, 2.8, 2.8, 9, 9, -9, -9, -2.8, -2.8, -9, -9];
yB = [-37.5, -37.5, -32, -19, 17, 30, 35, 35, 30, 17, -19, -32, -37.5];

plot(xB, yB, 'b-', 'LineWidth', 1.5, 'Color', [0.0000, 0.4470, 0.7410]); hold on
scatter(Xref(designs(1).idx), Yref(designs(1).idx), 15, [0.8500, 0.3250, 0.0980], 'filled');
axis equal

fontsize(15, "points")

%%

preproc.use_preprocessing = false;

times_all = 100:10:200;
data_all = load_time_block(dataFolder, times_all, designs(1).idx, true, preproc);

n_iter = 100;
train_frac = 0.50;
val_frac   = 0.25;
test_frac  = 0.25;

rng(1);

% Modify measurement error variancelast entry from 1e-4 to 1e-8
theta0 = [10, 10, 10, 1, 1e-7];

lb     = [1,  1,  1,  1, 1e-7];
ub     = [50, 50, 50, 1, 1e-7];

options = optimoptions('fmincon', 'Display', 'off');

results_con = struct([]);
results_uncon = struct([]);

N = size(data_all, 1);

count = 0;

for r = 1:n_iter
    fprintf("cGRF outperforms GRF %d / %d\n", count, r - 1);
    fprintf("Iteration %d / %d\n", r, n_iter);

    idx = randperm(N);

    n_train = round(train_frac * N);
    n_val   = round(val_frac * N);

    idx_train = idx(1:n_train);
    idx_val   = idx(n_train+1:n_train+n_val);
    idx_test  = idx(n_train+n_val+1:end);

    data_train = data_all(idx_train, :);
    data_val   = data_all(idx_val, :);
    data_test  = data_all(idx_test, :);

    %% cGRF vs GRF
    
    obj_uncon = @(theta) krigging(@cov_uncon, @mean_uncon, @GP, data_train, data_val, theta, "mse", true);
    obj_con = @(theta) krigging(@cov_con, @mean_con, @GP, data_train, data_val, theta, "mse", true);
    [params_best_uncon, ~] = fmincon(obj_uncon, theta0, [], [], [], [], lb, ub);
    [params_best_con, ~] = fmincon(obj_con, theta0, [], [], [], [], lb, ub);

    mse_uncon = krigging(@cov_uncon, @mean_uncon, @GP, data_train, data_test, params_best_uncon, "mse", true);
    nlpd_uncon = krigging(@cov_uncon, @mean_uncon, @GP, data_train, data_test, params_best_uncon, "nlpd", true);
    mse_con = krigging(@cov_con, @mean_con, @GP, data_train, data_test, params_best_con, "mse", true);
    nlpd_con = krigging(@cov_con, @mean_con, @GP, data_train, data_test, params_best_con, "nlpd", true);
    
    count = count + (mse_uncon > mse_con);

    %% Save results
    results_uncon(r).params_best = params_best_uncon;
    results_uncon(r).log_mse = mse_uncon;
    results_uncon(r).nlpd = nlpd_uncon;
    results_uncon(r).idx_train = idx_train;
    results_uncon(r).idx_val = idx_val;
    results_uncon(r).idx_test = idx_test;

    results_con(r).params_best = params_best_con;
    results_con(r).log_mse = mse_con;
    results_con(r).nlpd = nlpd_con;
    results_con(r).idx_train = idx_train;
    results_con(r).idx_val = idx_val;
    results_con(r).idx_test = idx_test;

    % save('random_split_results_sigma1e7.mat', 'results_con', 'results_uncon');
end



%% Visualization (Paper)

blue   = [0 0.4470 0.7410];
orange = [0.8500 0.3250 0.0980];

sig_exp = 4:8; sigma2 = 10.^(-sig_exp);

mu_con = zeros(size(sig_exp)); sd_con = zeros(size(sig_exp));
mu_uncon = zeros(size(sig_exp)); sd_uncon = zeros(size(sig_exp));

for k = 1:numel(sig_exp)
    i = sig_exp(k);
    fname = sprintf('results/random_split_results_sigma1e%d.mat', i);
    load(fname, 'results_con', 'results_uncon')

    log_mse_con = [results_con.log_mse];
    log_mse_uncon = [results_uncon.log_mse];

    mu_con(k) = mean(log_mse_con); sd_con(k) = std(log_mse_con);
    mu_uncon(k) = mean(log_mse_uncon); sd_uncon(k) = std(log_mse_uncon);
end

figure;
tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold on;

errorbar(sigma2, mu_con, 2*sd_con, '-o', 'Color', blue, 'MarkerFaceColor', 'none', ...
    'LineWidth', 2.5, 'CapSize', 15, 'MarkerSize', 10);
errorbar(sigma2, mu_uncon, 2*sd_uncon, '-o', 'Color', orange, 'MarkerFaceColor', 'none', ...
    'LineWidth', 2.5, 'CapSize', 15, 'MarkerSize', 10);

set(gca, 'XScale', 'log');
xlabel('\sigma^2'); ylabel('log MSPE');
xlim([5e-9 2e-4]); ylim([-17.5 -11.5]);
legend({'cGRF', 'GRF'}, 'Location', 'best');

fontsize(30, "points")


%% Function

function K = cov_uncon(data_j, data_k, params)
    lt = params(1); lx = params(2); ly = params(3); alpha = params(4);
    K = (1/alpha) * QQ(data_j(:,1), data_k(:,1), lt) .* QQ(data_j(:,2), data_k(:,2), lx) .* QQ(data_j(:,3), data_k(:,3), ly);
end

function mu = mean_uncon(data_i, params)
    mu = zeros(size(data_i,1), 1);
end

function K = cov_con(data_j, data_k, params)
    lt = params(1); lx = params(2); ly = params(3); alpha = params(4);
    K = (1/alpha) * QQ(data_j(:,1), data_k(:,1), lt) .* QQ(data_j(:,2), data_k(:,2), lx) .* cQQ(data_j(:,3), data_k(:,3), ly, -37.5, 35, [0, 1, 0, 1]);
    end

function mu = mean_con(data_i, params)
    ly = params(3);
    mu = cm(data_i(:,3), ly, -37.5, 35, [0, 1, 0, 1], zeros(size(data_i(:,1))), (data_i(:,1) - 100) * 0.0050);
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
        % Pointwise
        % se_best = -(log(normpdf(data_test_i(:,4), f_mu, sqrt(f_se.^2 + params_best(end)))));
        % mse_best = mean(se_best);
      
        % Multivariate
        mse_best = -logmvnpdf(data_test_i(:,4), f_mu, f_K + params_best(end) * eye(length(f_mu)));
    end
    
    if training
        result = mse_best;
    else
        result = {mse_best, se_best, f_mu, f_se};
    end
end

function logp = logmvnpdf(x, mu, Sigma)

    x  = x(:);
    mu = mu(:);
    N  = length(x);

    xc = x - mu;
    U = chol(Sigma);   % Sigma = U'*U
    alpha = U \ (U' \ xc);   % alpha = U^-1 * U'^-1 * xc = (U'*U)^-1 * xc = Sigma^-1 * xc
    quad = xc' * alpha;

    logdetSigma = 2 * sum(log(diag(U)));   % detSigma = det(U'*U) = det(U')det(U) = (\Prod_i U_ii)^2

    logp = -0.5 * (N * log(2*pi) + logdetSigma + quad);

end


%%%%% Data Loading %%%%%

function data = load_time_block(data_folder, times, idx_design, use_design, preproc)
    data_cells = cell(numel(times),1);
    for i = 1:numel(times)
        t = times(i);
        fname = fullfile(data_folder, sprintf('Poly-%04d_1.csv', t));
        [X, Y, Z, Vcorr] = read_and_preprocess_file(fname, preproc);

        if use_design
            idx_keep = idx_design(:);
        else
            idx_keep = (1:numel(X)).';
        end

        n = numel(idx_keep);
        data_cells{i} = [repmat(t, n, 1), X(idx_keep), Y(idx_keep), Vcorr(idx_keep), Z(idx_keep)];
    end
    data = vertcat(data_cells{:});
end

function [X, Y, Z, Vcorr, info] = read_and_preprocess_file(fname, preproc)
    T = readtable(fname);
    X = T{:,1}; Y = T{:,2}; Z = T{:,3}; V = T{:,4};
    
    if ~preproc.use_preprocessing
        Vcorr = V;
        info = struct('used', false);
        return;
    end

    % Add preprocessing steps here if needed
end

function idx = equal_spaced_design_indices(X, Y, nx, ny, thresh)
    xmin = min(X); xmax = max(X);
    ymin = min(Y); ymax = max(Y);

    xg = linspace(xmin, xmax, nx);
    yg = linspace(ymin, ymax, ny);
    [Xg, Yg] = meshgrid(xg, yg);
    target = [Xg(:), Yg(:)];
    obs = [X(:), Y(:)];

    % nearest observed points to target grid points
    [idx, dist] = knnsearch(obs, target);
    mask = dist < thresh;     
    idx  = idx(mask);

    idx = unique(idx, 'stable');
end

