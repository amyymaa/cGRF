clear;
addpath('/Users/yuema/Desktop/cGRF/Code/kernel/')
addpath('/Users/yuema/Desktop/cGRF/Code/PDE-discovery/')

%% Data Generation

x_star = linspace(-10,10,401); % gp_new below takes in -10,10
t_star = linspace(0,1,1001);
m = 0; % symmetry m
sol = pdepe(m, @burgerpde, @burgeric, @burgerbc, x_star, t_star);
surf(sol);

u_star = sol';
x_star = x_star';
t_star = t_star';

%% Parameter Initialization

% active learning
i = 100;            % time to collect data
N0 = 5;             % size of initial batch
n_s = 10;           % size of following batch
ns_max = 72;        % maximum number of interations
criterion = "BIC";  % criterion for sparse regression
eps = 0.1;          % error std of sparse regression, not GPR nugget
thresh = 0.001;     % threshold to break interation

% boundary condition
bc = "dir";
if strcmp(bc, "dir")
    alphaL = 1; betaL = 0; alphaR = 1; betaR = 0; 
elseif strcmp(bc, "neu")
    alphaL = 0; betaL = 1; alphaR = 0; betaR = 1;
elseif strcmp(bc, "rob")
    alphaL = 1; betaL = 1; alphaR = 1; betaR = 1;
end

% save
folder = "results_uncon/";
param = bc+"_11_400_001_"+num2str(eps);

%% Prepare Data (Response)

dt = t_star(i+1) - t_star(i);
du1 = (u_star(:,i+1)-u_star(:,i))/dt;
du1 = du1(2:end-1);
du_true = du1(2:end-1);
x_step = x_star(2)-x_star(1);

%% prepare derivative obervations (Design Matrix)

y_true = u_star(:,i);
dy_true = (y_true(3:end)-y_true(1:end-2))./(x_star(3:end)-x_star(1:end-2));
d2y_true = (y_true(3:end)-2*y_true(2:end-1)+y_true(1:end-2))/(x_step*x_step);
x_star = x_star(2:end-1);
N_star = size(x_star,1);
y_true = y_true(2:end-1);
Theta_true = pool_data(x_star,y_true,dy_true,d2y_true);
[m_x,m_y] = size(Theta_true);


%% true value

rng(1)
% rng(2)

Xi_true = zeros(20,1);
Xi_true(8,1) = -1;
Xi_true(11,1) = 1;
Xi_log = zeros(1,21);
Xi_log(1,8) = 1;
Xi_log(1,11) = 1;

warning('off', 'all')
time = 100;
error = zeros(time,1);
error_l0 = zeros(time,1);
sample = zeros(time,1);

coef = zeros(time,22);

tic
for times=1:time

    chosen_col_last = ones(21,1);
    du = du1 + eps*randn(size(du1)); % randn returns standard normal
    
    %% first bunch of data generated at random
    chosen_index = randsample(N_star, N0);
    x0 = x_star(chosen_index,:);
    u0 = y_true(chosen_index);
    % u0mean = mean(u0);
    % u0 = u0-u0mean;
    s = 0;
    
    while(1)
        
        s = s+1;

        %% GP
        
        if strcmp(folder, 'results_con/')
            gL = alphaL*burgeric(-10)+betaL*Dburgeric(-10); 
            % gL = alphaL*(burgeric(-10)-u0mean)+betaL*Dburgeric(-10); 
            gR = alphaR*burgeric(10)+betaR*Dburgeric(10);      
            % gR = alphaR*(burgeric(10)-u0mean)+betaR*Dburgeric(10); 
            [mu,mu_1,mu_2,sigma] = gp_new(x0,u0,-10,10,[betaL,alphaL,betaR,alphaR,gL,gR]);
        elseif strcmp(folder, 'results_uncon/')
            [mu,mu_1,mu_2,sigma] = gp_new_0(x0,u0);
        end
      

        %% pool data
        % sigma = sigma/(std(u0)^2);
        % y = mu(x_star)+u0mean;
        y = mu(x_star);
        dy = mu_1(x_star);
        d2y = mu_2(x_star);
   
        [Theta] = pool_data(x_star,y,dy,d2y);
        
        %% sparse regression
        Theta_chosen = Theta_true(chosen_index,:); % design matrix
        eta = du(chosen_index,:); % response

        mdl = stepwiselm(Theta_chosen,eta,'Criterion',criterion);
        chosen_col = mdl.Formula.InModel;
        Theta1 = Theta_chosen(:,chosen_col);
        mdl = fitlm(Theta1,eta);
        % tol = mdl.RMSE;
        % tol = tol/std(eta);
        tol = mdl.MSE;
        cof = table2array(mdl.Coefficients(2:end,1));
      
        chosen_col_1 = double(chosen_col)';
        error_1 = calError(chosen_col_1,chosen_col_last);

        if error_1 == 0 % same cols with similar coefs
            error_2 = norm(cof - cof_last,2)/norm(cof_last,2);
            if error_2 < thresh
                break
            end
        end
        
        chosen_col_last = chosen_col_1;
        cof_last = cof;

        
        %% max sample size
        if s > ns_max
            break
        end
        
        
        %% optimal design
        [chosen_index]=optimal_design(Theta_true,Theta,chosen_index,n_s,sigma,tol,x_star);
        x0 = x_star(chosen_index,:);
        u0 = y_true(chosen_index);
        % u0mean = mean(u0);
        % u0 = u0-u0mean;
   
        
        
    end   

    Xi = zeros(m_y,1);
    Xi(chosen_col) = cof;
    error(times,1) = calError(chosen_col,Xi_log);
    error_l0(times,1) = norm(Xi-Xi_true,2);
    sample(times,1) = size(chosen_index,1);
  
    coef(times,chosen_col) = cof;
    coef(times,end) = mdl.RMSE;
    
    disp(times)
    
end
toc


%% Save

% dlmwrite(folder+"result_"+param+".txt", ...
%     [error,error_l0,sample]);
% 
% dlmwrite(folder+"summary_"+param+".txt", ...
%     [mean(error),mean(error_l0),mean(sample); ...
%     std(error),std(error_l0),std(sample); ...
%     mean(error==0),NaN,NaN]);
% 
% dlmwrite(folder+"coef_"+param+".txt",coef);

%% Visualization of percentage of falsely idenfitied items

tiledlayout(1,3, 'TileSpacing', 'compact', 'Padding', 'compact')
epss = [0.1, 0.2, 0.4]; bcs = ["dir", "neu", "rob"];

for bc = bcs
    con = zeros(1,3);
    uncon = zeros(1,3);
    i = 1;
    for eps = epss
        df_con = readmatrix("results_con/summary_"+bc+"_11_400_001_"+num2str(eps)+".txt");
        df_uncon = readmatrix("results_uncon/summary_"+bc+"_11_400_001_"+num2str(eps)+".txt");
        % disp([1-df_con(3,1)', 1-df_uncon(3,1)'])
        con(i) = 1-df_con(3,1);
        uncon(i) = 1-df_uncon(3,1);
        i = i + 1;
    end
    nexttile
    hold on;
    plot([1,2,3], con, 'LineWidth', 1.5, 'Marker','o')
    plot([1,2,3], uncon, 'LineWidth', 1.5, 'Marker','o')
    hold off
    title(bc)
    xlabel("\sigma")
    if bc == "dir"
        ylabel('% False Discovery')
        legend('cGRF', 'GRF', 'Location', 'northwest');
    end
    xticks([1 2 3]); xticklabels([0.1 0.2 0.4]); xlim([0.9 3.1])
    yticks([0 0.3 0.6]); ylim([0, 0.65])
end
fontsize(15, "points")

%% Visualization of number of falsely idenfitied items

tiledlayout(3,3, 'TileSpacing', 'compact', 'Padding', 'compact')
epss = [0.1, 0.2, 0.4]; bcs = ["dir", "neu", "rob"];
for bc = bcs
    for eps = epss
        nexttile
        df_con = readmatrix("results_con/result_"+bc+"_11_400_001_"+num2str(eps)+".txt");
        df_uncon = readmatrix("results_uncon/result_"+bc+"_11_400_001_"+num2str(eps)+".txt");
       
        con_counts = accumarray(df_con(:,1)+1, 1, [4,1]);
        uncon_counts = accumarray(df_uncon(:,1)+1, 1, [4,1]);
        x = 0:3;
        
        bar_width = 0.3;
        offset = bar_width / 2;
        x1 = x - offset; 
        x2 = x + offset; 

        hold on;
        bar(x1, con_counts, bar_width, 'FaceColor', [0 0.4470 0.7410]);
        bar(x2, uncon_counts, bar_width, 'FaceColor', [0.8500 0.3250 0.0980]);
        hold off
        legend('cGRF', 'GRF');
        ylim([0 100]);
        xlim([-0.5 3.5]);

        if bc ~= "rob"
            set(gca, 'XTick', [])
        else
            xlabel('Number of False Discoveries')
        end
        if bc == "dir"
            title("\sigma="+num2str(eps))
        end

        if eps == 0.1
            ylabel("Frequency");
        else
            set(gca, 'YTick', [])
        end
        if eps == 0.4
            yyaxis right
            set(gca, 'YTick', [], 'YColor', 'k')
            h = ylabel(bc);
            set(h,'Rotation',-90,'VerticalAlignment','bottom','HorizontalAlignment','center');
        end
       

        % ylim([0 4]) % 2 outliers for GRF and 1 for cGRF are excluded
    end
end    
fontsize(15, "points")

%% Visualization of MSE of coefficient estimation

tiledlayout(3,3, 'TileSpacing', 'compact', 'Padding', 'compact')
epss = [0.1, 0.2, 0.4]; bcs = ["dir", "neu", "rob"];
for bc = bcs
    for eps = epss
        nexttile
        df_con = readmatrix("results_con/result_"+bc+"_11_400_001_"+num2str(eps)+".txt");
        df_uncon = readmatrix("results_uncon/result_"+bc+"_11_400_001_"+num2str(eps)+".txt");
        bx = boxplot([log(df_con(:,2)), log(df_uncon(:,2))], 'OutlierSize', 2);
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

        if bc ~= "rob" 
            set(gca, 'XTick', [])
        else
            xticklabels(["cGRF", "GRF"])
        end
        if bc == "dir"
            title("\sigma="+num2str(eps))
        end

        if eps == 0.1
            ylabel("log(MSE(\beta))");
        else
            set(gca, 'YTick', [])
        end

        if eps == 0.4
            yyaxis right
            set(gca, 'YTick', [], 'YColor', 'k')
            h = ylabel(bc);
            set(h,'Rotation',-90,'VerticalAlignment','bottom','HorizontalAlignment','center');
        end
        ylim([-7 2]) % 2 outliers for GRF and 1 for cGRF are excluded
    end
end    
fontsize(15, "points")

%% Functions

function [c,f,s] = burgerpde(x,t,u,dudx)
    c = 1;
    f = (1) * dudx;
    s = (-1) * u * dudx;
end


function u0 = burgeric(x)
    u0 = 2*exp(-15*(x-9.5)^2) + 1.5*exp(-15*(x+1)^2) + 1*exp(-25*(x+9.5)^2);
end

function du0 = Dburgeric(x)
    du0 = -60 * (x - 9.5) .* exp(-15 * (x - 9.5).^2) ...
        - 45 * (x + 1) .* exp(-15 * (x + 1).^2) ...
        - 50 * (x + 9.5) .* exp(-25 * (x + 9.5).^2);
end

function [pl, ql, pr, qr] = burgerbc(xl, ul, xr, ur, t)
    % Left boundary conditions (Robin form)
    alphaL = 1;   % Coefficient for u (Dirichlet: 1, Neumann: 0)
    betaL = 0;    % Coefficient for du/dx (Dirichlet: 0, Neumann: 1)
    gL = alphaL*burgeric(xl)+betaL*Dburgeric(xl);       % Boundary value (u or du/dx depending on type)

    pl = alphaL * ul - gL;  % p at left boundary
    ql = betaL;             % q at left boundary

    % Right boundary conditions (Robin form)
    alphaR = 1;   % Coefficient for u (Dirichlet: 1, Neumann: 0)
    betaR = 0;    % Coefficient for du/dx (Dirichlet: 0, Neumann: 1)
    gR = alphaR*burgeric(xr)+betaR*Dburgeric(xr);       % Boundary value (u or du/dx depending on type)

    pr = alphaR * ur - gR;  % p at right boundary
    qr = betaR;             % q at right boundary
end