%% 
clear; clc;
addpath("/Users/yuema/Desktop/cGRF/Code/kernel/")
addpath("/Users/yuema/Desktop/cGRF/Code/PDE-solver/")
warning('off','all')

%% Initializations

rep = 100;

% Spatial discretization
x_start               = 0;
x_end                 = 1;
NumOfSpatialPoints    = 16;
NumOfSpatialPoints_   = NumOfSpatialPoints*3;
int_step_x            = ((x_end-0.01)-(x_start+0.01))/(NumOfSpatialPoints-1); % start from interior
Spatial               = (x_start+0.01):int_step_x:(x_end-0.01); 
int_step_x_           = (x_end-x_start)/(NumOfSpatialPoints_-1); % finer grid for prediction
Spatial_              = x_start:int_step_x_:x_end; 

% Time discretization
t_start               = 0;
t_end                 = 0.25;
NumOfTimePoints       = 100;
int_step              = t_end/NumOfTimePoints;
Time                  = t_start:int_step:((NumOfTimePoints-1)*int_step);

% Solution
StateData             = zeros(NumOfSpatialPoints*rep, NumOfTimePoints);
FirstDerivData        = zeros(NumOfSpatialPoints*rep, NumOfTimePoints);
Deriv_Data            = zeros(NumOfSpatialPoints, NumOfTimePoints);
Deriv_Data(:, 1)      = -cos(Spatial*pi/x_end)*(pi/x_end)^2;

%% Covariance Function

% Time covariance functions
HP_Time    = [7*int_step, 1];
CC_Time    = (1/HP_Time(2)) * cRR(Time, Time, HP_Time(1), t_start, t_end, [0, 1, 0, 0]);
CI_Time    = (1/HP_Time(2)) * cQR(Time, Time, HP_Time(1), t_start, t_end, [0, 1, 0, 0])';
II_Time    = (1/HP_Time(2)) * cQQ(Time, Time, HP_Time(1), t_start, t_end, [0, 1, 0, 0]);

% Spatial covariance functions
HP_Space   = [7*int_step_x, 1];
coeffs     = [1, 1, 1, 0];
CC_Space   = (1/HP_Space(2)) * cQQ(Spatial, Spatial, HP_Space(1), x_start, x_end, coeffs);
CD2_Space  = (1/HP_Space(2)) * cQS(Spatial, Spatial, HP_Space(1), x_start, x_end, coeffs);
D2D2_Space = (1/HP_Space(2)) * cSS(Spatial, Spatial, HP_Space(1), x_start, x_end, coeffs);

% Spatial indices
Idx        = 1:NumOfSpatialPoints;

%% Multiple solver simulations

for j = 1:rep
disp(j)
%%% Sequential updating

% Initial noise covariance
NoiseCov = zeros(NumOfSpatialPoints, 1)*1e-8;

for i = 2:NumOfTimePoints
    
    if i == 2
        % calculate exactly
        Deriv_Data(:,2) = -cos(Spatial*pi/x_end)*(pi/x_end)^2;

    else
        % Estimate using GP predictive distribution
        A = kron(CC_Time(1:i-1, 1:i-1), CC_Space(Idx, Idx));
        C = kron(CI_Time(1:i-1, i), CD2_Space(Idx, Idx));
        B = D2D2_Space(Idx, Idx).*II_Time(i,i);

        Obs_Values = Deriv_Data(Idx,1:i-1); % get all
        Obs_Values = Obs_Values(:)'; % concatenate
        
        GP_D2_Mean = C' * ( (A+diag(sparse(NoiseCov'))) \ Obs_Values' );
        GP_D2_Var  = B - C' / (A+diag(sparse(NoiseCov'))) * C;
        
        Deriv_Data(:,i) = Deriv_Data(:, 1) + GP_D2_Mean + ...
            ( randn(1, NumOfSpatialPoints)*chol(GP_D2_Var+eye(NumOfSpatialPoints)*1e-8) )';

    end
    
    % Predict variance of derivative observations at next time point - mean and variance  
    A = kron(CC_Time(1:i-1, 1:i-1), CC_Space(Idx, Idx));
    C = kron(CC_Time(1:i-1, i), CC_Space(Idx, Idx));
    B = CC_Space(Idx, Idx).*CC_Time(i,i);
    
    GP_C_Var  = B - C' / (A+diag(sparse(NoiseCov'))) * C;
    
    NewNoiseCov = diag(GP_C_Var); 
    NoiseCov = [NoiseCov; NewNoiseCov];
        
end

%%% Solution

CC_Space_A = (1/HP_Space(2)) * cQQ(Spatial, Spatial, HP_Space(1), x_start, x_end, coeffs);
CC_Space_B = (1/HP_Space(2)) * [cQQ(Spatial_, Spatial_, HP_Space(1), x_start, x_end, coeffs), ...
                                cQR(Spatial_, Spatial_, HP_Space(1), x_start, x_end, coeffs); ...
                                cQR(Spatial_, Spatial_, HP_Space(1), x_start, x_end, coeffs)',...
                                cRR(Spatial_, Spatial_, HP_Space(1), x_start, x_end, coeffs)];
CC_Space_C = (1/HP_Space(2)) * [cQQ(Spatial, Spatial_, HP_Space(1), x_start, x_end, coeffs), ...
                                cQR(Spatial, Spatial_, HP_Space(1), x_start, x_end, coeffs)];

A = kron(CC_Time, CC_Space_A);
C = kron(CI_Time, CC_Space_C);
B = kron(II_Time, CC_Space_B);

Obs_Values = Deriv_Data(:)';

GP_I_Mean = C' * ( (A+diag(sparse(NoiseCov'))) \ Obs_Values' );
GP_I_Var  = B - C' / (A+diag(sparse(NoiseCov'))) * C;

GP_I_Mean = reshape(GP_I_Mean, NumOfSpatialPoints_*2, NumOfTimePoints);
GP_I_Mean_state = GP_I_Mean(1:NumOfSpatialPoints_, :) + ...
    repmat((2+cos(Spatial_*pi/x_end))', 1, NumOfTimePoints);
GP_I_Mean_deriv = GP_I_Mean((NumOfSpatialPoints_+1):end, :) + ...
    repmat((-sin(Spatial_*pi/x_end)*pi/x_end)', 1, NumOfTimePoints);

surf(GP_I_Mean_state);
surf(GP_I_Mean_deriv);

StateData( ((j-1)*NumOfSpatialPoints_+1):j*NumOfSpatialPoints_, : ) = GP_I_Mean_state;
FirstDerivData( ((j-1)*NumOfSpatialPoints_+1):j*NumOfSpatialPoints_, : ) = GP_I_Mean_deriv;

end

%%
GP_I_Mean = reshape(GP_I_Mean, NumOfSpatialPoints_*2, NumOfTimePoints);
GP_I_Mean_state = GP_I_Mean(1:NumOfSpatialPoints_, :) + ...
    repmat((2+cos(Spatial_*pi/x_end))', 1, NumOfTimePoints);
GP_I_Mean_deriv = GP_I_Mean((NumOfSpatialPoints_+1):end, :) + ...
    repmat((-sin(Spatial_*pi/x_end)*pi/x_end)', 1, NumOfTimePoints);

surf(GP_I_Mean_state);
surf(GP_I_Mean_deriv);
surf(GP_I_Mean_state + GP_I_Mean_deriv);

%%

save("heat_Rob_16100_state.mat", "StateData")
save("heat_Rob_16100_deriv.mat", "FirstDerivData")

%% Numeric solver (MATLAB)

m = 0; % symmetry m
sol = pdepe(m, @heatpde, @heatic, @heatbc, Spatial_, Time);
% sol = sol';
sol_deriv = sol;
for k = 1:NumOfTimePoints
    sol_deriv(k,:) = gradient(squeeze(sol(k,:)), Spatial_);
end

sol = sol';
sol_deriv = sol_deriv';
% Manual correction needed for numerical solver
sol_deriv(NumOfSpatialPoints_, :) = 0;
sol(1, :) = 3 - sol_deriv(1, :);

%%

StateData_16100 = load("results/heat_Rob_16100_state.mat");
FirstDerivData_16100 = load("results/heat_Rob_16100_deriv.mat");
rep = 100;

%%

State = zeros(NumOfSpatialPoints_, NumOfTimePoints);
for i=0:(rep-1)
    State = State + StateData_16100.StateData(i*48+1:(i+1)*48,:);
end
State = State/rep;


%% Visualization

blue   = [0 0.4470 0.7410];
orange = [0.8500 0.3250 0.0980];
xrange1 = [0, 0.1];
xrange2 = [0.9 1];
yrange1 = [1.9 3.1];
yrange2 = [-0.9 0.1];
quantiles = [0 1];

tiledlayout(3,3, 'TileSpacing', 'compact', 'Padding', 'compact')

i = 9; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
StateData_i = reshape(StateData_16100.StateData(:, i), NumOfSpatialPoints_, rep);
FirstDerivData_i = reshape(FirstDerivData_16100.FirstDerivData(:, i), NumOfSpatialPoints_, rep);

nexttile
lb = quantile(StateData_i + FirstDerivData_i, quantiles(1), 2);
ub = quantile(StateData_i + FirstDerivData_i, quantiles(2), 2);
avg = quantile(StateData_i + FirstDerivData_i, 0.5, 2);
func_num_matlab = sol(:, i) + sol_deriv(:, i);
plot_helper(Spatial_, lb, ub, avg, func_num_matlab, xrange1, yrange1)
ylabel("u+u_x")
title("Robin")
ax = gca;
ax.XTick = [];

nexttile([3 1])
imagesc(State)
hold on
xline(9, "Color", "k", "LineWidth", 2, "LineStyle", "--")
xline(49, "Color", "k", "LineWidth", 2, "LineStyle", "--")
xline(89, "Color", "k", "LineWidth", 2, "LineStyle", "--")
xline(0.5, "Color", "r", "LineWidth", 3)
yline(0.5, "Color", "r", "LineWidth", 3)
yline(48.4, "Color", "r", "LineWidth", 3)
hold off
ylabel("x") 
title("u")
ax = gca;
ax.YDir = 'normal';
ax.YAxisLocation = 'right';
yticks([1 4*3 7*3 10*3 13*3 16*3])
yticklabels({'0','0.2','0.4','0.6','0.8','1'})
ax.XTick = [];

nexttile
lb = quantile(FirstDerivData_i, quantiles(1), 2);
ub = quantile(FirstDerivData_i, quantiles(2), 2);
avg = quantile(FirstDerivData_i, 0.5, 2);
func_num_matlab = sol_deriv(:, i);
plot_helper(Spatial_, lb, ub, avg, func_num_matlab, xrange2, yrange2)
h = ylabel("u_x");
title("Neumann")
ax = gca;
ax.YAxisLocation = 'right';
set(h,'Rotation',-90,'VerticalAlignment','bottom','HorizontalAlignment','center');
ax = gca;
ax.XTick = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i = 49; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
StateData_i = reshape(StateData_16100.StateData(:, i), NumOfSpatialPoints_, rep);
FirstDerivData_i = reshape(FirstDerivData_16100.FirstDerivData(:, i), NumOfSpatialPoints_, rep);

nexttile
lb = quantile(StateData_i + FirstDerivData_i, quantiles(1), 2);
ub = quantile(StateData_i + FirstDerivData_i, quantiles(2), 2);
avg = quantile(StateData_i + FirstDerivData_i, 0.5, 2);
func_num_matlab = sol(:, i) + sol_deriv(:, i);
plot_helper(Spatial_, lb, ub, avg, func_num_matlab, xrange1, yrange1)
ylabel("u+u_x")
ax = gca;
ax.XTick = [];
nexttile
lb = quantile(FirstDerivData_i, quantiles(1), 2);
ub = quantile(FirstDerivData_i, quantiles(2), 2);
avg = quantile(FirstDerivData_i, 0.5, 2);
func_num_matlab = sol_deriv(:, i);
plot_helper(Spatial_, lb, ub, avg, func_num_matlab, xrange2, yrange2)
h = ylabel("u_x");
ax = gca;
ax.YAxisLocation = 'right';
set(h,'Rotation',-90,'VerticalAlignment','bottom','HorizontalAlignment','center');
ax = gca;
ax.XTick = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i = 89; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
StateData_i = reshape(StateData_16100.StateData(:, i), NumOfSpatialPoints_, rep);
FirstDerivData_i = reshape(FirstDerivData_16100.FirstDerivData(:, i), NumOfSpatialPoints_, rep);

nexttile
lb = quantile(StateData_i + FirstDerivData_i, quantiles(1), 2);
ub = quantile(StateData_i + FirstDerivData_i, quantiles(2), 2);
avg = quantile(StateData_i + FirstDerivData_i, 0.5, 2);
func_num_matlab = sol(:, i) + sol_deriv(:, i);
plot_helper(Spatial_, lb, ub, avg, func_num_matlab, xrange1, yrange1)
xlabel("x")
ylabel("u+u_x")
nexttile
lb = quantile(FirstDerivData_i, quantiles(1), 2);
ub = quantile(FirstDerivData_i, quantiles(2), 2);
avg = quantile(FirstDerivData_i, 0.5, 2);
func_num_matlab = sol_deriv(:, i);
plot_helper(Spatial_, lb, ub, avg, func_num_matlab, xrange2, yrange2)
xlabel("x")
h = ylabel("u_x");
ax = gca;
ax.YAxisLocation = 'right';
set(h,'Rotation',-90,'VerticalAlignment','bottom','HorizontalAlignment','center');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fontsize(20, "points")

%% Functions

function [c,f,s] = heatpde(x,t,u,dudx)
c = 1;
f = dudx;
s = 0;
end

function u0 = heatic(x)
u0 = 2+cos(x*pi);
end

function [pl,ql,pr,qr] = heatbc(xl,ul,xr,ur,t)
pl = ul-3;
ql = 1;
pr = 0;
qr = 1;
end

function plot_helper(x, lb, ub, avg, utrue, xlims, ylims)
% blue   = [0 0.4470 0.7410];
orange = [0.8500 0.3250 0.0980];

fill([x fliplr(x)], [lb' fliplr(ub')], "b", 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on
plot(x, avg, "LineWidth", 2, "Color", orange, "LineStyle", "-");
plot(x, utrue, "LineWidth", 2, "Color", "k", "LineStyle", "--");
hold off
if ~isnan(xlims(1))
    xlim(xlims)
end
if ~isnan(ylims(1))
    ylim(ylims)
end
end