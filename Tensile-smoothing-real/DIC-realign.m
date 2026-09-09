%% Settings

% old folder with original data
src_dir = "/Users/yuema/Desktop/cGRF/Code/Tensile-smoothing-real/data"; 
% new folder with realigned data
out_dir = "/Users/yuema/Desktop/cGRF/Code/Tensile-smoothing-real/data1";     

times = 100:1:200;

fname_fun = @(t) fullfile(src_dir, sprintf('Poly-%04d_0.csv', t));

base_time = 100;

%% 1. Use time 200 as dense reference design

obs_fname = fname_fun(200);
Dref = readmatrix(obs_fname);

Xref = Dref(:,1);
Yref = Dref(:,2);
Zref = Dref(:,3);

obs_ref = [Xref, Yref, Zref];

%% 2. Load baseline displacement at time 100

D100 = readmatrix(fname_fun(base_time));

X100 = D100(:,1);
Y100 = D100(:,2);
Z100 = D100(:,3);
V100 = D100(:,5);

obs100 = [X100, Y100, Z100];

[idx100, dist100] = knnsearch(obs100, obs_ref);

tol = 1e-8;
if max(dist100) > tol
    warning("Some reference points do not exactly match time 100. Max distance = %.3g", max(dist100));
end

V_base = V100(idx100);

%% 3. For each time, compute relative displacement wrt time 100

for t = times
    Dt = readmatrix(fname_fun(t));

    Xt = Dt(:,1);
    Yt = Dt(:,2);
    Zt = Dt(:,3);
    Vt = Dt(:,5);

    obst = [Xt, Yt, Zt];

    % Match current time to reference design
    [idxt, distt] = knnsearch(obst, obs_ref);

    if max(distt) > tol
        warning("Time %d: some points do not exactly match ref design. Max distance = %.3g", t, max(distt));
    end

    V_current = Vt(idxt);

    V_rel = V_current - V_base;

    % Save new dense data
    out = [Xref, Yref, Zref, V_rel];

    out_fname = fullfile(out_dir, sprintf('Poly-%04d_1.csv', t));
    writematrix(out, out_fname);
end