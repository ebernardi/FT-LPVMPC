%% HE MPC
clc; clear; close all;
yalmip('clear');

% ODE options 'RelTol', 1e-6, 'AbsTol', 1e-6
ode_options = odeset ('RelTol', 1e-6, 'AbsTol', 1e-6, ...
    'NormControl', 'on', 'InitialStep', 1.0e-2, 'MaxStep', 1.0);

%% Load polytope model
load polyModel

%% Simulation parameters
Ts = 0.05;                  % Sample time [min]
Time = 15;                 % Simulation end time 
Nsim = Time/Ts;        % Simulation steps
t = 0:Ts:Time-Ts;       % Simulation time

% %% LPV Models
% Theta_1s_min = 495;              % Minimum output fluid 1 temperature (K)
% Theta_1s_mid = 497.32;         % Middle output fluid 1 temperature (K)
% Theta_1s_max = 500;             % Maximum output fluid 1 temperature (K)
% 
% Theta_2s_min = 680;              % Minimum output fluid 2 temperature (K)
% Theta_2s_mid = 695.915;       % Middle output fluid 2 temperature (K)
% Theta_2s_max = 710;             % Maximum output fluid 2 temperature (K)
% 
% N = 2;                                      % Number of parameters
% L = 3;                                       % Linearization points per parameter
% M = L^N;                                 % Number of models
% run HE_polytope;                     % M models
% 
% % Save models' data
% save polyModel.mat

%% MPC controller
% Constraints
% Constraints
xmin = [495; 650; 530];
xmax = [500; 750; 590];
umin = [90; 7];
umax = [150; 20];

% Wheight matrix
Qx = sys(1).Cd'*sys(1).Cd;
lambda = 0.1;
Ru = lambda*diag([1 10]);

Theta_1s = Theta_1s_mid;     % Output fluid 1 temperature (K)
Theta_2s = Theta_2s_mid;     % Output fluid 2 temperature (K)
run HE_linear;
x0 = [Theta_1s; Theta_2s; Theta_p];
u0 = [Q1; Q2];

Theta_1s = Theta_1s_min;     % Output fluid 1 temperature (K)
Theta_2s = Theta_2s_max;     % Output fluid 2 temperature (K)
run HE_linear;
xsp = [Theta_1s; Theta_2s; Theta_p];

%% Constraint sets
Z = Polyhedron('lb', [xmin; umin], 'ub', [xmax; umax]); % Extended set

X = projection(Z, 1:nx); X = minHRep(X);
U = projection(Z, nx+1:nx+nu); U = minHRep(U);

for i = 1:M
    [sys(i).Klqr, sys(i).Plqr] = dlqr(sys(i).Ad, sys(i).Bd, Qx, Ru);
end

%% MHE
N_MHE = 5;
run MHE

%% MPC
N_MPC = 5;
run MPC

%% Simulation
% Vector initialization for plots
state = x0; input = [0; 0]; Y = C*x0; X = x0; tsim = 0;
Xsp = xsp;

X_MHE = repmat(x0, 1, N_MHE+1); U_MHE = repmat(u0, 1, N_MHE); 

mu_mhe = zeros(M, Nsim+1);
mu_fuzzy = zeros(M, Nsim);
obj = zeros(1, Nsim);                  % Objective cost

disp('Iniciando...')
elapsed_time = zeros(Nsim, 1) ;  % initialize the elapsed times 
for j = 1:Nsim
    tk = j*Ts;
    
	if tk == 4
        Theta_1s = Theta_1s_mid;     % Output fluid 1 temperature (K)
        Theta_2s = Theta_2s_mid;     % Output fluid 2 temperature (K)
        run HE_linear;
        xsp = [Theta_1s; Theta_2s; Theta_p];
    elseif tk == 9
        Theta_1s = Theta_1s_mid;     % Output fluid 1 temperature (K)
        Theta_2s = Theta_2s_min;     % Output fluid 2 temperature (K)
        run HE_linear;
        xsp = [Theta_1s; Theta_2s; Theta_p];
	end

    t_tic = tic ;   % to get time evaluated 
        
    mu_fuzzy(:, j) =membership(X(:, j), Theta_1s_min, Theta_1s_mid, Theta_1s_max, Theta_2s_min, Theta_2s_mid, Theta_2s_max);

    [sol, diag] = mhe{X_MHE, U_MHE, mu_mhe(:, j)};
    if diag
        msg = ['Infeasible MHE at t = ', num2str(tk)];
        disp(msg)
        return;
    end
    mu_mhe(:, j+1) = sol;

    [sol, diag] = mpc{X(:, j), xsp, mu_mhe(:, j+1)};
    if diag
        msg = ['Infeasible MPC at t = ', num2str(tk)];
        disp(msg)
        return;
    end
    umpc(:, j) = sol{1}; 
    obj(j) = sol{2};

    t_tic = toc(t_tic) ;              % get time elapsed
    elapsed_time(j) = t_tic ;   % store the time elapsed for the run

    % Continuous-time simulation (reality)
	[tc, x] = ode45(@(x, u) HE(X(:, j), umpc(:, j)), [0 Ts], X(:, j), ode_options);
    X(:, j+1) = x(end, :)';           % Discrete state vector

    Y = [Y C*x(2:end, :)'];       % Output vector
    Xsp = [Xsp xsp];               % State vector
    tsim = [tsim tk];               % Time vector
    
    X_MHE = [X_MHE(:, 2:end) X(:, j)];
    U_MHE = [U_MHE(:, 2:end) umpc(:, j)];
end
mu_mhe = mu_mhe(:, 2:end);

%%
run enPlotHE