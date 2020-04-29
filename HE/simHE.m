%% HE MPC
clc; clear; close all;
yalmip('clear');

% ODE options 'RelTol', 1e-6, 'AbsTol', 1e-6
ode_options = odeset ('RelTol', 1e-6, 'AbsTol', 1e-6, ...
    'NormControl', 'on', 'InitialStep', 1.0e-2, 'MaxStep', 1.0);

%% Load polytope model and observer matrices
% Type of observer(N° observer). Sub-observer(N° sub-observer). Matrix 
RUIO = struct;
UIOO = struct;
load polyModelObs

%% Simulation parameters
Ts = 0.05;                  % Sample time [min]
Time = 15;                 % Simulation end time 
Nsim = Time/Ts;        % Simulation steps
t = 0:Ts:Time-Ts;       % Simulation time
Fail_Q1 = 5; Fail_Q2 = 0.43;    % Actuator fault magnitude [5%, 5%]
Fail_S1 = 2.5; Fail_S2 = -3.5;	% Sensor fault magnitude [0.5% 0.5%]

%% Polytope model and observers
% % This section is commented to reduce simulation time (using pre-calculated observer matrices)
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
% % Observers start point
% Theta_1s = Theta_1s_mid;     % Output fluid 1 temperature (K)
% Theta_2s = Theta_2s_min;     % Output fluid 2 temperature (K)
% run HE_linear;
% x0_obs = [Theta_1s; Theta_2s; Theta_p];
% 
% % System start point
% Theta_1s = Theta_1s_min;     % Output fluid 1 temperature (K)
% Theta_2s = Theta_2s_mid;     % Output fluid 2 temperature (K)
% run HE_linear;
% x0 = [Theta_1s; Theta_2s; Theta_p];
% 
% % Reduced-order unknown input observer
% run HE_DLPV_RUIO;
% 
% % Unknown input output observer
% run HE_DLPV_UIOO;
% 
% % Save models' data
% save polyModelObs.mat

%% Noise
sig = 3e-3*([1 1 1])';         % Ouput noise sigma

rng default;                        % Random seed start
v = sig*randn(1, Nsim);    % Measurement noise v~N(0, sig)

%% Error detection threshold
Tau = 10;                % Convergence period
mag_1 = 1.3e-1;     % Value Q1
mag_2 = 5e-2;        % Value Q2
mag_3 = 6e-2;        % Value O1
mag_4 = 2e-1;     % Value O2

threshold = zeros(4, Nsim);

for k = 1:Nsim
    threshold(1, k) = mag_1 + 1000*exp(-(k-1)/Tau);  % Q1
    threshold(2, k) = mag_2 + 900*exp(-(k-1)/Tau);    % Q2
    threshold(3, k) = mag_3 + 100*exp(-(k-1)/Tau);    % O1
    threshold(4, k) = mag_4 + 1000*exp(-(k-1)/Tau);  % O2
end

%% MPC controller
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

%% Simulation Setup
% Fault Tolerant Control System (1 = disable; 2 = enable).Matrix 
FTCS = struct;

%% Simulation
disp('Simulating...')
FTC_OFF = 1; FTC_ON = 2;
for FT = 1:2    % 1 - FT is off; 2 -  FT is on
    % RUIOs
    for k = 1:N
        RUIO(k).Phi = zeros(N, Nsim+1);     % Observer states
        RUIO(k).X = zeros(nx, Nsim);           % Estimated states
        RUIO(k).error = zeros(1, Nsim);       % Error
        RUIO(k).Fact = zeros(1, Nsim);        % Estimated control Input
        RUIO(k).FQ = zeros(1, Nsim);           % Fault detect Q
        RUIO(k).delay = 0;                            % Detection delay
    end

    % UIOOs
    for k = 1:N
        UIOO(k).Z = zeros(nx, Nsim+1);      % Observer states
        UIOO(k).Ymon = zeros(N, Nsim);     % Monitorated outputs
        UIOO(k).X = zeros(nx, Nsim);           % Estimated states
        UIOO(k).Y = zeros(nx, Nsim);           % Estimated outputs
        UIOO(k).res = zeros(nx, Nsim);        % Residue
        UIOO(k).error = zeros(1, Nsim);       % Error
        UIOO(k).Fsen = zeros(1, Nsim);       % Estimated sensor fault
        UIOO(k).FO = zeros(1, Nsim);          % Fault detect S
    end
    
    % FTCS matrices
    FTCS(FT).U = zeros(nu, Nsim);                   % Control Input
    FTCS(FT).Ufail = zeros(nu, Nsim);              % Faulty control Input
    FTCS(FT).Uff = zeros(nu, Nsim+1);            % Feedforward control Input
    FTCS(FT).Ufails = zeros(nu, Nsim);             % Fails of control inputs
    FTCS(FT).X = zeros(nx, Nsim+1);               % States
    FTCS(FT).Y = zeros(ny, Nsim);                    % Measure outputs
    FTCS(FT).Xsp = zeros(nx, Nsim);                % Set-points  
    FTCS(FT).Y_hat = zeros(ny, Nsim);             % Estimated outputs
    FTCS(FT).Yfail = zeros(ny, Nsim);               % Faulty measure outputs
    FTCS(FT).Obj = zeros(1, Nsim);                  % Objective cost
    FTCS(FT).mu_mhe = zeros(M, Nsim);         % Membership using MHE
    FTCS(FT).mu_fuzzy = zeros(M, Nsim);        % Membership using fuzzy

    FTCS(FT).X_MHE = repmat(x0, 1, N_MHE+1);    % MHE states
    FTCS(FT).U_MHE = repmat(u0, 1, N_MHE);        % MHE inputs
    FTCS(FT).elapsed_time = zeros(Nsim, 1) ;          % MPC elapsed times 

    % Initial states and inputs
    FTCS(FT).X(:, 1) = x0;
    FTCS(FT).Y(:, 1) = C*x0;
    FTCS(FT).Y_hat(:, 1) = C*x0;
    RUIO(FT).X(:, 1) = x0;
    UIOO(FT).X(:, 1) = x0;
    
    % Vector initialization for plots
    FTCS(FT).Yc_sim = C*x0;
    
    if FT == FTC_ON
        disp('Fault tolerant = ON')
    else
        disp('Fault tolerant = OFF')
    end

	for k = 1:Nsim
        tk = k*Ts; % Simulation time
        
        % Set-point changes
        if tk == 4
            Theta_1s = Theta_1s_min;     % Output fluid 1 temperature (K)
            Theta_2s = Theta_2s_mid;     % Output fluid 2 temperature (K)
            run HE_linear;
            xsp = [Theta_1s; Theta_2s; Theta_p];
        elseif tk == 9
            Theta_1s = Theta_1s_max;     % Output fluid 1 temperature (K)
            Theta_2s = Theta_2s_min;     % Output fluid 2 temperature (K)
            run HE_linear;
            xsp = [Theta_1s; Theta_2s; Theta_p];
        end
        FTCS(FT).Xsp(:, k) = xsp;
        
        if k > 1
            FTCS(FT).mu_mhe(:, k) = FTCS(FT).mu_mhe(:, k-1);
            FTCS(FT).Y(:, k) = FTCS(FT).Y(:, k-1);
        end

        FTCS(FT).mu_fuzzy(:, k) =membership(FTCS(FT).Y(:, k), Theta_1s_min, Theta_1s_mid, Theta_1s_max, Theta_2s_min, Theta_2s_mid, Theta_2s_max);
        
        t_tic = tic ;   % to get time evaluated 
        
        [sol, diag] = mhe{FTCS(FT).X_MHE, FTCS(FT).U_MHE, FTCS(FT).mu_mhe(:, k)};
        if diag
            msg = ['Infeasible MHE at t = ', num2str(tk)];
            disp(msg)
            return;
        end
        FTCS(FT).mu_mhe(:, k) = sol;
                
        [sol, diag] = mpc{FTCS(FT).Y(:, k), FTCS(FT).Xsp(:, k), FTCS(FT).mu_mhe(:, k)};
        if diag
            msg = ['Infeasible MPC at t = ', num2str(tk)];
            disp(msg)
            return;
        end
        FTCS(FT).U(:, k) = sol{1}; FTCS(FT).Obj(k) = sol{2};
        
        t_tic = toc(t_tic) ;              % get time elapsed
        FTCS(FT).elapsed_time(k) = t_tic ;   % store the time elapsed for the run
        
        %% Process simulation with ODE
        [tsim, x] = ode45(@(x, u) HE(FTCS(FT).X(:, k), FTCS(FT).U(:, k)) , [0 Ts], FTCS(FT).X(:, k), ode_options);
        FTCS(FT).X(:, k+1) = x(end, :)';
        FTCS(FT).Y(:, k) = C*FTCS(FT).X(:, k);
        
        % MHE horizon update
        FTCS(FT).X_MHE = [FTCS(FT).X_MHE(:, 2:end) FTCS(FT).X(:, k)];
        FTCS(FT).U_MHE = [FTCS(FT).U_MHE(:, 2:end) FTCS(FT).U(:, k)];
    
	end
    
    time_avg = mean(FTCS(FT).elapsed_time) ;
    msg = ['Mean time = ', num2str(time_avg)];
    disp(msg)
    time_avg = max(FTCS(FT).elapsed_time) ;
    msg = ['Max time = ', num2str(time_avg)];
    disp(msg)
    time_avg = min(FTCS(FT).elapsed_time) ;
    msg = ['Min time = ', num2str(time_avg)];
    disp(msg)    

    fig = figure('Name', 'States');
    subplot(311)
    plot(t, FTCS(FT).Xsp(1, :), 'r-.', 'LineWidth', 1.5);
    hold on
    plot(t, FTCS(FT).Y(1, :), 'g--', 'LineWidth', 1.5); hold off
    xlabel('Time [min]'); ylabel('\theta_{1_s} [K]'); grid on
    % axis([0 inf 494 499])
    % leg = legend('Setpoint', 'Estimated', 'Measured', 'Location', 'SouthEast');
    % set(leg, 'Position', [0.748 0.764 0.148 0.109], 'FontSize', 8);
    % leg.ItemTokenSize = [20, 15];
    subplot(312)
    plot(t, FTCS(FT).Xsp(2, :), 'r-.', 'LineWidth', 1.5);
    hold on
    plot(t, FTCS(FT).Y(2, :), 'g--', 'LineWidth', 1.5); hold off
    xlabel('Time [min]'); ylabel('\theta_{2_s} [K]'); grid on
    % axis([0 inf 675 705])
    subplot(313)
    plot(t, FTCS(FT).Y(3, :), 'g--', 'LineWidth', 1.5);
    xlabel('Time [min]'); ylabel('\theta_p [K]'); grid on
    % axis([0 inf 554 562])
        
end

%%
% run enPlotHE