%% CSTR LPV_FT-MPC
clc; clear; close all;
yalmip('clear');

% ODE options 'RelTol', 1e-6, 'AbsTol', 1e-6
ode_options = odeset ('RelTol', 1e-6, 'AbsTol', 1e-6, ...
    'NormControl', 'on', 'InitialStep', 1.0e-2, 'MaxStep', 1.0);

%% Load polytope model and observer matrices
% Type of observer(N° observer). Sub-observer(N° sub-observer). Matrix 
RUIO = struct;
UIOO = struct;

% Fault Tolerant Control System (1 = disable; 2 = enable).Matrix 
FTCS = struct;
sys = struct;

% Load models saved
load polyModelObs

%% Simulation parameters
Ts = 0.05;                  % Sample time [min]
Time = 80;                 % Simulation end time 
Nsim = Time/Ts;        % Simulation steps
t = 0:Ts:Time-Ts;       % Simulation time

Fail_Q1 = -5; Fail_Q2 = 5;      % Actuator fault magnitude [5%, 5%]
Fail_S1 = 2; Fail_S3 = -1.5;    % Sensor fault magnitude [2% 0.5%]

%% Polytope model and observers
% This section is commented to reduce simulation time (using pre-calculated observer matrices)
% V_min = 90;             % Minimum volume (m^3)
% V_mid = 98;             % Middle volume (m^3)
% V_max = 110;          % Maximum volume (m^3)
% 
% Tr_min = 440;           % Mimimum temperature (°K)
% Tr_mid = 445;           % Middle temperature (°K)
% Tr_max = 450;          % Maximum temperature (°K)
% 
% N = 2;                       % Number of parameters
% L = 3;                        % Linearization points per parameter
% M = L^N;                  % Number of models
% run CSTR_polytope;  % M Models
% 
% % Observers start point
% Vr = V_mid;                  % [l] Reactor volume
% Tr = Tr_min;                  % [K] Output temperature
% run CSTR_linear;
% x0_obs = [Vr; Ca; Tr];
% 
% % Reduced-order unknown input observer
% run CSTR_DLPV_RUIO;
% 
% % Unknown input output observer
% run CSTR_DLPV_UIOO;
% 
% % Save models' data
% save polyModelObs.mat

%% Noise
sig = 5e-15*([1 5 2])';       % Ouput noise sigma

rng default;                        % Random seed start
v = sig*randn(1, Nsim);     % Measurement noise v~N(0, sig)

%% Error detection threshold
Tau = 2;               % Convergence period
mag_1 = 9e-2;     % Value Q1
mag_2 = 5e-4;     % Value Q2
mag_3 = 2e-7;     % Value O1
mag_4 = 5e-3;     % Value O2

threshold = zeros(4, Nsim);

for k = 1:Nsim
    threshold(1, k) = mag_1 + 450*exp(-(k-1)/Tau);  % Q1
    threshold(2, k) = mag_2 + 100*exp(-(k-1)/Tau);  % Q2
    threshold(3, k) = mag_3 + 0.6*exp(-(k-1)/Tau);   % O1
    threshold(4, k) = mag_4 + 400*exp(-(k-1)/Tau);  % O2
end

%% MPC controller
% Constraints
xmin = [90; 0.03; 430];
xmax = [110; 0.17; 460];
umin = [90; 85];
umax = [110; 105];

% Wheight matrix
Qx = Cd'*Cd;
lambda = 0.1;
Ru = lambda*diag([2 1]);

%% Constraint sets
Z = Polyhedron('lb', [xmin; umin], 'ub', [xmax; umax]); % Extended set

X = projection(Z, 1:nx); X = minHRep(X);
U = projection(Z, nx+1:nx+nu); U = minHRep(U);

for i = 1:M
    [sys(i).Klqr, sys(i).Plqr] = dlqr(sys(i).Ad, sys(i).Bd, Qx, Ru);
end

%% MHE
N_MHE = 8;
run MHE

%% MPC
N_MPC = 8;
run MPC

%% Simulation
disp('Simulating...')
FTC_OFF = 1; FTC_ON = 2;
for FT = FTC_OFF:FTC_ON    % 1 - FT is off; 2 -  FT is on
    
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
    
    % Initial state
    Vr = V_mid;                  % [l] Reactor volume
    Tr = Tr_min;                  % [K] Output temperature
    run CSTR_linear;
    x0 = [Vr; Ca; Tr];
    u0 = [Qs; Qc];
    
    % Initial set-point
    Vr = V_min;                  % [l] Reactor volume
    Tr = Tr_max;                  % [K] Output temperature
	Ca = CAe/(1+(k0*(Vr/qe)*exp(-E_R/Tr)));
    xsp = [Vr; Ca; Tr];
    
    % FTCS matrices
    FTCS(FT).U = zeros(nu, Nsim);                   % Control Input
    FTCS(FT).Ufail = zeros(nu, Nsim);              % Faulty control Input
    FTCS(FT).Ufails = zeros(nu, Nsim);             % Fails of control inputs
    FTCS(FT).Uff = zeros(nu, Nsim+1);            % Feedforward control Input
    FTCS(FT).X = zeros(nx, Nsim+1);               % States
    FTCS(FT).Y = zeros(ny, Nsim);                    % Measure outputs
    FTCS(FT).Xsp = zeros(nx, Nsim);                % Set-points  
    FTCS(FT).Y_hat = zeros(ny, Nsim);             % Estimated outputs
    FTCS(FT).Yfail = zeros(ny, Nsim);               % Faulty measure outputs
    FTCS(FT).Yfails = zeros(ny, Nsim);              % Fails of measure outputs
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
    RUIO(FT).X(:, 1) = x0_obs;
    UIOO(FT).X(:, 1) = x0_obs;
    
    % Display status of the FTCS
    if FT == FTC_ON
        disp('Fault tolerant = ON')
    else
        disp('Fault tolerant = OFF')
    end

    % Simulation loop
	for k = 1:Nsim

        tk = k*Ts;      % Simulation time
        
%     %% Setpoint
%     Xsp(1, k) = V_mid-2;
%     Xsp(3, k) = Tr_mid;
% 
%     if tk < 160
%         Xsp(3, k) = Tr_min;
%     elseif tk >= 160 && tk < 200
%         Xsp(3, k) = Tr_min+((Tr_mid-Tr_min)*(tk-160)/40);
%     elseif tk >= 200 && tk < 340
%         Xsp(3, k) = Tr_mid;   
%     elseif tk >= 340 && tk < 380
%         Xsp(1, k) = V_mid-2+((V_max-V_mid+2)*(tk-340)/40);
%     elseif tk >= 380 && tk < 520
%         Xsp(1, k) = V_max;
%     elseif tk >= 520 && tk < 560
%         Xsp(3, k) = Tr_mid+((Tr_max-Tr_mid)*(tk-520)/40);
%         Xsp(1, k) = V_max;
%     elseif tk >= 560 && tk < 740
%         Xsp(3, k) = Tr_max;
%         Xsp(1, k) = V_max;
%     end
        
        % Set-point changes
        if tk == 4
            Vr = V_mid-2;
            Tr = Tr_min;
            Ca = CAe/(1+(k0*(Vr/qe)*exp(-E_R/Tr)));
            xsp = [Vr; Ca; Tr];
        elseif tk >= 9 && tk < 19
            Vr = V_mid-2;
            Tr = Tr_min+((Tr_mid-Tr_min)*(tk-9)/10);
            Ca = CAe/(1+(k0*(Vr/qe)*exp(-E_R/Tr)));
            xsp = [Vr; Ca; Tr];
        end
        FTCS(FT).Xsp(:, k) = xsp;

        % Update data
        if k > 1
            FTCS(FT).mu_mhe(:, k) = FTCS(FT).mu_mhe(:, k-1);
            FTCS(FT).Y_hat(:, k) = FTCS(FT).Y_hat(:, k-1);
            FTCS(FT).Uff(:, k) = FTCS(FT).Uff(:, k-1);
        end
        
        t_tic = tic;    % To get time evaluated 
        
        %% MHE QP
        [sol, diag] = mhe{FTCS(FT).X_MHE, FTCS(FT).U_MHE, FTCS(FT).mu_mhe(:, k)};
        if diag
            msg = ['Infeasible MHE at t = ', num2str(tk)];
            disp(msg)
            return;
        end
        FTCS(FT).mu_mhe(:, k) = sol;
        
        %% MPC QP
        [sol, diag] = mpc{FTCS(FT).Y_hat(:, k), FTCS(FT).Xsp(:, k), FTCS(FT).Uff(:, k), FTCS(FT).mu_mhe(:, k)};
        if diag
            msg = ['Infeasible MPC at t = ', num2str(tk)];
            disp(msg)
            return;
        end
        FTCS(FT).U(:, k) = sol{1}; FTCS(FT).Obj(k) = sol{2};
        
        t_tic = toc(t_tic) ;              % Get time elapsed
        FTCS(FT).elapsed_time(k) = t_tic ;   % Store the time elapsed for the run

        %% Actuator fault income
        % No fault
        FTCS(FT).Ufails(:, k) = [0; 0];
        FTCS(FT).Ufail(:, k) = FTCS(FT).U(:, k);
 
        % Q1 fault
        if tk > 35 && tk < 45
            FTCS(FT).Ufails(:, k) = [-Fail_Q1; 0];
            FTCS(FT).Ufail(:, k) = FTCS(FT).U(:, k) + FTCS(FT).Ufails(:, k);
        end
 
        % Q2 fault
        if tk > 20 && tk < 30
            FTCS(FT).Ufails(:, k) = [0; -Fail_Q2+Fail_Q2*(exp(-2*(tk-20)/1))];
            FTCS(FT).Ufail(:, k) = FTCS(FT).U(:, k) + FTCS(FT).Ufails(:, k);
        elseif tk >= 30 && tk < 32
            FTCS(FT).Ufails(:, k) = [0; -Fail_Q2*(exp(-8*(tk-30)/1))];
            FTCS(FT).Ufail(:, k) = FTCS(FT).U(:, k) + FTCS(FT).Ufails(:, k);
        end

        % Natural system saturation
        for j = 1:nu
            if FTCS(FT).Ufail(j, k) >= umax(j)
                FTCS(FT).Ufail(j, k) = umax(j);
            elseif FTCS(FT).Ufail(j, k) <= umin(j)
                FTCS(FT).Ufail(j, k) = umin(j);
            end
        end
        
        %% Process simulation with ODE
        [tsim, x] = ode45(@(x, u) CSTR(FTCS(FT).X(:, k), FTCS(FT).Ufail(:, k)) , [0 Ts], FTCS(FT).X(:, k), ode_options);
        FTCS(FT).X(:, k+1) = x(end, :)' + v(:, k);  % State evolution with measure noise

        % Natural state limits (maybe should be erased)
        for j = 1:nu
            if FTCS(FT).X(j, k) >= xmax(j)
                FTCS(FT).X(j, k) = xmax(j);
            elseif FTCS(FT).X(j, k) <= xmin(j)
                FTCS(FT).X(j, k) = xmin(j);
            end
        end
        
        FTCS(FT).Y(:, k) = Cd*FTCS(FT).X(:, k);     % Discrete-time output        

        %% Sensor fault income
        FTCS(FT).Yfails(:, k) = [0; 0; 0];
        FTCS(FT).Yfail(:, k) = FTCS(FT).Y(:, k);
        
        if tk > 65 && tk < 75
            FTCS(FT).Yfails(:, k) = [0; 0; Fail_S3-Fail_S3*(exp(-3*(tk-65)/1))];
            FTCS(FT).Yfail(:, k) = FTCS(FT).Y(:, k) + [0; 0; Fail_S3-Fail_S3*(exp(-3*(tk-65)/1))];
        elseif tk >= 75 && tk < 77
            FTCS(FT).Yfails(:, k) = [0; 0; Fail_S3*(exp(-5*(tk-75)/1))];
            FTCS(FT).Yfail(:, k) = FTCS(FT).Y(:, k) + [0; 0; Fail_S3*(exp(-5*(tk-75)/1))];
        end

        if tk > 50 && tk < 61
            FTCS(FT).Yfails(:, k) = [Fail_S1-Fail_S1*(exp(-5*(tk-50)/1)); 0; 0];
            FTCS(FT).Yfail(:, k) = FTCS(FT).Y(:, k) + [Fail_S1-Fail_S1*(exp(-5*(tk-50)/1)); 0; 0];
        end
        
        %% membership
        FTCS(FT).mu_fuzzy(:, k) = membership(FTCS(FT).Y(:, k), V_min, V_mid, V_max, Tr_min, Tr_mid, Tr_max);

        %% LPV-RUIO 1
        RUIO(1).Phi(:, k+1) = zeros(N, 1);
        RUIO(1).X(:, k) = zeros(nx, 1);
        RUIO(1).Fact(k) = zeros(1, 1);
        for i = 1:M
            RUIO(1).Phi(:, k+1) = RUIO(1).Phi(:, k+1) + FTCS(FT).mu_fuzzy(i, k)*(RUIO(1).O(i).K*RUIO(1).Phi(:, k) + RUIO(1).O(i).L_ast*FTCS(FT).Yfail(:, k) + RUIO(1).O(i).B_bar_1*FTCS(FT).U(:, k) + RUIO(1).O(i).delta_bar_1);
            RUIO(1).X(:, k) = RUIO(1).X(:, k) + FTCS(FT).mu_fuzzy(i, k)*(RUIO(1).O(i).T*[RUIO(1).Phi(:, k); RUIO(1).O(i).U_1*FTCS(FT).Yfail(:, k)-RUIO(1).O(i).U_1*RUIO(1).O(i).C_tilde_1*RUIO(1).Phi(:, k)]);

            RUIO(1).Fact(k) = RUIO(1).Fact(k) ...
                                        + FTCS(FT).mu_fuzzy(i, k)*(RUIO(1).O(i).U_1*(FTCS(FT).X(:, k+1) ...
                                        - RUIO(1).O(i).C_tilde_1*RUIO(1).Phi(:, k+1)) ...
                                        + RUIO(1).O(i).A_bar_22*RUIO(1).O(i).U_1*(RUIO(1).O(i).C_tilde_1*RUIO(1).Phi(:, k) ...
                                        - FTCS(FT).Yfail(:, k)) - RUIO(1).O(i).A_bar_21*RUIO(1).Phi(:, k) ...
                                        - RUIO(1).O(i).B_bar_2*FTCS(FT).U(:, k) -  RUIO(1).O(i).delta_bar_2);
        end

        % Error norm 1
        RUIO(1).error(k) = sqrt((RUIO(1).X(1, k)-FTCS(FT).Yfail(1, k))^2 + (RUIO(1).X(2, k)-FTCS(FT).Yfail(2, k))^2 + (RUIO(1).X(3, k)-FTCS(FT).Yfail(3, k))^2);

        if RUIO(1).error(k) > threshold(1, k)
            RUIO(1).FQ(k) = true;
        else
            RUIO(1).FQ(k) = false;
        end

        %% LPV-RUIO 2
        RUIO(2).Phi(:, k+1) = zeros(N, 1);
        RUIO(2).X(:, k) = zeros(nx, 1);
        RUIO(2).Fact(k) = zeros(1, 1);
        for i = 1:M
            RUIO(2).Phi(:, k+1) = RUIO(2).Phi(:, k+1) ...
                                               + FTCS(FT).mu_fuzzy(i, k)*(RUIO(2).O(i).K*RUIO(2).Phi(:, k) ...
                                               + RUIO(2).O(i).L_ast*FTCS(FT).Yfail(:, k) ...
                                               + RUIO(2).O(i).B_bar_1*FTCS(FT).U(:, k) + RUIO(2).O(i).delta_bar_1);
            RUIO(2).X(:, k) = RUIO(2).X(:, k) ...
                                        + FTCS(FT).mu_fuzzy(i, k)*(RUIO(2).O(i).T*[RUIO(2).Phi(:, k); ...
                                        RUIO(2).O(i).U_1*FTCS(FT).Yfail(:, k)-RUIO(2).O(i).U_1*RUIO(2).O(i).C_tilde_1*RUIO(2).Phi(:, k)]);

            RUIO(2).Fact(k) = RUIO(2).Fact(k) ...
                                        + FTCS(FT).mu_fuzzy(i, k)*(RUIO(2).O(i).U_1*(FTCS(FT).X(:, k+1) ...
                                        - RUIO(2).O(i).C_tilde_1*RUIO(2).Phi(:, k+1)) ...
                                        + RUIO(2).O(i).A_bar_22*RUIO(2).O(i).U_1*(RUIO(2).O(i).C_tilde_1*RUIO(2).Phi(:, k) ...
                                        - FTCS(FT).Yfail(:, k)) - RUIO(2).O(i).A_bar_21*RUIO(2).Phi(:, k) ...
                                        - RUIO(2).O(i).B_bar_2*FTCS(FT).U(:, k) -  RUIO(2).O(i).delta_bar_2);
        end

        % Error norm 2
        RUIO(2).error(k) = sqrt((RUIO(2).X(1, k)-FTCS(FT).Yfail(1, k))^2 ...
                                              + (RUIO(2).X(2, k)-FTCS(FT).Yfail(2, k))^2 ...
                                              + (RUIO(2).X(3, k)-FTCS(FT).Yfail(3, k))^2);

        if RUIO(2).error(k) > threshold(2, k)
            RUIO(2).FQ(k) = true;
        else
            RUIO(2).FQ(k) = false;
        end
        
        %% LPV-UIOO 1
        UIOO(1).Ymon(:, k) = UIOO(1).T2*FTCS(FT).Yfail(:, k);
        UIOO(1).Z(:, k+1) = zeros(nx, 1);      
        for i = 1:M
            UIOO(1).Z(:, k+1) = UIOO(1).Z(:, k+1) ...
                                            + FTCS(FT).mu_fuzzy(i, k)*(UIOO(1).O(i).N*UIOO(1).Z(:, k) ...
                                            + UIOO(1).O(i).L*UIOO(1).Ymon(:, k) ...
                                            + UIOO(1).O(i).G*FTCS(FT).U(:, k) + UIOO(1).O(i).Tg);
        end

        UIOO(1).X(:, k) = UIOO(1).Z(:, k) - UIOO(1).E*UIOO(1).Ymon(:, k);
        UIOO(1).Y(:, k) = Cd*UIOO(1).X(:, k);

        % Residue 1
        UIOO(1).res(:, k) = FTCS(FT).Yfail(:, k) - UIOO(1).Y(:, k);

        % Error norm 1
        UIOO(1).error(k) = sqrt(UIOO(1).res(1, k)^2);

        if UIOO(1).error(k) > threshold(3, k)
            UIOO(1).FO(k) = true;
        else
            UIOO(1).FO(k) = false;
        end        

        %% LPV-UIOO 2
        UIOO(2).Ymon(:, k) = UIOO(2).T2*FTCS(FT).Yfail(:, k);
        UIOO(2).Z(:, k+1) = zeros(nx, 1);      
        for i = 1:M
            UIOO(2).Z(:, k+1) = UIOO(2).Z(:, k+1) ...
                                            + FTCS(FT).mu_fuzzy(i, k)*(UIOO(2).O(i).N*UIOO(2).Z(:, k) ...
                                            + UIOO(2).O(i).L*UIOO(2).Ymon(:, k) ...
                                            + UIOO(2).O(i).G*FTCS(FT).U(:, k) + UIOO(2).O(i).Tg);
        end

        UIOO(2).X(:, k) = UIOO(2).Z(:, k) - UIOO(2).E*UIOO(2).Ymon(:, k);
        UIOO(2).Y(:, k) = Cd*UIOO(2).X(:, k);

        % Residue 2
        UIOO(2).res(:, k) = FTCS(FT).Yfail(:, k) - UIOO(2).Y(:, k);

        % Error norm 2
        UIOO(2).error(k) = sqrt(UIOO(2).res(3, k)^2);

        if UIOO(2).error(k) > threshold(4, k)
            UIOO(2).FO(k) = true;
        else
            UIOO(2).FO(k) = false;
        end        
        
        %% Actuator fault estimation
        % Actuator fault 1
        if RUIO(1).FQ(k) && ~RUIO(2).FQ(k) && UIOO(1).FO(k) && UIOO(2).FO(k)
            if RUIO(2).delay > 1
                RUIO(2).Fact(k) = RUIO(2).Fact(k);
            else
                RUIO(2).delay = RUIO(2).delay + 1;
                RUIO(2).Fact(k) = 0;
            end
        else
            RUIO(2).delay = 0;
            RUIO(2).Fact(k) = 0;
        end
        
        % Actuator fault 2
        if ~RUIO(1).FQ(k) && RUIO(2).FQ(k) && UIOO(1).FO(k) && ~UIOO(2).FO(k)
            if RUIO(1).delay > 1
                RUIO(1).Fact(k) = RUIO(1).Fact(k);
            else
                RUIO(1).delay = RUIO(1).delay + 1;
                RUIO(1).Fact(k) = 0;
            end
        else
            RUIO(1).delay = 0;
            RUIO(1).Fact(k) = 0;       
        end
        
        %% Sensor fault estimation
        % Sensor fault 1
        if RUIO(1).FQ(k) && RUIO(2).FQ(k) && UIOO(1).FO(k) && ~UIOO(2).FO(k)
            UIOO(1).Fsen(k) = UIOO(2).res(1, k);
        else
            UIOO(1).Fsen(k) = zeros(size(UIOO(2).res(1, k)));
        end

        % Sensor fault 2
        if RUIO(1).FQ(k) && RUIO(2).FQ(k) && ~UIOO(1).FO(k) && UIOO(2).FO(k)
            UIOO(2).Fsen(k) = UIOO(1).res(3, k);
        else
            UIOO(2).Fsen(k) = zeros(size(UIOO(1).res(3, k)));
        end        
        
        % If FT-MPC is enabled
        if FT == FTC_ON
            FTCS(FT).Uff(:, k) = [RUIO(1).Fact(k); RUIO(2).Fact(k)];
            FTCS(FT).Y_hat(:, k) = [FTCS(FT).Yfail(1, k)-UIOO(1).Fsen(k); FTCS(FT).Yfail(2, k); FTCS(FT).Yfail(3, k)-UIOO(2).Fsen(k)];
        else
            FTCS(FT).Uff(:, k) = [0; 0];
            FTCS(FT).Y_hat(:, k) = FTCS(FT).Yfail(:, k);            
        end
                
        % MHE horizon update
        FTCS(FT).X_MHE = [FTCS(FT).X_MHE(:, 2:end) FTCS(FT).Y_hat(:, k)];                         % Update sensor fault compensation
        FTCS(FT).U_MHE = [FTCS(FT).U_MHE(:, 2:end) FTCS(FT).U(:, k)+FTCS(FT).Uff(:, k)]; % Update actuator fault compensation
        
        % Store data for plot
        FTCS(FT).RUIO(1).error(k) = RUIO(1).error(k);
        FTCS(FT).RUIO(1).Fact(k) = RUIO(1).Fact(k);
        FTCS(FT).RUIO(2).error(k) = RUIO(2).error(k);
        FTCS(FT).RUIO(2).Fact(k) = RUIO(2).Fact(k);
        
        FTCS(FT).UIOO(1).error(k) = UIOO(1).error(k);
        FTCS(FT).UIOO(1).Fsen(k) = UIOO(1).Fsen(k);
        FTCS(FT).UIOO(2).error(k) = UIOO(2).error(k);
        FTCS(FT).UIOO(2).Fsen(k) = UIOO(2).Fsen(k);        
    
	end
        
end

yalmip('clear')
disp('Saving...')
save FTCS.mat

%%
run enPlotCSTR