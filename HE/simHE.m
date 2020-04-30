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
Time = 50;                 % Simulation end time 
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
Tau = 2;                % Convergence period
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
xmin = [490; 675; 545];
xmax = [505; 715; 570];
umin = [70; 7];
umax = [130; 11];

% Wheight matrix
Qx = sys(1).Cd'*sys(1).Cd;
lambda = 0.1;
Ru = lambda*diag([1 10]);

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
    
    Theta_1s = Theta_1s_mid;     % Output fluid 1 temperature (K)
    Theta_2s = Theta_2s_mid;     % Output fluid 2 temperature (K)
    run HE_linear;
    x0 = [Theta_1s; Theta_2s; Theta_p];
    u0 = [Q1; Q2];

    Theta_1s = Theta_1s_min;     % Output fluid 1 temperature (K)
    Theta_2s = Theta_2s_max;     % Output fluid 2 temperature (K)
    run HE_linear;
    xsp = [Theta_1s; Theta_2s; Theta_p];    
    
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
    FTCS(FT).mu_fuzzy(:, k) = membership(FTCS(FT).Y(:, k), Theta_1s_min, Theta_1s_mid, Theta_1s_max, Theta_2s_min, Theta_2s_mid, Theta_2s_max);
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
            Theta_1s = Theta_1s_mid;     % Output fluid 1 temperature (K)
            Theta_2s = Theta_2s_mid-1;     % Output fluid 2 temperature (K)
            run HE_linear;
            xsp = [Theta_1s; Theta_2s; Theta_p];
        elseif tk == 17
            Theta_1s = Theta_1s_max-1;     % Output fluid 1 temperature (K)
            Theta_2s = Theta_2s_mid+1;     % Output fluid 2 temperature (K)
            run HE_linear;
            xsp = [Theta_1s; Theta_2s; Theta_p];
        end
        FTCS(FT).Xsp(:, k) = xsp;
        
        if k > 1
            FTCS(FT).mu_mhe(:, k) = FTCS(FT).mu_mhe(:, k-1);
            FTCS(FT).Y(:, k) = FTCS(FT).Y(:, k-1);
            FTCS(FT).Uff(:, k) = FTCS(FT).Uff(:, k-1);
        end
        
%         FTCS(FT).mu_fuzzy(:, k) = membership(FTCS(FT).Y(:, k), Theta_1s_min, Theta_1s_mid, Theta_1s_max, Theta_2s_min, Theta_2s_mid, Theta_2s_max);
      
        t_tic = tic ;   % to get time evaluated 
        
        [sol, diag] = mhe{FTCS(FT).X_MHE, FTCS(FT).U_MHE, FTCS(FT).mu_mhe(:, k)};
        if diag
            msg = ['Infeasible MHE at t = ', num2str(tk)];
            disp(msg)
            return;
        end
        FTCS(FT).mu_mhe(:, k) = sol;
%         FTCS(FT).mu_mhe(:, k) = [0; 0; 0; 0; 0; 0; 0; 0; 1];
                
        [sol, diag] = mpc{FTCS(FT).Y(:, k), FTCS(FT).Xsp(:, k), FTCS(FT).Uff(:, k), FTCS(FT).mu_mhe(:, k)};
        if diag
            msg = ['Infeasible MPC at t = ', num2str(tk)];
            disp(msg)
            return;
        end
        FTCS(FT).U(:, k) = sol{1}; FTCS(FT).Obj(k) = sol{2};
        
        t_tic = toc(t_tic) ;              % get time elapsed
        FTCS(FT).elapsed_time(k) = t_tic ;   % store the time elapsed for the run
        
   %% Actuator fault income
        % No fault
        FTCS(FT).Ufails(:, k) = [0; 0];
        FTCS(FT).Ufail(:, k) = FTCS(FT).U(:, k);
        FTCS(FT).Umax(:, k) = umax;
        FTCS(FT).Umin(:, k) = umin;
 
        % Q1 fault
        if tk > 10 && tk < 15
            FTCS(FT).Ufails(:, k) = [Fail_Q1; 0];
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
            if FTCS(FT).Ufail(j, k) >= FTCS(FT).Umax(j, k)
                FTCS(FT).Ufail(j, k) = FTCS(FT).Umax(j, k);
            elseif FTCS(FT).Ufail(j, k) <= FTCS(FT).Umin(j, k)
                FTCS(FT).Ufail(j, k) = FTCS(FT).Umin(j, k);
            end
        end        
        
        %% Process simulation with ODE
        [tsim, x] = ode45(@(x, u) HE(FTCS(FT).X(:, k), FTCS(FT).Ufail(:, k)) , [0 Ts], FTCS(FT).X(:, k), ode_options);
        FTCS(FT).X(:, k+1) = x(end, :)';

        % Natural state limits (maybe should be erased)
        for j = 1:nu
            if FTCS(FT).X(j, k) >= xmax(j)
                FTCS(FT).X(j, k) = xmax(j);
            elseif FTCS(FT).X(j, k) <= xmin(j)
                FTCS(FT).X(j, k) = xmin(j);
            end
        end
        
        FTCS(FT).Y(:, k) = C*FTCS(FT).X(:, k);                   % Discrete-time output

        %% Sensor fault income
        FTCS(FT).Yfail(:, k) = FTCS(FT).Y(:, k);
        
        if tk > 35 && tk < 45
            FTCS(FT).Yfail(:, k) = FTCS(FT).Y(:, k) + [0; Fail_S2-Fail_S2*(exp(-5*(tk-35)/1)); 0];
        elseif tk >= 45 && tk < 47
            FTCS(FT).Yfail(:, k) = FTCS(FT).Y(:, k) + [0; +Fail_S2*(exp(-8*(tk-45)/1)); 0];
        end

        if tk > 50 && tk < 61
            FTCS(FT).Yfail(:, k) = FTCS(FT).Y(:, k) + [Fail_S1-Fail_S1*(exp(-5*(tk-50)/1)); 0; 0];
        end
        
        %% membership
        FTCS(FT).mu_fuzzy(:, k) = membership(FTCS(FT).Y(:, k), Theta_1s_min, Theta_1s_mid, Theta_1s_max, Theta_2s_min, Theta_2s_mid, Theta_2s_max);

        %% LPV-RUIO 1
        RUIO(1).Phi(:, k+1) = zeros(N, 1);
        RUIO(1).X(:, k) = zeros(nx, 1);
        RUIO(1).Fact(k) = zeros(1, 1);
        for i = 1:M
            RUIO(1).Phi(:, k+1) = RUIO(1).Phi(:, k+1) + FTCS(FT).mu_fuzzy(i, k)*(RUIO(1).O(i).K*RUIO(1).Phi(:, k) + RUIO(1).O(i).L_ast*FTCS(FT).Yfail(:, k) + RUIO(1).O(i).B_bar_1*FTCS(FT).U(:, k) + RUIO(1).O(i).delta_bar_1);
            RUIO(1).X(:, k) = RUIO(1).X(:, k) + FTCS(FT).mu_fuzzy(i, k)*(RUIO(1).O(i).T*[RUIO(1).Phi(:, k); RUIO(1).O(i).U_1*FTCS(FT).Yfail(:, k)-RUIO(1).O(i).U_1*RUIO(1).O(i).C_tilde_1*RUIO(1).Phi(:, k)]);

            RUIO(1).Fact(k) = RUIO(1).Fact(k) + FTCS(FT).mu_fuzzy(i, k)*(RUIO(1).O(i).U_1*(FTCS(FT).X(:, k+1) - RUIO(1).O(i).C_tilde_1*RUIO(1).Phi(:, k+1)) + RUIO(1).O(i).A_bar_22*RUIO(1).O(i).U_1*(RUIO(1).O(i).C_tilde_1*RUIO(1).Phi(:, k) - FTCS(FT).Yfail(:, k)) - RUIO(1).O(i).A_bar_21*RUIO(1).Phi(:, k) - RUIO(1).O(i).B_bar_2*FTCS(FT).U(:, k) -  RUIO(1).O(i).delta_bar_2);
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
            RUIO(2).Phi(:, k+1) = RUIO(2).Phi(:, k+1) + FTCS(FT).mu_fuzzy(i, k)*(RUIO(2).O(i).K*RUIO(2).Phi(:, k) + RUIO(2).O(i).L_ast*FTCS(FT).Yfail(:, k) + RUIO(2).O(i).B_bar_1*FTCS(FT).U(:, k) + RUIO(2).O(i).delta_bar_1);
            RUIO(2).X(:, k) = RUIO(2).X(:, k) + FTCS(FT).mu_fuzzy(i, k)*(RUIO(2).O(i).T*[RUIO(2).Phi(:, k); RUIO(2).O(i).U_1*FTCS(FT).Yfail(:, k)-RUIO(2).O(i).U_1*RUIO(2).O(i).C_tilde_1*RUIO(2).Phi(:, k)]);

            RUIO(2).Fact(k) = RUIO(2).Fact(k) + FTCS(FT).mu_fuzzy(i, k)*(RUIO(2).O(i).U_1*(FTCS(FT).X(:, k+1) - RUIO(2).O(i).C_tilde_1*RUIO(2).Phi(:, k+1)) + RUIO(2).O(i).A_bar_22*RUIO(2).O(i).U_1*(RUIO(2).O(i).C_tilde_1*RUIO(2).Phi(:, k) - FTCS(FT).Yfail(:, k)) - RUIO(2).O(i).A_bar_21*RUIO(2).Phi(:, k) - RUIO(2).O(i).B_bar_2*FTCS(FT).U(:, k) -  RUIO(2).O(i).delta_bar_2);
        end

        % Error norm 2
        RUIO(2).error(k) = sqrt((RUIO(2).X(1, k)-FTCS(FT).Yfail(1, k))^2 + (RUIO(2).X(2, k)-FTCS(FT).Yfail(2, k))^2 + (RUIO(2).X(3, k)-FTCS(FT).Yfail(3, k))^2);

        if RUIO(2).error(k) > threshold(2, k)
            RUIO(2).FQ(k) = true;
        else
            RUIO(2).FQ(k) = false;
        end
        
        %% LPV-UIOO 1
        UIOO(1).Ymon(:, k) = UIOO(1).T2*FTCS(FT).Yfail(:, k);
        UIOO(1).Z(:, k+1) = zeros(nx, 1);      
        for i = 1:M
            UIOO(1).Z(:, k+1) = UIOO(1).Z(:, k+1) + FTCS(FT).mu_fuzzy(i, k)*(UIOO(1).O(i).N*UIOO(1).Z(:, k) + UIOO(1).O(i).L*UIOO(1).Ymon(:, k) + UIOO(1).O(i).G*FTCS(FT).U(:, k) + UIOO(1).O(i).Tg);
        end

        UIOO(1).X(:, k) = UIOO(1).Z(:, k) - UIOO(1).E*UIOO(1).Ymon(:, k);
        UIOO(1).Y(:, k) = C*UIOO(1).X(:, k);

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
            UIOO(2).Z(:, k+1) = UIOO(2).Z(:, k+1) + FTCS(FT).mu_fuzzy(i, k)*(UIOO(2).O(i).N*UIOO(2).Z(:, k) + UIOO(2).O(i).L*UIOO(2).Ymon(:, k) + UIOO(2).O(i).G*FTCS(FT).U(:, k) + UIOO(2).O(i).Tg);
        end

        UIOO(2).X(:, k) = UIOO(2).Z(:, k) - UIOO(2).E*UIOO(2).Ymon(:, k);
        UIOO(2).Y(:, k) = Cd*UIOO(2).X(:, k);

        % Residue 2
        UIOO(2).res(:, k) = FTCS(FT).Yfail(:, k) - UIOO(2).Y(:, k);

        % Error norm 2
        UIOO(2).error(k) = sqrt(UIOO(2).res(2, k)^2);

        if UIOO(2).error(k) > threshold(4, k)
            UIOO(2).FO(k) = true;
        else
            UIOO(2).FO(k) = false;
        end        
        
        %% Actuator fault estimation
        if RUIO(1).FQ(k) && ~RUIO(2).FQ(k)% && ~UIOO(1).FO(k) && UIOO(2).FO(k) % Actuator fault 1
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
        if ~RUIO(1).FQ(k) && RUIO(2).FQ(k)% && UIOO(1).FO(k) && ~UIOO(2).FO(k) % Actuator fault 2
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
        
        % If FT-MPC is enabled
        if FT == FTC_ON
            FTCS(FT).Uff(:, k) = [RUIO(1).Fact(k); RUIO(2).Fact(k)];
        else
            FTCS(FT).Uff(:, k) = [0; 0];
        end
        
        % MHE horizon update
        FTCS(FT).X_MHE = [FTCS(FT).X_MHE(:, 2:end) FTCS(FT).X(:, k)];
        FTCS(FT).U_MHE = [FTCS(FT).U_MHE(:, 2:end) FTCS(FT).U(:, k)+FTCS(FT).Uff(:, k)]; % Update fault compensation
        
        % Save data to plot
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

%%
run enPlotHE