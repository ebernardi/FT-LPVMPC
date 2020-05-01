%% HE Model
syms rho1 rho2 rhop Cp1 Cp2 Cpp Ar h1 h2 V1 V2 Vp theta1e theta2e q1 q2 theta1s theta2s thetap
sys = struct;

% Fluid 1 is cold process
% Fluid 2 is hot stream

%% Parameters
Rho1 = 1;               % Fluid 1 density (kg/l)
Rho2 = 1;               % Fluid 2 density (kg/l)
Rhop = 7.874;        % Wall density (kg/l)
Cp_1 = 1000;         % Heat capacity fluid 1 (cal/kg K)
Cp_2 = 1000;         % Heat capacity fluid 2 (cal/kg K)
Cp_p = 1075.53;    % Wall specific heat (cal/kg K)
a = 0.881;              % Area HE (m^2)
h_1 = 32374;         % Heat transfer fluid 1 (cal/min K m^2)
h_2 = 14716.6667;% Heat transfer fluid 2 (cal/min K m^2)
V_1 = 16;               % Tube Volume (l)
V_2 = 2.11;            % Case Volume (l)
V_p = 1.19;            % Wall Volume (l)
Theta_1e = 480;    % Input Temp. fluid 1 (K)
Theta_2e = 900;    % Input Temp. fluid 2 (K)

%% Non-Linear model
system = [(q1*rho1*Cp1*(theta1e-theta1s)-Ar*h1*(theta1s-thetap))/(rho1*V1*Cp1);
                  (q2*rho2*Cp2*(theta2e-theta2s)+Ar*h2*(thetap-theta2s))/(rho2*V2*Cp2);
                  (Ar*h1*(theta1s-thetap)-Ar*h2*(thetap-theta2s))/(rhop*Cpp*Vp)];

inputs = [q1 q2];
outputs = [theta1s theta2s thetap];
states = [theta1s theta2s thetap];
nx = length(states); nu = length(inputs); ny = length(outputs);
C = eye(ny, nx);       % Output matrix
D = zeros(ny, nu);    % Input/Output matrix

%% Linealization
% Symbolic matrices
A_sym = jacobian(system, states);
B_sym = jacobian(system, inputs);

for i = 1:M
    switch i
        case 1
            % System 1 (Theta_1s = Theta_1s_min [K]; Theta_2s = Theta_2s_min [K])
            Theta_1s = Theta_1s_min;
            Theta_2s = Theta_2s_min;
        case 2
            % System 2 (Theta_1s = Theta_1s_min [K]; Theta_2s = Theta_2s_mid [K])
            Theta_1s = Theta_1s_min;
            Theta_2s = Theta_2s_mid;
        case 3
            % System 3 (Theta_1s = Theta_1s_min [K]; Theta_2s = Theta_2s_max [K])
            Theta_1s = Theta_1s_min;
            Theta_2s = Theta_2s_max;
        case 4
            % System 4 (Theta_1s = Theta_1s_mid [K]; Theta_2s = Theta_2s_min [K])
            Theta_1s = Theta_1s_mid;
            Theta_2s = Theta_2s_min;
        case 5
            % System 5 (Theta_1s = Theta_1s_mid [K]; Theta_2s = Theta_2s_mid [K])
            Theta_1s = Theta_1s_mid;
            Theta_2s = Theta_2s_mid;
        case 6
            % System 6 (Theta_1s = Theta_1s_mid [K]; Theta_2s = Theta_2s_max [K])
            Theta_1s = Theta_1s_mid;
            Theta_2s = Theta_2s_max;
        case 7
            % System 7 (Theta_1s = Theta_1s_max [K]; Theta_2s = Theta_2s_min [K])
            Theta_1s = Theta_1s_max;
            Theta_2s = Theta_2s_min;
        case 8
            % System 8 (Theta_1s = Theta_1s_max [K]; Theta_2s = Theta_2s_mid [K])
            Theta_1s = Theta_1s_max;
            Theta_2s = Theta_2s_mid;
        case 9
            % System 9 (Theta_1s = Theta_1s_max [K]; Theta_2s = Theta_2s_max [K])
            Theta_1s = Theta_1s_max;
            Theta_2s = Theta_2s_max;
        otherwise
    end
    
    % Wall temperature
    Theta_p = (h_2*Theta_2s + h_1*Theta_1s) / (h_2 + h_1);

    % Linear states
    Xinit = [Theta_1s; Theta_2s; Theta_p];

    % Fluid 1 flow rate
    Q1 = ( a*h_1*(Theta_1s-Theta_p) ) / ( Rho1*Cp_1*(Theta_1e - Theta_1s) );

    % Fluid 2 flow rate
    Q2 = ( -a*h_2*(Theta_p - Theta_2s) ) / ( Rho2*Cp_2*(Theta_2e - Theta_2s) );

    % Test
    % Theta_2s = ( (h_1+h_2)*Theta_p - h_1*Theta_1s ) / h_2
    % Theta_2s = (Q2*Rho2*Cp_2*Theta_2e + (a*h_2*h_1*Theta_1s / (h_1+h_2)) )/ ...
    %                       (Q2*Rho2*Cp_2 + a*h_2 - (a*h_2^2 / (h_1+h_2)))
    % Theta_p = ( (Q1*Rho1*Cp_1 + a*h_1)*Theta_1s - Q1*Rho1*Cp_1*Theta_1e ) / (a*h_1)

    % Inputs
    Uinit = [Q1; Q2];

    % Linear systems matrices
    A = subs(A_sym, {rho1, rho2, rhop, Cp1, Cp2, Cpp, Ar, h1, h2, V1, V2, Vp, theta1e, theta2e, theta1s, theta2s, thetap, q1, q2}, ...
              {Rho1, Rho2, Rhop, Cp_1, Cp_2, Cp_p, a, h_1, h_2, V_1, V_2, V_p, Theta_1e, Theta_2e, Xinit(1), Xinit(2), Xinit(3), Uinit(1), Uinit(2)});
    B = subs(B_sym, {rho1, rho2, rhop, Cp1, Cp2, Cpp, Ar, h1, h2, V1, V2, Vp, theta1e, theta2e, theta1s, theta2s, thetap, q1, q2}, ...
              {Rho1, Rho2, Rhop, Cp_1, Cp_2, Cp_p, a, h_1, h_2, V_1, V_2, V_p, Theta_1e, Theta_2e, Xinit(1), Xinit(2), Xinit(3), Uinit(1), Uinit(2)});
    A = double(A);
    B = double(B);

    f = subs(system, {rho1, rho2, rhop, Cp1, Cp2, Cpp, Ar, h1, h2, V1, V2, Vp, theta1e, theta2e, theta1s, theta2s, thetap, q1, q2}, ...
              {Rho1, Rho2, Rhop, Cp_1, Cp_2, Cp_p, a, h_1, h_2, V_1, V_2, V_p, Theta_1e, Theta_2e, Xinit(1), Xinit(2), Xinit(3), Uinit(1), Uinit(2)});
    f = double(f);

    % Constant term
    delta = f - (A*Xinit+B*Uinit);

	% Continuous system
    sys(i).A = A; sys(i).B = B; sys(i).delta = delta;
    sys(i).C = C; sys(i).D = D;
    
    % Steady state
    sys(i).Xinit = Xinit; sys(i).Uinit = Uinit;
    
    % Euler discretization method
    sys(i).Ad = (A*Ts) + eye(nx); sys(i).Bd = B*Ts; sys(i).deltad = delta*Ts;
    
%     % Zero order holder discretization method
%     Abar = [A delta; zeros(1, nx) 0];
%     Bbar = [B; zeros(1, nu)];
%     Cbar = [C zeros(ny, 1)];
%     Dbar = D;
%     [sys(i).Abard, sys(i).Bbard, sys(i).Cbard, sys(i).Dbard] = c2dm(Abar, Bbar, Cbar, Dbar, Ts, 'zoh');
%     sys(i).Ad = sys(i).Abard(1:nx, 1:nx); sys(i).Bd = sys(i).Bbard(1:nx, 1:nu); sys(i).deltad = sys(i).Abard(1:nx, nx+1);

    % Common matrices
    sys(i).Cd = C; sys(i).Dd = D;
    
end
