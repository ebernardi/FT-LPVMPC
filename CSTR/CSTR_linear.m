%% CSTR Model
syms qs qc V CA T

% Parameters
E_R = 1e4;                  % Activation energy [°K]
Te = 350;                    % Feed temperature [°K]
Tce = 350;                  % Inlet coolant temp. [°K]
dH = -2e5;                  % Heat of reaction [cal/mol]
Cp = 1;                        % Specific heats [cal/g °K]
rho = 1e3;                   % Liquid densities [g/l]
CAe = 1;                      % Feed concentration [mol/l]
ha = 7e5;                     % Heat transfer term [cal/min °K]
k0 = 7.2e10;                % Reaction rate constant [l/min]
k1 = dH*k0/(rho*Cp);
k2 = rho*Cp/(rho*Cp);
k3 = ha/(rho*Cp);
k4 = 10;                       % Valve constant [l/min m^3/2]
qe = 100;                     % Feed flow rate [l/min]

%% Non-Linear model
system = [(qe - qs);
                  ((qe/V)*(CAe - CA) - k0*CA*exp(-E_R/T));
                  ((qe/V)*(Te - T) - k1*CA*exp(-E_R/T) + k2*(qc/V)*(1 - exp(-k3/qc))*(Tce - T))];

states = [V CA T];
outputs = [V CA T];
inputs = [qs qc];
nx = length(states); nu = length(inputs); ny = length(outputs);
C = eye(ny, nx);       % Output matrix
D = zeros(ny, nu);    % Input/Output matrix

%% Linealization
% Symbolic matrices
A_sym = jacobian(system, states);
B_sym = jacobian(system, inputs);

% Reactor temperature
% Tr = -E_R/log(-(q*(Ca-CAe))/(k0*Ca*Vr));

% Output concentration
Ca = CAe/(1+(k0*(Vr/qe)*exp(-E_R/Tr)));

% States
X_lin = [Vr; Ca; Tr];

% Output flow rate
Qs = double(solve(qe - qs));

% Coolant flow rate
Qc = double(solve(qe/Vr*(Te-Tr) - k1*Ca*exp(-E_R/Tr) + k2*(qc/Vr)*(1-exp(-k3/qc))*(Tce-Tr) == 0));

% Inputs
U_lin = [Qs; Qc];

% Linear systems matrices
A = subs(A_sym, {CA, T, qs, qc, V}, {Ca, Tr, Qs, Qc, Vr});
B = subs(B_sym, {CA, T, qs, qc, V}, {Ca, Tr, Qs, Qc, Vr});
A = double(A);
B = double(B);

f = subs(system, {CA, T, qs, qc, V}, {Ca, Tr, Qs, Qc, Vr});
f = double(f);

% Constant term
delta = f - (A*X_lin+B*U_lin);

% Euler discretization method
Ad = (A*Ts) + eye(nx); Bd = B*Ts; Cd = C; Dd = D; deltad = delta*Ts;