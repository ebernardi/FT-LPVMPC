%% LPV MPC
yalmip('clear')

%% Controller definition
x = sdpvar(nx*ones(1, N_MPC+1), ones(1, N_MPC+1));
u = sdpvar(nu*ones(1, N_MPC), ones(1, N_MPC));
xs = sdpvar(nx, 1);
uf = sdpvar(nu, 1);
mu = sdpvar(M, 1);

% Artificial variables
xa = sdpvar(nx, 1); 
ua = sdpvar(nu, 1);

objective = 0; constraints = [];

A = [sys.Ad]*kron(mu, eye(nx));
B = [sys.Bd]*kron(mu, eye(nu));
deltad = [sys.deltad]*mu;

gamma = 10*Plqr;

% Stage constraints and objective
for k = 1:N_MPC
    % Objective
    objective = objective + (x{k}-xa)'*Qx*(x{k}-xa); 
    objective = objective + (u{k}+uf-ua)'*Ru*(u{k}+uf-ua);

    % Dynamic constraint    
    constraints = [constraints, x{k+1} == A*x{k} + B*(u{k}+uf) + deltad];

    % Box-type constraint
    constraints = [constraints, umin <= u{k}+uf <= umax];
    constraints = [constraints, xmin <= x{k} <= xmax];
end

% Terminal constraints
objective = objective + (xa-xs)'*gamma*(xa-xs);                                % Artificial variable terminal cost
objective = objective + (x{N_MPC+1}-xa)'*Plqr*(x{N_MPC+1}-xa);  % Terminal cost

constraints = [constraints, x{N_MPC+1} == xa];                               % Equality terminal constraint
constraints = [constraints, xa == A*xa + B*ua + deltad];                  % Artificial variables equilibirum condition
constraints = [constraints, Xf.A*x{N_MPC+1} <= Xf.b];                     % Inequality terminal constraint

% Defining the parameters in, and the solution
parameters = {x{1}, xs, uf, mu};
solution = {u{1}, objective};

% Options for Optimizer  
options = sdpsettings('solver', 'gurobi');
options.verbose = 0;                             % 0 to none, 2 to debug, 2+ for more

mpc = optimizer(constraints, objective, options, parameters, solution);