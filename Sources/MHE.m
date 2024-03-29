%% MHE
yalmip('clear')

%% MHE definition
x = sdpvar(nx*ones(1, N_MHE+1), ones(1, N_MHE+1));
xru = sdpvar(nx*ones(1, N_MHE+1), ones(1, N_MHE+1));
nux = sdpvar(nx*ones(1, N_MHE+1), ones(1, N_MHE+1));
u = sdpvar(nu*ones(1, N_MHE), ones(1, N_MHE));
mu = sdpvar(M, 1);
mu_last = sdpvar(M, 1);
numu = sdpvar(M, 1);
e = sdpvar(nx*ones(1, N_MHE+1), ones(1, N_MHE+1));

objective = 0; constraints = [];

A = [sys.Ad]*kron(mu, eye(nx));
B = [sys.Bd]*kron(mu, eye(nu));
deltad = [sys.deltad]*mu;

Qe = 1e4; Qnux = 8e2; Qnumu = 1e1; % Qe = 1e4; Qnux = 8e2; Qnumu = 1e1;

% Stage constraints and objective
for k = 1:N_MHE
    % Objective
    objective = objective + (e{k})'*Qe*(e{k}) + (nux{k})'*Qnux*(nux{k}); 

    % Dynamic constraint    
    constraints = [constraints, e{k+1} == xru{k+1} - (A*x{k} + B*u{k} + deltad)];
    constraints = [constraints, xru{k} == x{k} + nux{k}];   
end
objective = objective + (numu)'*Qnumu*(numu);

constraints = [constraints, mu == mu_last + numu];
constraints = [constraints, sum(mu) == 1];
constraints = [constraints, zeros(M, 1) <= mu <= ones(M, 1)];

% Defining the parameters in, and the solution
parameters = {[x{:}], [u{:}], mu_last}; % 
solution = {mu};

% Options for Optimizer  
options = sdpsettings('solver', 'gurobi');
options.verbose = 0;                             % 0 to none, 2 to debug, 2+ for more

mhe = optimizer(constraints, objective, options, parameters, solution);