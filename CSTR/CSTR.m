%% CSTR model
function out = CSTR(x, u)
	%
	% CSTR function
	%
	% x: State matrix [V CA T]
	% u: Inputs matrix [qs qc]
	% out: Output matrix [dV dCA dT]

	% States
	V = x(1);
	CA = x(2);
	T = x(3);

	% Inputs
	qs = u(1);
	qc = u(2);
	
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
    qe = 100;                     % Feed flow rate [l/min]
	
	% Matter balance -> Volume
	dV = qe - qs;
	% Matter balance -> Concentration 
	dCA = (qe/V)*(CAe-CA) - k0*CA*exp(-E_R/T);
	% Energy balance -> Temperature
	dT = (qe/V)*(Te-T) - k1*CA*exp(-E_R/T) + k2*(qc/V)*(1-exp(-k3/qc))*(Tce-T);

 	out = [dV; dCA; dT];
end