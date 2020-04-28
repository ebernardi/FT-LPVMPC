%% Heat Exchanger model
function out = HE(x, u)
	%
	% Heat exchanger function
	%
	% x: state matrix [Theta_1s Theta_2s Theta_p]
	% u: input matrix [q1 q2]
	% out: output matrix [dTheta_1s dTheta_2s dTheta_p]

	% States
	Theta_1s = x(1);
	Theta_2s = x(2);
	Theta_p = x(3);

	% Inputs
	q1 = u(1);
	q2 = u(2);
    
    % Fluid 1 is cold process
    % Fluid 2 is hot stream
    	
    % Parameters
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
	
	% Energy balance -> Process temperature
	dTheta_1s = (q1*Rho1*Cp_1*(Theta_1e-Theta_1s) - a*h_1*(Theta_1s-Theta_p))/(Rho1*V_1*Cp_1);
	% Energy balance -> Wall temperature
	dTheta_2s = (q2*Rho2*Cp_2*(Theta_2e-Theta_2s) + a*h_2*(Theta_p-Theta_2s))/(Rho2*V_2*Cp_2);
	% Energy balance -> Hot stream temperature
	dTheta_p = (a*h_1*(Theta_1s-Theta_p) - a*h_2*(Theta_p-Theta_2s))/(Rhop*Cp_p*V_p);

 	out = [dTheta_1s; dTheta_2s; dTheta_p];
end
