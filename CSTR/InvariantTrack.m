function [Oinft, M] = InvariantTrack(A, B, C, D, Klqr, X, U)

[nx, nu] = size(B);
ny = length(C(:, 1));

if isempty(X)
    nrx = 0;
else
	nrx = length(X(:, 1));
end
if isempty(U)
	nru = 0;
else
	nru = length(U(:, 1));
end

%% Calculos preliminares
E = [A-eye(nx, nx) B; C D];
F = [zeros(nx, ny); eye(ny, ny)];
[Ue, Se, Ve] = svd(E);
Vt = null(Ve');
if(isempty(Vt))
	Vt = zeros(length(Ve(:, 1)), 1);
end
Ut = null(Ue');
if(isempty(Ut))
	Ut = zeros(length(Ue(:, 1)), 1);
end
rango = rank(E);
if(rango < (nx+ny))
%     UL = zeros(nx+ny, nx+ny-rango);
    FUL = F'*Ut;
    FULt = null(FUL');
    G = FULt;
    G = ones(length(FULt(:, 1)), 1);
    if(isempty(G))
        G = zeros(length(FULt(:, 1)), 1);
    end
else
	G = eye(ny, ny);
end
M = Ue'*F*G;
M = Se\M;
M=Ve*M;
if(rango < (nx+nu))
    M = [M Vt];
    N = [G zeros(ny, nx+nu-rango)];
else
	N = G;
end
Mx = M(1:nx, :);
Mu = M(nx+1:nx+nu, :);
L = [-Klqr eye(nu)]*M;

display('Conjunto Terminal para Seguimiento');

lambda = 0.99;  % Arbitrariamente cercano a 1
Acl = [A+B*Klqr B*L; zeros(nu, nx) eye(nu)];    % Matriz A del sistema extendido de lazo cerrado

% Conjunto de restricciones del sistema extendido
W = zeros(2*nrx+2*nru, nx+nu+1); 
W(1:nrx, :) = [X(:, 1:end-1), zeros(nrx, nu), X(:, end)];
W(nrx+1:nrx+nru, :) = [U(:, 1:end-1)*Klqr, U(:, 1:end-1)*L, U(:, end)];
W(nrx+nru+1:2*nrx+nru, :) = [zeros(nrx, nx), X(:, 1:end-1)*Mx, lambda*X(:, end)];
W(2*nrx+nru+1:end, :) = [zeros(nru, nx), U(:, 1:end-1)*Mu, lambda*U(:, end)];

maxiteration = 100;
Oi = W;

maxiteration = maxiteration+1;

for i = 2:maxiteration
    fprintf('Iteracion %d \n', i-1);
    Oit = Oi;
    Oi = remred([Oit; [W(:, 1:end-1)*Acl^(i-1), W(:, end)]]);
    if issubset(Oit, Oi)
        fprintf('Invariante con numero de pasos = %d \n', i);
        break;
    end
end

% Oi = proyx(Oi, nx);
Oinft = Oi;