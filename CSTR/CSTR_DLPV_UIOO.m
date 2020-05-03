% Type of observer(N° observer). Sub-observer(N° sub-observer). Matrix 
UIOO = struct;

% Dimension of system matrices
n = 3;
p = 2;

for j = 1:N
    UIOO(j).H = zeros(size(sys(1).C));
    if j == 1
        UIOO(j).H(3, :) = sys(1).C(3, :);
        UIOO(j).alpha = 1;
    else
        UIOO(j).H(1, :) = sys(1).C(1, :);
        UIOO(j).alpha = 1;
    end
    UIOO(j).T2 = null(UIOO(j).H, 'r')';
    UIOO(j).J = UIOO(j).T2*sys(1).C;

    yalmip('clear');
    X = sdpvar(n);
    S = sdpvar(n, p);

    const = [];
    const = [const, X >= 0];

    for i = 1:M
        if j == 1
            UIOO(j).O(i).F = sys(i).Bd(:, 1);
        else
            UIOO(j).O(i).F = sys(i).Bd;
        end

        UIOO(j).O(i).W = sdpvar(n, p);

        UIOO(j).O(i).LMI = [2*UIOO(j).alpha*X, (sys(i).Ad'*X+sys(i).Ad'*UIOO(j).J'*S'-UIOO(j).J'*UIOO(j).O(i).W');
                                      (X*sys(i).Ad+S*UIOO(j).J*sys(i).Ad-UIOO(j).O(i).W*UIOO(j).J) 2*UIOO(j).alpha*X];

        const = [const, UIOO(j).O(i).LMI <= 0, (X+S*UIOO(j).J)*UIOO(j).O(i).F == 0];

    end

    diagnostics = optimize(const);

    string = yalmiperror(diagnostics.problem);
    if diagnostics.problem == 0
        disp('Factible!')
    else
        disp('Infactible!')
        return
    end

    % Matrices to be determinated
    X = double(X);
    S = double(S);

    UIOO(j).E = X\S;
    UIOO(j).T1 = (eye(n) + UIOO(j).E*UIOO(j).J);

    for i = 1:M
        UIOO(j).O(i).W = double(UIOO(j).O(i).W);

        UIOO(j).O(i).K = X\UIOO(j).O(i).W;
        UIOO(j).O(i).G = UIOO(j).T1*sys(i).Bd;
        UIOO(j).O(i).Tg = UIOO(j).T1*sys(i).deltad;
        UIOO(j).O(i).N = UIOO(j).T1*sys(i).Ad - UIOO(j).O(i).K*UIOO(j).J;
        UIOO(j).O(i).L = UIOO(j).O(i).K - UIOO(j).O(i).N*UIOO(j).E;
    end
end