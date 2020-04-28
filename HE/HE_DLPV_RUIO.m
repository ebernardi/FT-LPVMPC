% Type of observer(N° observer). Sub-observer(N° sub-observer). Matrix 
RUIO = struct;

% Dimension of system matrices
n = 2; p = 3;

for j = 1:N
    
    if j == 1
        RUIO(j).N = [0 0; 1 0; 0 1];
        RUIO(j).Q = [0 0; 0 1; 1 0];
        RUIO(j).alpha = -0.1;
    else
        RUIO(j).N = [0 0.1; 0 0; 1 0];
        RUIO(j).Q = [1 0; 0 0; 0 1];
        RUIO(j).alpha = -0.1;
    end
    
    % LMIs
    setlmis([]);
    X = lmivar(1, [n 1]);

    % LMI #1: X > 0
    lmiterm([-1 1 1 X], 1, 1);                  % LMI #1: X

    for i = 1:M

        % Transformation matrix T for each model
        RUIO(j).O(i).D = sys(i).Bd(:, j);
        RUIO(j).O(i).B = sys(i).Bd;
        RUIO(j).O(i).N = RUIO(j).N;
        RUIO(j).O(i).T = [RUIO(j).N RUIO(j).O(i).D];
        RUIO(j).O(i).A_bar = RUIO(j).O(i).T\sys(i).Ad*RUIO(j).O(i).T;
        RUIO(j).O(i).B_bar = RUIO(j).O(i).T\RUIO(j).O(i).B;
        RUIO(j).O(i).B_bar_1 = RUIO(j).O(i).B_bar(1:2, :);
        RUIO(j).O(i).B_bar_2 = RUIO(j).O(i).B_bar(3, :);
        RUIO(j).O(i).D_bar = RUIO(j).O(i).T\RUIO(j).O(i).D;
        RUIO(j).O(i).delta_bar = RUIO(j).O(i).T\sys(i).deltad;
        RUIO(j).O(i).delta_bar_1 = RUIO(j).O(i).delta_bar(1:2, :);
        RUIO(j).O(i).delta_bar_2 = RUIO(j).O(i).delta_bar(3, :);

        % Split state vector
        RUIO(j).O(i).A_bar_11 = RUIO(j).O(i).A_bar(1:2, 1:2);
        RUIO(j).O(i).A_bar_12 = RUIO(j).O(i).A_bar(1:2, 3);
        RUIO(j).O(i).A_bar_21 = RUIO(j).O(i).A_bar(3, 1:2);
        RUIO(j).O(i).A_bar_22 = RUIO(j).O(i).A_bar(3, 3);

        % Matrix U
        RUIO(j).O(i).rank_CD = rank(sys(i).C*RUIO(j).O(i).D);
        RUIO(j).O(i).Q = RUIO(j).Q;
        RUIO(j).O(i).U = [sys(i).C*RUIO(j).O(i).D RUIO(j).O(i).Q];
        RUIO(j).O(i).inv_U = inv(RUIO(j).O(i).U);
        RUIO(j).O(i).U_1 = RUIO(j).O(i).inv_U(1, :);
        RUIO(j).O(i).U_2 = RUIO(j).O(i).inv_U(2:3, :);

        % More matrices and check observability
        RUIO(j).O(i).A_tilde_1 = RUIO(j).O(i).A_bar_11 - RUIO(j).O(i).A_bar_12*RUIO(j).O(i).U_1*sys(i).C*RUIO(j).O(i).N;
        RUIO(j).O(i).C_tilde_1 = sys(i).C*RUIO(j).O(i).N;
        RUIO(j).O(i).E_1 = RUIO(j).O(i).A_bar_12*RUIO(j).O(i).U_1;

        RUIO(j).O(i).rank_C_tilde_1 = rank(RUIO(j).O(i).C_tilde_1);
        RUIO(j).O(i).O_M = [RUIO(j).O(i).C_tilde_1' (RUIO(j).O(i).C_tilde_1*RUIO(j).O(i).A_tilde_1)' (RUIO(j).O(i).C_tilde_1*RUIO(j).O(i).A_tilde_1*RUIO(j).O(i).A_tilde_1)']';
        RUIO(j).O(i).rank_Obs_M = rank(RUIO(j).O(i).O_M);

        % LMI
        RUIO(j).O(i).W = lmivar(2, [n p]);

        % LMI #i+1: M < 0
        lmiterm([i+1 1 1 X], 2*RUIO(j).alpha, 1);  % LMI #i+1: 2*alpha*X; −left hand side
        lmiterm([i+1 2 1 RUIO(j).O(i).W], -1, RUIO(j).O(i).C_tilde_1, 's'); % LMI #i+1: -W*C_tilde_1; −left hand side
        lmiterm([i+1 2 1 X], 1, RUIO(j).O(i).A_tilde_1, 's'); % LMI #i+1: X*A_tilde_1; −left hand side
        lmiterm([i+1 2 2 X], 2*RUIO(j).alpha, 1);  % LMI #i+1: 2*alpha*X; −left hand side
    end

    LMIs = getlmis;

    [~, xfeas] = feasp(LMIs);

    X = dec2mat(LMIs, xfeas, X);

    for i = 1:M
        RUIO(j).O(i).W = dec2mat(LMIs, xfeas, RUIO(j).O(i).W);
        RUIO(j).O(i).L = X\RUIO(j).O(i).W;
        RUIO(j).O(i).K = RUIO(j).O(i).A_tilde_1 - RUIO(j).O(i).L*RUIO(j).O(i).C_tilde_1;
        RUIO(j).O(i).L_ast = RUIO(j).O(i).L + RUIO(j).O(i).E_1;
    end
end