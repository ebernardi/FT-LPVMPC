function mu = membership(x, V_min, V_mid, V_max, Tr_min, Tr_mid, Tr_max)
    % Estados
    % x(1) => V
    % x(2) => CA

    V = x(1);
    T = x(3);
    % Saturación de estados
    if V < V_min
        V = V_min;
    elseif V > V_max
        V = V_max;
    end
	if T < Tr_min
        T = Tr_min;
    elseif T > Tr_max
        T = Tr_max;
	end
    
    % Membresía estado V
    if V < V_mid
        M1_1 = (V_mid-V)/(V_mid-V_min);
        M1_2 = (V-V_min)/(V_mid-V_min);
        M1_3 = 0;
    else
        M1_1 = 0;
        M1_2 = (V_max-V)/(V_max-V_mid);
        M1_3 = (V-V_mid)/(V_max-V_mid);
    end
    
    % Membresía estado CA
    if T < Tr_mid
        M2_1 = (Tr_mid-T)/(Tr_mid-Tr_min);
        M2_2 = (T-Tr_min)/(Tr_mid-Tr_min);
        M2_3 = 0;
    else
        M2_1 = 0;
        M2_2 = (Tr_max-T)/(Tr_max-Tr_mid);
        M2_3 = (T-Tr_mid)/(Tr_max-Tr_mid);        
    end

    mu1 = M1_1*M2_1;
    mu2 = M1_1*M2_2;
    mu3 = M1_1*M2_3;
    mu4 = M1_2*M2_1;
    mu5 = M1_2*M2_2;
    mu6 = M1_2*M2_3;
    mu7 = M1_3*M2_1;
    mu8 = M1_3*M2_2;
    mu9 = M1_3*M2_3;
    
    mu = [mu1 mu2 mu3 mu4 mu5 mu6 mu7 mu8 mu9];
end
