function mu = membership(x, Theta_1s_min, Theta_1s_mid, Theta_1s_max, Theta_2s_min, Theta_2s_mid, Theta_2s_max)
    % Estados
    % x(1) => Theta_1s
    % x(2) => Theta_2s
    % x(3) => Theta_p

    Theta_1s = x(1);
    Theta_2s = x(2);
    
    % Saturación de estados
    if Theta_1s < Theta_1s_min
        Theta_1s = Theta_1s_min;
    elseif Theta_1s > Theta_1s_max
        Theta_1s = Theta_1s_max;
    end
	if Theta_2s < Theta_2s_min
        Theta_2s = Theta_2s_min;
    elseif Theta_2s > Theta_2s_max
        Theta_2s = Theta_2s_max;
    end
    
    % Membresía estado Theta_1s
    if Theta_1s < Theta_1s_mid
        M1_1 = (Theta_1s_mid-Theta_1s)/(Theta_1s_mid-Theta_1s_min);
        M1_2 = (Theta_1s-Theta_1s_min)/(Theta_1s_mid-Theta_1s_min);
        M1_3 = 0;
    else
        M1_1 = 0;
        M1_2 = (Theta_1s_max-Theta_1s)/(Theta_1s_max-Theta_1s_mid);
        M1_3 = (Theta_1s-Theta_1s_mid)/(Theta_1s_max-Theta_1s_mid);
    end
    
        % Membresía estado Theta_2s
    if Theta_2s < Theta_2s_mid
        M2_1 = (Theta_2s_mid-Theta_2s)/(Theta_2s_mid-Theta_2s_min);
        M2_2 = (Theta_2s-Theta_2s_min)/(Theta_2s_mid-Theta_2s_min);
        M2_3 = 0;
    else
        M2_1 = 0;
        M2_2 = (Theta_2s_max-Theta_2s)/(Theta_2s_max-Theta_2s_mid);
        M2_3 = (Theta_2s-Theta_2s_mid)/(Theta_2s_max-Theta_2s_mid);
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