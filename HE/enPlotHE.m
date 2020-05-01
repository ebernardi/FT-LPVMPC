% clc; clear; close all;

% % Load data
% load runHE

% When generates flat figures
set(0, 'DefaultFigureRenderer', 'painters');

vecrojo = [0.7; 0; 0]; vecverde = [0; 0.8; 0]; vecazul = [0; 0; 0.6]; negro = [.1; .1; .1]; gris = [.5; .7; .5];

azul = [0 0.4470 0.7410];
naranja = [0.8500 0.3250 0.0980];
amarillo = [0.9290 0.6940 0.1250];
violeta = [0.4940 0.1840 0.5560];
verde = [0.4660 0.6740 0.1880];
celeste = [0.3010 0.7450 0.9330];
bordo = [0.6350 0.0780 0.1840];

orange_red = [255 69 0]/255;
forest_green = [34 139 34]/255;
royal_blue = [65 105 225]/255;
dark_blue = [0 0 139]/255;
gold = [255 215 0]/255;
chocolate = [210 105 30]/255;
arrow = [212 55 144]/255;

disp('Plotting')
for FT = 1:2    % 1 - FT is off; 2 -  FT is on
    
	if FT == FTC_ON
        disp('Fault tolerant = ON')
	else
        disp('Fault tolerant = OFF')
    end
    
	error = abs(FTCS(FT).Y(2, :) - FTCS(FT).Xsp(2, :));
    IAE = trapz(t, abs(error));
    ISE = trapz(t, error.^2);
    ITAE = trapz(t, t.*abs(error));
    RMS = sqrt(mean(error.^2));
    msg = ['IAE = ', num2str(IAE)];
    disp(msg)
    msg = ['ITAE = ', num2str(ITAE)];
    disp(msg)
    msg = ['ISE = ', num2str(ISE)];
    disp(msg)
    msg = ['RMS = ', num2str(RMS)];
    disp(msg)    
    
    time_avg = mean(FTCS(FT).elapsed_time) ;
    msg = ['Mean time = ', num2str(time_avg)];
    disp(msg)
    time_avg = max(FTCS(FT).elapsed_time) ;
    msg = ['Max time = ', num2str(time_avg)];
    disp(msg)
    time_avg = min(FTCS(FT).elapsed_time) ;
    msg = ['Min time = ', num2str(time_avg)];
    disp(msg)    

    %% Outputs
    fig = figure('Name', 'Outputs');
    subplot(311)
    plot(t, FTCS(FT).Xsp(1, :), 'r-.', 'LineWidth', 1.5);
    hold on
    plot(t, FTCS(FT).Y(1, :), 'g--', t, FTCS(FT).Yfail(1, :), ':k', 'LineWidth', 1.5); hold off
    xlabel('Time [min]'); ylabel('\theta_{1_s} [K]'); grid on
    axis([0 inf 494 502])
    leg = legend('Setpoint', 'Actual', 'Measured', 'Location', 'SouthEast');
    set(leg, 'Position', [0.748 0.764 0.148 0.109], 'FontSize', 8);
    leg.ItemTokenSize = [20, 15];
    subplot(312)
    plot(t, FTCS(FT).Xsp(2, :), 'r-.', 'LineWidth', 1.5);
    hold on
    plot(t, FTCS(FT).Y(2, :), 'g--', t, FTCS(FT).Yfail(2, :), ':k', 'LineWidth', 1.5); hold off
    xlabel('Time [min]'); ylabel('\theta_{2_s} [K]'); grid on
    axis([0 inf 675 715])
    subplot(313)
    plot(t, FTCS(FT).Y(3, :), 'g--', 'LineWidth', 1.5);
    xlabel('Time [min]'); ylabel('\theta_p [K]'); grid on
    axis([0 inf 556 568])
    
%     print -dsvg figs/FTCS_HE_outputs.svg
   
    %% State space
    fig = figure('Name', 'State space');
    plot(FTCS(FT).X(3, :), FTCS(FT).X(2, :), '-', 'Color', azul, 'LineWidth', 1.5);
    hold on
    plot(FTCS(FT).X(3, 1), FTCS(FT).X(2, 1), 'o', 'Color', verde, 'LineWidth', 1.5);
    plot(FTCS(FT).Xsp(3, :), FTCS(FT).Xsp(2, :), '*', 'Color', bordo, 'LineWidth', 1.5);
    hold off

    %% Membership
    fig = figure('Name', 'Membership');
    hold on
    for l=1:M
        plot(t, FTCS(FT).mu_mhe(l, :));
        legendInfo{l} = ['\mu' num2str(l)];
        plot(t, FTCS(FT).mu_fuzzy(l, :), '--');
        legendInfo{l+1} = ['\mu' num2str(l+1)];
    end
    ylabel('\mu'), xlabel('Time [min]')
    legend(legendInfo);
    xlim([0 Time])
    hold off
    
    %% Input
    fig = figure('Name', 'Inputs');
    subplot(2, 1, 1)
    stairs(t, FTCS(FT).U(1, :), 'Color', orange_red, 'LineWidth', 1.5);
    xlabel('Time [s]'); ylabel('q_1 [l/m]'); grid on
    xlim([0 Time])
    subplot(2, 1, 2)
    stairs(t, FTCS(FT).U(2, :), 'Color', orange_red, 'LineWidth', 1.5);
    xlabel('Time [min]'); ylabel('q_2 [l/m]'); grid on
    xlim([0 Time])

%     print -dsvg figs/FTCS_HE_input.svg    

    %% RUIO error detection
    fig = figure('Name', 'RUIO Error');
    subplot(211)
    plot(t, FTCS(FT).RUIO(1).error, 'b', 'LineWidth', 1.5)
    hold on; grid on
    plot(t, threshold(1, :),  'r--', 'LineWidth', 1.5)
    hold off
    axis([0 inf 0 2.5])
    xlabel('Time [min]'); ylabel('|e_q|');

    subplot(212)
    plot(t, FTCS(FT).RUIO(2).error, 'b', 'LineWidth', 1.5)
    hold on; grid on
    plot(t, threshold(2, :),  'r--', 'LineWidth', 1.5)
    hold off
    axis([0 inf 0 1.2])
    xlabel('Time [min]'); ylabel('|e_q|');
%     legend('MPC', 'FTMPC', 'Threshold', 'Location', 'NorthWest');

    %% UIOO error detection
    fig = figure('Name', 'UIOO Error');
    subplot(211)
    plot(t, FTCS(FT).UIOO(1).error, 'b', 'LineWidth', 1.5)
    hold on; grid on
    plot(t, threshold(3, :),  'r--', 'LineWidth', 1.5)
    hold off
    axis([0 inf 0 2])
    xlabel('Time [min]'); ylabel('|e_x|');

    subplot(212)
    plot(t, FTCS(FT).UIOO(2).error, 'b', 'LineWidth', 1.5)
    hold on; grid on
    plot(t, threshold(4, :),  'r--', 'LineWidth', 1.5)
    hold off
    axis([0 inf 0 4])
    xlabel('Time [min]'); ylabel('|e_x|');
%     legend('MPC', 'FTMPC', 'Threshold', 'Location', 'NorthWest');

    %% Actuator fault estimation
    fig = figure('Name', 'Actuator fault');
    subplot(211)
    stairs(t, FTCS(FT).RUIO(1).Fact, 'b', 'LineWidth', 1.5)
    hold on; grid on
    stairs(t, FTCS(FT).Ufails(1, :), 'm--', 'LineWidth', 1.5)
    hold off
    axis([0 inf 0 6])
    xlabel('Time [min]'); ylabel('Q_1 [l/min]');
%     legend('MPC', 'Actual fault', 'FTMPC', 'Location', 'NorthEast');

    subplot(212)
    stairs(t, FTCS(FT).RUIO(2).Fact, 'b', 'LineWidth', 1.5)
    hold on; grid on
    stairs(t, FTCS(FT).Ufails(2, :), 'm--', 'LineWidth', 1.5)
    hold off
    axis([0 inf -0.5 0])
    xlabel('Time [min]'); ylabel('Q_2 [l/min]');

    %% Sensor fault estimation
    fig = figure('Name', 'Sensor fault');
    subplot(211)
    stairs(t, FTCS(FT).UIOO(1).Fsen, 'b', 'LineWidth', 1.5)
    hold on; grid on
    stairs(t, FTCS(FT).Yfails(1, :), 'm--', 'LineWidth', 1.5)
    hold off
    axis([0 inf 0 4])
    xlabel('Time [min]'); ylabel('\theta_1 [K]');
%     legend('MPC', 'Actual fault', 'FTMPC', 'Location', 'NorthEast');

    subplot(212)
    stairs(t, FTCS(FT).UIOO(2).Fsen, 'b', 'LineWidth', 1.5)
    hold on; grid on
    stairs(t, FTCS(FT).Yfails(2, :), 'm--', 'LineWidth', 1.5)
    hold off
    axis([0 inf -4 0])
    xlabel('Time [min]'); ylabel('\theta_2 [K]');
    
    %% Objective function
    fig = figure('Name', 'Objective function');
    plot(t, FTCS(FT).Obj(:))
    xlim([0 Time])
    xlabel('Time [min]'); ylabel('Cost');
    
end