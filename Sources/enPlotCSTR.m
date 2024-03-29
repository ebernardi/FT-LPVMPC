clc; clear; close all;

% % Load data
load FTCS

% When generates flat figures
set(0, 'DefaultFigureRenderer', 'painters');

% Indices
disp('MHE-MPC - Volume')
error = abs(FTCS(FTC_OFF).Y(1, :) - FTCS(FTC_OFF).Xsp(1));
IAE = trapz(t, abs(error));
ISE = trapz(t, error.^2);
ITAE = trapz(t, t.*abs(error));
ITSE = trapz(t, t.*error.^2);
msg = ['IAE = ', num2str(IAE)];
disp(msg)
msg = ['ITAE = ', num2str(ITAE)];
disp(msg)
msg = ['ISE = ', num2str(ISE)];
disp(msg)
msg = ['ITSE = ', num2str(ITSE)];
disp(msg)
disp(' ')

disp('MHE-MPC - Ca')
error = abs(FTCS(FTC_OFF).Y(2, :) - FTCS(FTC_OFF).Xsp(2));
IAE = trapz(t, abs(error));
ISE = trapz(t, error.^2);
ITAE = trapz(t, t.*abs(error));
ITSE = trapz(t, t.*error.^2);
msg = ['IAE = ', num2str(IAE)];
disp(msg)
msg = ['ITAE = ', num2str(ITAE)];
disp(msg)
msg = ['ISE = ', num2str(ISE)];
disp(msg)
msg = ['ITSE = ', num2str(ITSE)];
disp(msg)
disp(' ')

disp('MHE-MPC - Temperature')
error = abs(FTCS(FTC_OFF).Y(3, :) - FTCS(FTC_OFF).Xsp(3));
IAE = trapz(t, abs(error));
ISE = trapz(t, error.^2);
ITAE = trapz(t, t.*abs(error));
ITSE = trapz(t, t.*error.^2);
msg = ['IAE = ', num2str(IAE)];
disp(msg)
msg = ['ITAE = ', num2str(ITAE)];
disp(msg)
msg = ['ISE = ', num2str(ISE)];
disp(msg)
msg = ['ITSE = ', num2str(ITSE)];
disp(msg)
disp(' ')

FTCS(FTC_OFF).input_shifted = [0; FTCS(FTC_OFF).U(1, 1:end-1)'];
msg = ['MHE-MPC - TV = ', num2str(sum(abs(FTCS(FTC_OFF).input_shifted-FTCS(FTC_OFF).U'))/60)];
disp(msg)
disp(' ')

disp('AFTCS - Volume')
error = abs(FTCS(FTC_ON).Y(1, :) - FTCS(FTC_ON).Xsp(1));
IAE = trapz(t, abs(error));
ISE = trapz(t, error.^2);
ITAE = trapz(t, t.*abs(error));
ITSE = trapz(t, t.*error.^2);
msg = ['IAE = ', num2str(IAE)];
disp(msg)
msg = ['ITAE = ', num2str(ITAE)];
disp(msg)
msg = ['ISE = ', num2str(ISE)];
disp(msg)
msg = ['ITSE = ', num2str(ITSE)];
disp(msg)
disp(' ')

disp('AFTCS - Ca')
error = abs(FTCS(FTC_ON).Y(2, :) - FTCS(FTC_ON).Xsp(2));
IAE = trapz(t, abs(error));
ISE = trapz(t, error.^2);
ITAE = trapz(t, t.*abs(error));
ITSE = trapz(t, t.*error.^2);
msg = ['IAE = ', num2str(IAE)];
disp(msg)
msg = ['ITAE = ', num2str(ITAE)];
disp(msg)
msg = ['ISE = ', num2str(ISE)];
disp(msg)
msg = ['ITSE = ', num2str(ITSE)];
disp(msg)
disp(' ')

disp('AFTCS - Temperature')
error = abs(FTCS(FTC_ON).Y(3, :) - FTCS(FTC_ON).Xsp(3));
IAE = trapz(t, abs(error));
ISE = trapz(t, error.^2);
ITAE = trapz(t, t.*abs(error));
ITSE = trapz(t, t.*error.^2);
msg = ['IAE = ', num2str(IAE)];
disp(msg)
msg = ['ITAE = ', num2str(ITAE)];
disp(msg)
msg = ['ISE = ', num2str(ISE)];
disp(msg)
msg = ['ITSE = ', num2str(ITSE)];
disp(msg)
disp(' ')

FTCS(FTC_ON).input_shifted = [0; FTCS(FTC_ON).U(1, 1:end-1)'];
msg = ['AFTCS - TV = ', num2str(sum(abs(FTCS(FTC_ON).input_shifted-FTCS(FTC_ON).U'))/60)];
disp(msg)

% error = abs(X(2, 1:5000) - Xsp(2));
% IAEtrack = trapz(t(1:5000), abs(error));
% error = abs(X(2, 5001:end) - Xsp(2));
% IAEpert = trapz(t(5001:end), abs(error));
% msg = ['IAE tracking = ', num2str(IAEtrack)];
% disp(msg)
% msg = ['IAE disturbance = ', num2str(IAEpert)];
% disp(msg)
% disp(' ')

vecrojo = [0.7; 0; 0]; vecverde = [0; 0.8; 0]; vecazul = [0; 0; 0.6]; negro = [.1; .1; .1]; gris = [.5; .7; .5];

red = [0.7; 0; 0];
black = [.1; .1; .1];
gray = [.5; .7; .5];
blue = [0 0.4470 0.7410];
orange = [0.8500 0.3250 0.0980];
yellow = [0.9290 0.6940 0.1250];
violet = [0.4940 0.1840 0.5560];
green = [0.4660 0.6740 0.1880];
lightblue = [0.3010 0.7450 0.9330];
purple = [0.6350 0.0780 0.1840];

orange_red = [255 69 0]/255;
forest_green = [34 139 34]/255;
royal_blue = [65 105 225]/255;
dark_blue = [0 0 139]/255;
gold = [255 215 0]/255;
chocolate = [210 105 30]/255;
arrow = [212 55 144]/255;

disp('Plotting')
% for FT = FTC_OFF:FTC_ON    % 1 - FT is off; 2 -  FT is on
%     
% 	if FT == FTC_ON
%         disp('Fault tolerant = ON')
% 	else
%         disp('Fault tolerant = OFF')
% 	end
%     
% 	error = abs(FTCS(FT).Y(2, :) - FTCS(FT).Xsp(2, :));
%     IAE = trapz(t, abs(error));
%     ISE = trapz(t, error.^2);
%     ITAE = trapz(t, t.*abs(error));
%     RMS = sqrt(mean(error.^2));
%     msg = ['IAE = ', num2str(IAE)];
%     disp(msg)
%     msg = ['ITAE = ', num2str(ITAE)];
%     disp(msg)
%     msg = ['ISE = ', num2str(ISE)];
%     disp(msg)
%     msg = ['RMS = ', num2str(RMS)];
%     disp(msg)    
%     
%     time_avg = mean(FTCS(FT).elapsed_time) ;
%     msg = ['Mean time = ', num2str(time_avg)];
%     disp(msg)
%     time_avg = max(FTCS(FT).elapsed_time) ;
%     msg = ['Max time = ', num2str(time_avg)];
%     disp(msg)
%     time_avg = min(FTCS(FT).elapsed_time) ;
%     msg = ['Min time = ', num2str(time_avg)];
%     disp(msg)    
% 
%     %% Outputs
%     fig = figure('Name', 'Outputs');
%     subplot(211)
%     plot(t, FTCS(FT).Xsp(1, :), 'r-.', 'LineWidth', 1.5);
%     hold on
%     plot(t, FTCS(FT).Y(1, :), 'g--', t, FTCS(FT).Yfail(1, :), ':k', 'LineWidth', 1.5); hold off
%     xlabel('Time [min]'); ylabel('V [l]');
%     grid on
%     axis([0 inf 90 110])
%     leg = legend('Setpoint', 'Actual', 'Measured', 'Location', 'NorthWest');
%     leg.ItemTokenSize = [20, 15];
% %     subplot(312)
% %     plot(t, FTCS(FT).Y(2, :), 'g--', 'LineWidth', 1.5);
% %     xlabel('Time [min]'); ylabel('C_A [mol/l]'); grid on
% %     axis([0 inf 0.04 0.12])
%     subplot(212)
%     plot(t, FTCS(FT).Xsp(3, :), 'r-.', 'LineWidth', 1.5);
%     hold on
%     plot(t, FTCS(FT).Y(3, :), 'g--', t, FTCS(FT).Yfail(3, :), ':k', 'LineWidth', 1.5); hold off
%     xlabel('Time [min]'); ylabel('T [K]'); grid on
%     axis([0 inf 438 452])
%     
% %     print -dsvg ../Figs/FTCS_CSTR_outputs.svg
% 
%     %% State space
%     fig = figure('Name', 'State space');
%     plot(FTCS(FT).X(3, :), FTCS(FT).X(2, :), '-', 'Color', azul, 'LineWidth', 1.5);
%     hold on
%     plot(FTCS(FT).X(3, 1), FTCS(FT).X(2, 1), 'o', 'Color', verde, 'LineWidth', 1.5);
%     plot(FTCS(FT).Xsp(3, :), FTCS(FT).Xsp(2, :), '*', 'Color', bordo, 'LineWidth', 1.5);
%     hold off
%     
% %     %% State space
% %     fig = figure('Name', 'State space');
% %     plot3(FTCS(FT).X(1, :), FTCS(FT).X(2, :), FTCS(FT).X(3, :), '-', 'Color', azul, 'LineWidth', 1.5);
% %     hold on
% %     plot3(FTCS(FT).X(1, 1), FTCS(FT).X(2, 1), FTCS(FT).X(3, 1), 'o', 'Color', verde, 'LineWidth', 1.5);
% %     plot3(FTCS(FT).Xsp(1, :), FTCS(FT).Xsp(2, :), FTCS(FT).Xsp(3, :), '*', 'Color', bordo, 'LineWidth', 1.5);
% %     ellipse(Wbmi, xsp, 20, 'red', '-')
% %     plot(Xx)
% %     hold off
%     
%     %% Membership
%     fig = figure('Name', 'Membership');
%     hold on
%     for l=1:M
%         plot(t, FTCS(FT).mu_mhe(l, :));
%         legendInfo{l} = ['\mu' num2str(l)];
%         plot(t, FTCS(FT).mu_fuzzy(l, :), '--');
%         legendInfo{l+1} = ['\mu' num2str(l+1)];
%     end
%     ylabel('\mu'), xlabel('Time [min]')
%     legend(legendInfo);
%     xlim([0 Time])
%     hold off
% 
%     %% Input
%     fig = figure('Name', 'Inputs');
%     subplot(2, 1, 1)
%     stairs(t, FTCS(FT).U(1, :), 'Color', orange_red, 'LineWidth', 1.5);
%     xlabel('Time [min]'); ylabel('Qs [l/m]'); grid on
%     xlim([0 Time])
%     subplot(2, 1, 2)
%     stairs(t, FTCS(FT).U(2, :), 'Color', orange_red, 'LineWidth', 1.5);
%     xlabel('Time [min]'); ylabel('Qc [l/m]'); grid on
%     xlim([0 Time])
% 
% %     print -dsvg ../Figs/FTCS_CSTR_input.svg
% 
%     %% RUIO error detection
%     fig = figure('Name', 'RUIO Error');
%     subplot(211)
%     plot(t, FTCS(FT).RUIO(1).error, 'b', 'LineWidth', 1.5)
%     hold on; grid on
%     plot(t, threshold(1, :),  'r--', 'LineWidth', 1.5)
%     hold off
%     axis([0 inf 0 0.6])
%     xlabel('Time [min]'); ylabel('|e_q|');
% 
%     subplot(212)
%     plot(t, FTCS(FT).RUIO(2).error, 'b', 'LineWidth', 1.5)
%     hold on; grid on
%     plot(t, threshold(2, :),  'r--', 'LineWidth', 1.5)
%     hold off
%     axis([0 inf 0 0.3])
%     xlabel('Time [min]'); ylabel('|e_q|');
% %     legend('MPC', 'FTMPC', 'Threshold', 'Location', 'NorthWest');
% 
%     %% UIOO error detection
%     fig = figure('Name', 'UIOO Error');
%     subplot(211)
%     plot(t, FTCS(FT).UIOO(1).error, 'b', 'LineWidth', 1.5)
%     hold on; grid on
%     plot(t, threshold(3, :),  'r--', 'LineWidth', 1.5)
%     hold off
%     axis([0 inf 0 0.3])
%     xlabel('Time [min]'); ylabel('|e_x|');
% 
%     subplot(212)
%     plot(t, FTCS(FT).UIOO(2).error, 'b', 'LineWidth', 1.5)
%     hold on; grid on
%     plot(t, threshold(4, :),  'r--', 'LineWidth', 1.5)
%     hold off
%     axis([0 inf 0 0.3])
%     xlabel('Time [min]'); ylabel('|e_x|');
% %     legend('MPC', 'FTMPC', 'Threshold', 'Location', 'NorthWest');
% 
%     %% Actuator fault estimation
%     fig = figure('Name', 'Actuator fault');
%     subplot(211)
%     stairs(t, FTCS(FT).RUIO(1).Fact, 'b', 'LineWidth', 1.5)
%     hold on; grid on
%     stairs(t, FTCS(FT).Ufails(1, :), 'm--', 'LineWidth', 1.5)
%     hold off
% %     axis([0 inf -6 0])
%     xlabel('Time [min]'); ylabel('Qs [l/min]');
% %     legend('MPC', 'Actual fault', 'FTMPC', 'Location', 'NorthEast');
% 
%     subplot(212)
%     stairs(t, FTCS(FT).RUIO(2).Fact, 'b', 'LineWidth', 1.5)
%     hold on; grid on
%     stairs(t, FTCS(FT).Ufails(2, :), 'm--', 'LineWidth', 1.5)
%     hold off
% %     axis([0 inf -0.5 0])
%     xlabel('Time [min]'); ylabel('Qc [l/min]');
%     
%     %% Objective function
%     fig = figure('Name', 'Objective function');
%     plot(t, FTCS(FT).Obj(:))
%     xlim([0 Time])
%     xlabel('Time [min]'); ylabel('Cost');
%     
% end

%% Final figures
%% Output with and without fault-tolerant controller
fig = figure('Name', 'Outputs');
subplot(211)
plot(t, FTCS(FTC_ON).Xsp(1, :), 'r-.', 'LineWidth', 1.5);
hold on
plot(t, FTCS(FTC_OFF).Y(1, :), 'b', t, FTCS(FTC_ON).Y(1, :), 'k--', 'LineWidth', 1.5);
xlabel('Time [min]'); ylabel('V [l]'); grid on; hold off;
axis([0 Time 89 111])
leg = legend('Setpoint', 'MPC', 'FT-MPC', 'Location', 'NorthEast');
leg.ItemTokenSize = [20, 15];
% subplot(312)
% plot(t, FTCS(FTC_ON).Xsp(2, :), 'r-.', 'LineWidth', 1.5);
% hold on
% plot(t, FTCS(FTC_OFF).Y(2, :), 'b', t, FTCS(FTC_ON).Y(2, :), 'k--', 'LineWidth', 1.5);
% xlabel('Time [min]'); ylabel('C_A [mol/l]'); grid on; hold off;
% axis([0 Time 0.06 0.1])
subplot(212)
plot(t, FTCS(FTC_ON).Xsp(3, :), 'r-.', 'LineWidth', 1.5);
hold on
plot(t, FTCS(FTC_OFF).Y(3, :), 'b', t, FTCS(FTC_ON).Y(3, :), 'k--', 'LineWidth', 1.5); hold off
xlabel('Time [min]'); ylabel('T [K]'); grid on; hold off;
axis([0 Time 440 448])

% Create textarrow
annotation(fig, 'textarrow',[0.4089 0.3929], [0.738 0.707], ...
    'String', {'Q_s fault', 'income'}, 'LineWidth', 1, 'HorizontalAlignment', 'center', ...
    'HeadWidth', 6, 'HeadLength', 6, 'FontSize', 8);
annotation(fig, 'textarrow', [0.611 0.598], [0.3849 0.338], ...
    'String', {'FT-MPC recovery'}, 'LineWidth', 1, 'HorizontalAlignment', 'center', ...
    'HeadWidth', 6, 'HeadLength', 6, 'FontSize', 8);
annotation(fig, 'textarrow', [0.402 0.4315], [0.154 0.173], ...
    'String', {'MPC divergence'}, 'LineWidth', 1, 'HorizontalAlignment','center', ...
    'HeadWidth', 6, 'HeadLength', 6, 'FontSize', 8);
annotation(fig, 'textarrow', [0.232 0.22], [0.28 0.3142], ...
    'String', {'Temp. sensor', 'fault income'}, 'LineWidth', 1, 'HorizontalAlignment', 'center',...
    'HeadWidth', 6, 'HeadLength', 6, 'FontSize', 8);
annotation(fig, 'textarrow', [0.7525 0.741], [0.263 0.309], ...
    'String', {'Q_c fault', 'income'}, 'LineWidth', 1, 'HorizontalAlignment', 'center', ...
    'HeadWidth', 6, 'HeadLength', 6, 'FontSize', 8);
annotation(fig, 'textarrow', [0.591 0.5698], [0.809 0.746], ...
    'String', {'Volume sensor','fault income'}, 'LineWidth', 1, 'HorizontalAlignment', 'center', ...
    'HeadWidth', 6, 'HeadLength', 6, 'FontSize', 8);
    
print -dsvg ../Figs/FTCS_CSTR_outputs.svg

%% Actuator fault estimation
fig = figure('Name', 'Actuator fault estimation');
subplot(211)
stairs(t, FTCS(FTC_ON).Ufails(1, :), 'Color', orange, 'LineWidth', 1.5)
hold on;
stairs(t, FTCS(FTC_ON).RUIO(1).Fact, '--', 'Color', blue, 'LineWidth', 1.5)
axis([0 Time 0 5])
xlabel('Time [min]'); ylabel('Qs [l/min]'); grid on; hold off;
leg = legend('Fault', 'Estimation', 'Location', 'NorthWest');
leg.ItemTokenSize = [20, 15];

subplot(212)
stairs(t, FTCS(FTC_ON).Ufails(2, :), 'Color', orange, 'LineWidth', 1.5)
hold on;
stairs(t, FTCS(FTC_ON).RUIO(2).Fact, '--', 'Color', blue, 'LineWidth', 1.5)
axis([0 Time -6 0])
xlabel('Time [min]'); ylabel('Qc [l/min]'); grid on; hold off;

% Create axes
ax = axes('Parent', fig, 'Position', [0.598 0.666 0.267 0.211], 'FontSize', 8);
hold(ax, 'on');
stairs(t, FTCS(FTC_ON).Ufails(1, :), 'Color', orange, 'LineWidth', 1.5)
stairs(t, FTCS(FTC_ON).RUIO(1).Fact, '--', 'Color', blue, 'LineWidth', 1.5)
xlim(ax, [28 42]); ylim(ax, [3.99 4.01]);
box(ax, 'on'); grid(ax, 'on');

% Create axes
ax = axes('Parent', fig, 'Position', [0.217 0.197 0.309 0.212], 'FontSize', 8);
hold(ax, 'on');
stairs(t, FTCS(FTC_ON).Ufails(2, :), 'Color', orange, 'LineWidth', 1.5)
stairs(t, FTCS(FTC_ON).RUIO(2).Fact, '--', 'Color', blue, 'LineWidth', 1.5)
xlim(ax, [69.8 71.5]); ylim(ax, [-5 0]);
box(ax, 'on'); grid(ax, 'on');

% Create textarrow
annotation(fig, 'textarrow', [0.385 0.352], [0.336 0.335], ...
    'String', {'Threshold', 'effect'}, 'LineWidth', 1, ...
    'HorizontalAlignment', 'center', 'HeadWidth', 6, ...
    'HeadLength', 6, 'FontSize', 8);

print -dsvg ../Figs/FTCS_CSTR_actuator_estimation.svg

%% Sensor fault estimation
fig = figure('Name', 'Sensor fault estimation');
subplot(211)
stairs(t, FTCS(FTC_ON).Yfails(1, :), 'Color', orange, 'LineWidth', 1.5)
hold on; grid on
stairs(t, FTCS(FTC_ON).UIOO(1).Fsen, '--', 'Color', blue, 'LineWidth', 1.5)
hold off
axis([0 Time 0 2.5])
xlabel('Time [min]'); ylabel('V [l]');
leg = legend('Fault', 'Estimation', 'Location', 'NorthEast');
leg.ItemTokenSize = [20, 15];

subplot(212)
stairs(t, FTCS(FTC_ON).Yfails(3, :), 'Color', orange, 'LineWidth', 1.5)
hold on; grid on
stairs(t, FTCS(FTC_ON).UIOO(2).Fsen, '--', 'Color', blue, 'LineWidth', 1.5)
hold off
axis([0 Time -2.5 0])
xlabel('Time [min]'); ylabel('T [K]');

% Create axes
ax = axes('Parent', fig, 'Position', [0.47 0.191 0.344 0.204], 'FontSize', 8);
hold(ax, 'on');
stairs(t, FTCS(FTC_ON).Yfails(3, :), 'Color', orange, 'LineWidth', 1.5)
stairs(t, FTCS(FTC_ON).UIOO(2).Fsen, '--', 'Color', blue, 'LineWidth', 1.5)
xlim(ax, [10 20]); ylim(ax, [-2.002 -1.98]);
box(ax, 'on'); grid(ax, 'on');

% Create axes
ax = axes('Parent', fig, 'Position', [0.225 0.657 0.284 0.211], 'FontSize', 8);
hold(ax, 'on');
stairs(t, FTCS(FTC_ON).Yfails(1, :), 'Color', orange, 'LineWidth', 1.5)
stairs(t, FTCS(FTC_ON).UIOO(1).Fsen, '--', 'Color', blue, 'LineWidth', 1.5)
xlim(ax, [49.5 53]); ylim(ax, [1 2.5]);
box(ax, 'on'); grid(ax, 'on');

% Create textarrow
annotation(fig, 'textarrow', [0.37 0.337], [0.728 0.728], ...
    'String', {'Threshold', 'effect'}, 'LineWidth', 1, ...
    'HorizontalAlignment', 'center', 'HeadWidth', 6, ...
    'HeadLength', 6, 'FontSize', 8);

print -dsvg ../Figs/FTCS_CSTR_sensor_estimation.svg

%% Input
fig = figure('Name', 'Inputs');
subplot(2, 1, 1)
plot(t, umin(1)*ones(1, length(t)), 'r--', 'LineWidth', 1.5)
hold on
stairs(t, FTCS(FTC_OFF).U(1, :), 'b', 'LineWidth', 1.5)
stairs(t, FTCS(FTC_ON).U(1, :), 'k--', 'LineWidth', 1.5);
plot(t, umax(1)*ones(1, length(t)), 'r--', 'LineWidth', 1.5);
xlabel('Time [min]'); ylabel('Qs [l/m]'); grid on; hold off;
xlim([0 Time])
leg = legend('Constraints', 'MPC', 'FT-MPC', 'Location', 'SouthEast');
leg.ItemTokenSize = [20, 15];
set(leg, 'Orientation', 'horizontal');
subplot(2, 1, 2)
plot(t, umin(2)*ones(1, length(t)), 'r--', t, umax(2)*ones(1, length(t)), 'r--', 'LineWidth', 1.5)
hold on
stairs(t, FTCS(FTC_OFF).U(2, :), 'b', 'LineWidth', 1.5)
stairs(t, FTCS(FTC_ON).U(2, :), 'k--', 'LineWidth', 1.5);
xlabel('Time [min]'); ylabel('Qc [l/m]'); grid on; hold off;
xlim([0 Time])

annotation(fig, 'textarrow', [0.698 0.726], [0.33 0.304], ...
    'String', {'FT-MPC fault', 'compensation'}, 'LineWidth', 1, 'HorizontalAlignment','center', ...
    'HeadWidth', 6, 'HeadLength', 6, 'FontSize', 8);
annotation(fig, 'textarrow', [0.235 0.261], [0.352 0.22], ...
    'String', {'MPC mistaken', 'fault compensation'}, 'LineWidth', 1, 'HorizontalAlignment', 'center',...
    'HeadWidth', 6, 'HeadLength', 6, 'FontSize', 8);

print -dsvg ../Figs/FTCS_CSTR_input.svg

%% Period consumption
fig = figure(10);
axes1 = axes('Parent', fig, 'Position', [0.13 0.314 0.775 0.442]);
area(t, (FTCS(FTC_ON).time_MHE+FTCS(FTC_ON).time_MPC+FTCS(FTC_ON).time_FDD)*100/Ts, 'FaceColor', green); hold(axes1, 'on');
area(t, FTCS(FTC_ON).time_MPC*100/Ts, 'FaceColor', yellow, 'FaceAlpha', 0.5);
area(t, FTCS(FTC_ON).time_MHE*100/Ts, 'FaceColor', blue, 'FaceAlpha', 1);
area(t, FTCS(FTC_ON).time_FDD*100/Ts, 'FaceColor', red, 'FaceAlpha', 0.7);
xlabel('Time [min]'); ylabel('Period [%]'); grid on
axis([0 Time 0 100]);
leg = legend('All', 'MPC', 'MHE', 'FDD');
set(leg, 'Location', 'NorthEast');
set(axes1, 'FontSize', 8);
leg.ItemTokenSize = [20, 10];

print -dsvg ../Figs/periodLoad.svg
