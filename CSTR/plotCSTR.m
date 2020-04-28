%% plot LPV-MPC on CSTR
% clc; clear; close all;
% load('run-3k6.mat')

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

error = abs(X(2, :) - xsp(2));
IAE = trapz(tsim, abs(error));
ISE = trapz(tsim, error.^2);
ITAE = trapz(tsim, tsim.*abs(error));
msg = ['IAE = ', num2str(IAE)];
disp(msg)
msg = ['ITAE = ', num2str(ITAE)];
disp(msg)
msg = ['ISE = ', num2str(ISE)];
disp(msg)

time_avg = mean(elapsed_time) ;
msg = ['Mean time = ', num2str(time_avg)];
disp(msg)
time_avg = max(elapsed_time) ;
msg = ['Max time = ', num2str(time_avg)];
disp(msg)
time_avg = min(elapsed_time) ;
msg = ['Min time = ', num2str(time_avg)];
disp(msg)

%% Outputs
fig = figure(3);
subplot(311)
plot(tsim, Xsp(1, :), ':', 'Color', orange_red, 'LineWidth', 1.5);
hold on
plot(tsim, X(1, :), 'Color', chocolate, 'LineWidth', 1.5);
xlim([0 Time])
hold off; grid on
xlabel('Time [min]'); ylabel('V [l]');
subplot(312)
plot(tsim, Xsp(2, :), ':', 'Color', orange_red, 'LineWidth', 1.5);
hold on
plot(tsim, X(2, :), '-', 'Color', azul, 'LineWidth', 1.5);
xlim([0 Time])
hold off; grid on

subplot(313)
plot(tsim, Xsp(3, :), ':', 'Color', orange_red, 'LineWidth', 1.5);
hold on
plot(tsim, X(3, :), '-', 'Color', azul, 'LineWidth', 1.5);
xlim([0 Time])
hold off; grid on
xlabel('Time [min]'); ylabel('T [°K]');

print -dsvg figs/state.svg

%% State space
figure(2)
plot(X(3, :), X(2, :), '-', 'Color', azul, 'LineWidth', 1.5);
hold on
plot(x0(3), x0(2), 'o', 'Color', verde, 'LineWidth', 1.5);
plot(Xsp(3, :), Xsp(2, :), '*', 'Color', bordo, 'LineWidth', 1.5);
hold off

%% Input
figure(4)
subplot(2, 1, 1)
stairs(t, umpc(1, :), 'Color', orange_red, 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('u_1 [l/m]'); grid on
xlim([0 Time])
subplot(2, 1, 2)
stairs(t, umpc(2, :), 'Color', orange_red, 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('u_1 [l/m]'); grid on
xlim([0 Time])

print -dsvg figs/input.svg

%% Objective function
figure(6)
plot(t, obj(:))
xlim([0 Time])
xlabel('Muestra'); ylabel('objective');

%% Membership
figure(31)
hold on
for l=1:M
    plot(t, mu_mhe(l, :));
    legendInfo{l} = ['\mu' num2str(l)];
    plot(t, mu_fuzzy(l, :), '--');
    legendInfo{l+1} = ['\mu' num2str(l+1)];
end
ylabel('\mu'), xlabel('Time [s]')
legend(legendInfo);
xlim([0 Time])
% title('Pertenencia')
hold off