clc
clear
close all

trials = 250; % number of trials

% row 1: PI controller
% row 2: Optimal Controller with Observer
% row 3: Optimal Controller with Time-Varying Kalman Filter
% row 4: Optimal Controller with Steady-State Kalman Filter
% row 5: Model Predictive Controller
muE = zeros([5, trials]);
stdE = zeros([5, trials]);

% Optimal controller parameters
Q = 1;
R = 1;

% Standard Deviation of Gaussian noise
sigma = sqrt(0.1);
Wvar = sigma^2;
Vvar = sigma^2;
nu = 1;

for i = 1:trials
disp("i = " + i);
N = 100; % number of time steps

E = zeros([1 N]); % energy (in kWh) stored over time
E(1) = 0; % initial energy stored in kWh
Es = zeros([1 N]); % energy provided by solar system into battery (in kWh)

% Input and output constraints
Esmin = -2;
Esmax = 2;
dEsmax = 2;
E_min = 1;
E_max = 6;

% PI CONTROLLER

% PI controller gains
Kp = 1;
Ki = 1;
Ed = 3; % setpoint
intD = 0; % integral of difference
for k = 1:(N - 1)
    D = Ed - E(k); % difference
    intD = intD + D; % integral of difference

    Es(k) = Kp*D + Ki*intD;
    El = nu+randn*sigma; % disturbance/noise

    E(k + 1) = E(k) + Es(k) - El;
end

% Mean Absolute Error (MAE) and Root Mean Square Error (RMSE)
muE(1, i) = mean(abs(E - Ed));
stdE(1, i) = sqrt(mean((E - Ed).^2));

if (i == 1)
% PI Control Output Plot
figure
plot(1:N, E);
title("PI Control (Variance = " + sigma^2 + "), MATLAB");
grid on
axis([1, N, 0, 7])
hold on
plot([1, N], [E_min, E_min]);
plot([1, N], [E_max, E_max]);
plot([1, N], [Ed, Ed]);
legend(["E(k)", "E_{min}", "E_{max}", "E_{set}"], 'location', 'best');
xlabel("k")
ylabel("Energy [kWh]")
hold off

% PI Control Input Plot
figure
plot(1:N, Es);
title("PI Control Input (Wvar = " + Wvar + ", Wvar = " + Vvar + "), MATLAB");
grid on
hold on
axis([1, N, 1.1*Esmin, 1.1*Esmax])
plot([1, N], [Esmin, Esmin]);
plot([1, N], [Esmax, Esmax]);
legend(["Es(k)", "Es_{min}", "Es_{max}"], 'location', 'best');
xlabel("k")
ylabel("Input Energy [kWh]")
hold off

% PI Control Delta-Input Plot
figure
plot(2:N, abs(Es(2:end) - Es(1:(end-1))));
title("PI Control Input Rate of Change (Wvar = " + Wvar + ", Wvar = " + Vvar + "), MATLAB");
grid on
hold on
axis([2, N, 0, 1.1*dEsmax])
plot([2, N], [dEsmax, dEsmax]);
legend(["\DeltaEs(k)", "\DeltaEs_{max}"], 'location', 'best');
xlabel("k")
ylabel("Rate of Change of Input Energy [kWh]")
hold off
end

% Optimal Control with Observer (OCO)
% Observer gain
L = 1;

Sinf = (Q + sqrt(Q^2 + 4*Q*R))/2;
pinf = -Ed * Q * (R + Sinf)/Sinf;
Kinf = Sinf/(R + Sinf);
Kcinf = -pinf/(R + Sinf);

E_opt = zeros([1, N]);
E_opt(1) = E(1);

E_opt_hat = zeros([1, N]);
E_opt_hat(1) = E(1);

Es_opt = zeros([1, N]);
for k = 1:(N - 1)
    Es_opt(k) = -Kinf * E_opt_hat(k) + Kcinf;
    El = nu+randn*sigma;

    E_opt(k + 1) = E_opt(k) + Es_opt(k) - El;
    E_opt_hat(k + 1) = E_opt_hat(k) + Es_opt(k) + L*(E_opt(k) - E_opt_hat(k));
end

% Mean Absolute Error (MAE) and Root Mean Square Error (RMSE)
muE(2, i) = mean(abs(E_opt - Ed));
stdE(2, i) = sqrt(mean((E_opt - Ed).^2));

if(i == 1)
% OCO Output Plot
figure
plot(1:N, E_opt);
title("Optimal Control with Observer (OCO) (Variance = " + sigma^2 + "), MATLAB");
grid on
axis([1, N, 0, 7])
hold on
plot([1, N], [E_min, E_min]);
plot([1, N], [E_max, E_max]);
plot([1, N], [Ed, Ed]);
legend(["E(k)", "E_{min}", "E_{max}", "E_{set}"], 'location', 'best');
xlabel("k")
ylabel("Energy [kWh]")
hold off

% OCO Input Plot
figure
plot(1:N, Es_opt);
title("OCO Input (Wvar = " + Wvar + ", Wvar = " + Vvar + "), MATLAB");
grid on
hold on
axis([1, N, 1.1*Esmin, 1.1*Esmax])
plot([1, N], [Esmin, Esmin]);
plot([1, N], [Esmax, Esmax]);
legend(["Es(k)", "Es_{min}", "Es_{max}"], 'location', 'best');
xlabel("k")
ylabel("Input Energy [kWh]")
hold off

% OCO Delta-Input Plot
figure
plot(2:N, abs(Es_opt(2:end) - Es_opt(1:(end-1))));
title("OCO Input Rate of Change (Wvar = " + Wvar + ", Wvar = " + Vvar + "), MATLAB");
grid on
hold on
axis([2, N, 0, 1.1*dEsmax])
plot([2, N], [dEsmax, dEsmax]);
legend(["\DeltaEs(k)", "\DeltaEs_{max}"], 'location', 'best');
xlabel("k")
ylabel("Rate of Change of Input Energy [kWh]")
hold off
end

% Optimal Control with Time-varying Kalman Filter (OCTVKF)
P = zeros([N 2]); % column 1 = a priori, column 2 = a posteriori
Kk = zeros([1 N]);
E_opt2 = zeros([1, N]);
E_opt2(1) = E(1);
E_hat = zeros([N, 2]); % column 1 = a priori, column 2 = a posteriori
E_hat(1, 2) = E(1);
y = zeros([1, N]);

Es_opt = zeros([1 N]);
for k = 2:N
    % Kalman filtering
    P(k, 1) = P(k-1, 2) + Wvar;
    Kk(k) = P(k, 1)/(Vvar + P(k, 1));
    P(k, 2) = (1 - Kk(k))*P(k, 1);

    Es_opt(k-1) = -Kinf*E_hat(k-1, 2) + Kcinf;
    E_opt2(k) = E_opt2(k-1) + Es_opt(k-1) + nu+randn*sqrt(Wvar);

    % More Kalman filtering
    y(k) = E_opt2(k) + nu+randn*sqrt(Vvar);
    E_hat(k, 1) = E_hat(k-1, 2) + Es_opt(k-1);
    E_hat(k, 2) = E_hat(k, 1) + Kk(k)*(y(k) - E_hat(k, 1));
end

% Mean Absolute Error (MAE) and Root Mean Square Error (RMSE)
muE(3, i) = mean(abs(E_opt2 - Ed));
stdE(3, i) = sqrt(mean((E_opt2 - Ed).^2));

if (i == 1)
% OCTVKF Output Plot
figure
plot(1:N, E_opt2);
title("Optimal Control with Time-Varying Kalman Filter (Wvar = " + Wvar + ", Wvar = " + Vvar + "), MATLAB");
grid on
axis([1, N, 0, 7])
hold on
plot([1, N], [E_min, E_min]);
plot([1, N], [E_max, E_max]);
plot([1, N], [Ed, Ed]);
legend(["E(k)", "E_{min}", "E_{max}", "E_{set}"], 'location', 'best');
xlabel("k")
ylabel("Energy [kWh]")
hold off

% OCTVKF Input Plot
figure
plot(1:N, Es_opt);
title("OCTVKF Input (Wvar = " + Wvar + ", Wvar = " + Vvar + "), MATLAB");
grid on
hold on
axis([1, N, 1.1*Esmin, 1.1*Esmax])
plot([1, N], [Esmin, Esmin]);
plot([1, N], [Esmax, Esmax]);
legend(["Es(k)", "Es_{min}", "Es_{max}"], 'location', 'best');
xlabel("k")
ylabel("Input Energy [kWh]")
hold off

% OCTVKF Delta-Input Plot
figure
plot(2:N, abs(Es_opt(2:end) - Es_opt(1:(end-1))));
title("OCTVKF Input Rate of Change (Wvar = " + Wvar + ", Wvar = " + Vvar + "), MATLAB");
grid on
hold on
axis([2, N, 0, 1.1*dEsmax])
plot([2, N], [dEsmax, dEsmax]);
legend(["\DeltaEs(k)", "\DeltaEs_{max}"], 'location', 'best');
xlabel("k")
ylabel("Rate of Change of Input Energy [kWh]")
hold off
end

% Optimal Control with Steady-State Kalman Filter (OCSSKF)
Pinf = (Wvar + sqrt(Wvar + 4*Vvar))/2;
Kkinf = Pinf/(Vvar + Pinf);

E_opt3 = zeros([1, N]);
E_opt3(1) = E(1);

Es_opt = zeros([1 N]);
for k = 2:N
    Es_opt(k-1) = -Kinf*E_hat(k-1, 2) + Kcinf;
    E_opt3(k) = E_opt3(k-1) + Es_opt(k-1) + nu+randn*sqrt(Wvar);
    y(k) = E_opt3(k) + nu+randn*sqrt(Vvar);

    % Steady-state Kalman filtering
    E_hat(k, 1) = E_hat(k-1, 2) + Es_opt(k-1);
    E_hat(k, 2) = E_hat(k, 1) + Kkinf*(y(k) - E_hat(k, 1));
end

% Mean Absolute Error (MAE) and Root Mean Square Error (RMSE)
muE(4, i) = mean(abs(E_opt3 - Ed));
stdE(4, i) = sqrt(mean((E_opt3 - Ed).^2));

if (i == 1)
% OCSSKF Output Plot
figure
plot(1:N, E_opt3);
title("Optimal Control with Steady-State Kalman Filter (Wvar = " + Wvar + ", Wvar = " + Vvar + "), MATLAB");
grid on
axis([1, N, 0, 7])
hold on
plot([1, N], [E_min, E_min]);
plot([1, N], [E_max, E_max]);
plot([1, N], [Ed, Ed]);
legend(["E(k)", "E_{min}", "E_{max}", "E_{set}"], 'location', 'best');
xlabel("k")
ylabel("Energy [kWh]")
hold off

% OCSSKF Input Plot
figure
plot(1:N, Es_opt);
title("OCSSKF Input (Wvar = " + Wvar + ", Wvar = " + Vvar + "), MATLAB");
grid on
hold on
axis([1, N, 1.1*Esmin, 1.1*Esmax])
plot([1, N], [Esmin, Esmin]);
plot([1, N], [Esmax, Esmax]);
legend(["Es(k)", "Es_{min}", "Es_{max}"], 'location', 'best');
xlabel("k")
ylabel("Input Energy [kWh]")
hold off

% OCSSKF Delta-Input Plot
figure
plot(2:N, abs(Es_opt(2:end) - Es_opt(1:(end-1))));
title("OCSSKF Input Rate of Change (Wvar = " + Wvar + ", Wvar = " + Vvar + "), MATLAB");
grid on
hold on
axis([2, N, 0, 1.1*dEsmax])
plot([2, N], [dEsmax, dEsmax]);
legend(["\DeltaEs(k)", "\DeltaEs_{max}"], 'location', 'best');
xlabel("k")
ylabel("Rate of Change of Input Energy [kWh]")
hold off
end

% Model Predictive Controller (MPC)
Ts = 0.01; % Sample Time
mInfo = [Ts, 1, 1, 0, 1, 1, 0]; % first row of the system in MOD format
phi = 1;
gamma = [1, -1];
C = 1;
D = [0, 0];
model = ss2mod(phi, gamma, C, D, mInfo);

r = Ed*ones([N+1, 1]); % command
Hu = 3; % control horizon
Hp = 5; % prediction horizon
TF = N*Ts; % final time

% input and output constraints
Ucon = [Esmin*ones([N+1, 1]), Esmax*ones([N+1, 1]), dEsmax*ones([N+1, 1])];
Ycon = [E_min*ones([N+1, 1]), E_max*ones([N+1, 1])];

% estimator gain
Kest = smpcest(model, Wvar, Vvar);

% Gaussian noise disturbance
w = sqrt(Wvar).*(randn([N+1, 1]))+nu;

% MPC simulation
[yp, u, ~] = scmpc(model, model, Q, R, Hu, Hp, TF, r, Ucon, Ycon, Kest, 0, [], w);
%clc
%{
M = min(size(yp, 1), N+1);
yp = yp(1:M);
u = u(1:M);
%}

% Mean Absolute Error (MAE) and Root Mean Square Error (RMSE)
muE(5, i) = mean(abs(yp - Ed));
stdE(5, i) = sqrt(mean((yp - Ed).^2));

if (i == 1)
% MPC Output Plot
figure
plot(0:N, [yp, r]);
title("Model Predictive Controller (Wvar = " + Wvar + ", Wvar = " + Vvar + "), MATLAB");
grid on
axis([1, N, 0, 7])
hold on
plot([0, N], [E_min, E_min]);
plot([0, N], [E_max, E_max]);
legend(["E(k)", "E_{set}", "E_{min}", "E_{max}"], 'location', 'best');
xlabel("k")
ylabel("Energy [kWh]")
hold off

% MPC Input Plot
figure
plot(0:N, u);
title("MPC Control Input (Wvar = " + Wvar + ", Wvar = " + Vvar + "), MATLAB");
grid on
hold on
axis([0, N, 1.1*Esmin, 1.1*Esmax])
plot([0, N], [Esmin, Esmin]);
plot([0, N], [Esmax, Esmax]);
legend(["Es(k)", "Es_{min}", "Es_{max}"], 'location', 'best');
xlabel("k")
ylabel("Input Energy [kWh]")
hold off

% MPC Input Rate of Change Plot
figure
plot(1:N, abs(u(2:end) - u(1:(end-1))));
title("MPC Input Rate of Change (Wvar = " + Wvar + ", Wvar = " + Vvar + "), MATLAB");
grid on
hold on
axis([1, N, 0, 1.1*dEsmax])
plot([1, N], [dEsmax, dEsmax]);
legend(["\DeltaEs(k)", "\DeltaEs_{max}"], 'location', 'best');
xlabel("k")
ylabel("Rate of Change of Input Energy [kWh]")
hold off
end

end

% Mean Errors
fprintf("Average MAE for PI Controller: %f\n", mean(muE(1, :)));
fprintf("Average MAE for Optimal Controller with Observer: %f\n", mean(muE(2, :)));
fprintf("Average MAE for Optimal Controller with Time-Varying Kalman Filter: %f\n", mean(muE(3, :)));
fprintf("Average MAE for Optimal Controller with Steady-State Kalman Filter: %f\n", mean(muE(4, :)));
fprintf("Average MAE for Model Predictive Controller: %f\n\n", mean(muE(5, :)));

% RMSEs
fprintf("Average RMSE for PI Controller: %f\n", mean(stdE(1, :)));
fprintf("Average RMSE for Optimal Controller with Observer: %f\n", mean(stdE(2, :)));
fprintf("Average RMSE for Optimal Controller with Time-Varying Kalman Filter: %f\n", mean(stdE(3, :)));
fprintf("Average RMSE for Optimal Controller with Steady-State Kalman Filter: %f\n", mean(stdE(4, :)));
fprintf("Average RMSE for Model Predictive Controller: %f\n", mean(stdE(5, :)));