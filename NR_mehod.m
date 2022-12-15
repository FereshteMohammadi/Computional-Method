%% CMiM_10_assignment
% Simple kinematic analysis
% Fereshte Mohammadi Khoshoue

close all;
clear;
clc;

%% Initial value
a = 0.1; %m
b = 0.2;
omega = 1; %rad/s
error = 1e-9;
time = linspace(0, 1, 101);
i = 1;
phi_0 = deg2rad(30); %rad

%% Initial value of theta and d
syms theta d 

F1 = a * cos(phi_0) + b * cos(theta) - d;
F2 = a * sin(phi_0) + b * sin(theta);

[theta, d] = solve(F1, F2, [theta, d]);

theta_0 = double(theta(1));
d_0 = double(d(1));

fprintf("theta = %g rad\nd = %g\n",theta_0,d_0);

X = [theta_0; d_0]; % X_0 = [theta_0, d_0];

%%
for t = linspace(0, 1, 101)

    phi = pi/6 + omega * t;
    F = @(X) constarint(X, a, b, phi);
    J = @(X) jacobian(X, b);

    [X, iteration_counter] = NewtonRaphson_method(F, J, X, error);

     dFt = [-a * omega * sin(omega * t)
             a * omega * cos(omega * t)];

    dFq = [-b * sin(X(1)), -1
           -b * cos(X(1)), 0];

    d_X = - dFt \ dFq;
    
    theta(i,1) = X(1);
    d(i,1) = X(2);
    d_theta(i,1) = d_X(1);
    d_d(i,1) = d_X(2);

    i = i + 1;
end

%% Plot
figure(1)
subplot(2,1,1)
plot(time, theta, 'LineWidth', 2);
ylabel('\Theta');
subplot(2,1,2);
plot(time, d, 'LineWidth', 2);
ylabel('d');
xlabel('Time(s)');

figure(2)
subplot(2,1,1)
plot(time, d_theta, 'LineWidth', 2);
ylabel('\Theta derivation');
subplot(2,1,2);
plot(time, d_d, 'LineWidth', 2);
ylabel('d derivation');
xlabel('Time(s)');

%% Constraint Matrix
function C = constarint(X, a, b, phi)
    theta = X(1);
    d = X(2);
    C = [a * cos(phi) + b * cos(theta) - d
         a * sin(phi) - b * sin(theta)];
end

%% Jacobian Matrix
function J = jacobian(X, b)
    theta = X(1);
    J = [-b * sin(theta), -1
         -b * cos(theta), 0]; 
end

%% Newton-Raphson method
function [x, iteration_counter] = NewtonRaphson_method(F, J, x, tol)
    F_value = F(x);
    F_norm = norm(F_value);
    iteration_counter = 0;

    while F_norm > tol && iteration_counter < 100
        delta = J(x)\-F_value;
        x = x + delta;
        F_value = F(x);
        F_norm = norm(F_value);
        iteration_counter = iteration_counter + 1;
    end

    if F_norm > tol
        iteration_counter = -1;
    end
end


