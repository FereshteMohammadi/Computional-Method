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
X = [pi; 0]; % X_0 = [theta_0, d_0];
i = 1;

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
figure(3)
subplot(2,1,1)
plot(time, theta);
legend('\Theta');
subplot(2,1,2);
plot(time, d);
legend('displacement');
xlabel('Time(s)');

figure(2)
subplot(2,1,1)
plot(time, d_theta);
legend('d_\Theta');
subplot(2,1,2);
plot(time, d_d);
legend('d_displacement');
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


