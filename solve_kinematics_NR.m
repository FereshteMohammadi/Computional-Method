function [T, Q, Qd] = solve_kinematics_NR(sys)
%SOLVE_KINEMATICS_FSOLVE Solve the multibody system sys on kinematics using
%Newton-Raphson method
% Solve the system based on the system definition and solver settings
% provided in sys
% Returns
% T - vector of time where solution is calculated
% Q - matrix with coordinates at each time in T. Coordinates are stored
% in columns

q = initial_coordinates(sys);

n_steps = ceil(sys.solver.t_final / sys.solver.t_step) + 1;
T = linspace(0, sys.solver.t_final, n_steps);

Q = zeros(length(q), length(T));
Qd = zeros(length(q), length(T));
for ii = 1:length(T)
    t = T(ii);
    [q, iteration_counter] = NewtonRaphson_method(...
        @(q) constraints(sys, q, t), ...
        @(q) constraints_dq(sys, q), ...
        q, ...
        1e-8);
    if iteration_counter == -1
        error("Newton-Raphson method did not converged at time %g", T(ii))
    end
    Q(:, ii) = q;
    % velocity level analysis
    qd = -constraints_dq(sys, q) \ constraints_dt(sys, t);
    Qd(:, ii) = qd;
end

end