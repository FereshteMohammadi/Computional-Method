function [T, Q] = solve_kinematics_fsolve(sys)
%SOLVE_KINEMATICS_FSOLVE Solve the multibody system sys on kinematics using
%fsolve method
% Solve the system based on the system definition and solver settings
% provided in sys
% Returns
% T - vector of time where solution is calculated
% Q - matrix with coordinates at each time in T. Coordinates are stored
% in columns

q = initial_coordinates(sys);
opts = optimoptions("fsolve", "Display", "none");

n_steps = ceil(sys.solver.t_final / sys.solver.t_step) + 1;
T = linspace(0, sys.solver.t_final, n_steps);

Q = zeros(length(q), length(T));
for ii = 1:length(T)
    [q,~,exitflag,~] = fsolve(@(q) constraints(sys, q, T(ii)), q, opts);
    if exitflag < 1
        error("Check the fsolve solution as it exits with flag %d", exitflag)
    end
    Q(:, ii) = q;
end

end