function sys = set_solver_settings(sys, t_final, t_step)
arguments
    sys (1,1) struct
    t_final (1,1) double {mustBePositive} = 1.0
    t_step (1,1) double {mustBePositive} = 0.01
end

sys.solver.t_final = t_final;
sys.solver.t_step = t_step;