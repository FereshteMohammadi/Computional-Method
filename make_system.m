function sys = make_system(gravity)
%MAKE_SYSTEM Create a data structure to store complete multibody system
arguments
    gravity (2,1) double = [0; -9.80665]
end
sys = struct();

sys.gravity = gravity;

sys.bodies = struct([]);

sys.joints = struct('revolute', struct([]), ...
    'simple', struct([]), ...
    'simple_driving', struct([]), ...
    'translational', struct([]));

sys.solver = struct('t_final', 1, 't_step', 0.01);

end