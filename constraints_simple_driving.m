function C = constraints_simple_driving(sys, q, t)

C = zeros(length(sys.joints.simple_driving), 1);
c_id = 0;

for j = sys.joints.simple_driving
    qi = q(j.body_qidx);
    C(c_id + 1) = qi(j.coord_id) - j.coord_fun(t);
    c_id = c_id + 1;
end