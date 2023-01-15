function Ct = constraints_dt_simple_driving(sys, t)

Ct = zeros(length(sys.joints.simple_driving), 1);
c_id = 0;

for j = sys.joints.simple_driving
    Ct(c_id + 1) = -j.coord_fun_dt(t);
    c_id = c_id + 1;
end