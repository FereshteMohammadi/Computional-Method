function Cq = constraints_dq_simple_driving(sys, q)

Cq = zeros(length(sys.joints.simple_driving), length(q));
c_id = 0;

for j = sys.joints.simple_driving
    Cq(c_id + 1, j.body_qidx(j.coord_id)) = 1;
    c_id = c_id + 1;
end