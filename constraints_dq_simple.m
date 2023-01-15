function Cq = constraints_dq_simple(sys, q)

Cq = zeros(length(sys.joints.simple), length(q));
c_id = 0;

for j = sys.joints.simple
    Cq(c_id + 1, j.body_qidx(j.coord_id)) = 1;
    c_id = c_id + 1;
end