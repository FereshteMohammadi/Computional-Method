function Cq = constraints_dq_revolute(sys, q)
%CONSTRAINTS_REVOLUTE Compute constraints for the revolute joints
Cq = zeros(2 * length(sys.joints.revolute), length(q));
c_id = 0;

O = Omega();
I = eye(2);
for j = sys.joints.revolute
    qi = q(j.body_i_qidx);
    qj = q(j.body_j_qidx);
    Ai = rot(qi(3));
    Aj = rot(qj(3));
    Cq(c_id + (1:2), j.body_i_qidx) = [I, O * Ai * j.s_i];
    Cq(c_id + (1:2), j.body_j_qidx) = -[I, O * Aj * j.s_j];
    c_id = c_id + 2;
end

end