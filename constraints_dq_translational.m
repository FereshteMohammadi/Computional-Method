function Cq = constraints_dq_translational(sys, q)

Cq = zeros(2 * length(sys.joints.translational), length(q));
c_id = 0;

% for j = sys.joints.translational
%     Cq(c_id + 1, j.body_qidx(j.coord_id)) = 1;
%     c_id = c_id + 1;
% end

O = Omega();
I = eye(2);
for j = sys.joints.translational
    qi = q(j.body_i_qidx);
    qj = q(j.body_j_qidx);
    Ai = rot(qi(3));
    Aj = rot(qj(3));
    Cq(c_id + (1:2), j.body_i_qidx) = [I, O * Ai * j.s_i];
    Cq(c_id + (1:2), j.body_j_qidx) = -[I, O * Aj * j.s_j];
    c_id = c_id + 2;
end

end