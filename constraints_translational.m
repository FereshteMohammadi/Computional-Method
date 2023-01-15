function C = constraints_translational(sys, q)

C = zeros(2 * length(sys.joints.translational), 1);
c_id = 0;

% for j = sys.joints.translational
%     qi = q(j.body_qidx);
%     C(c_id + 1) = qi(j.coord_id) - j.coord_value;
%     c_id = c_id + 1;
% end

for j = sys.joints.translational
    qi = q(j.body_i_qidx);
    qj = q(j.body_j_qidx);
    Ai = rot(qi(3));
    Aj = rot(qj(3));
    C(c_id + (1:2)) = qi(1:2) + Ai * j.s_i ...
        - qj(1:2) - Aj * j.s_j;
    c_id = c_id + 2;
end