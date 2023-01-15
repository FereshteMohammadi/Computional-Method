function sys = add_translational_joint(sys, body_i, ...
    body_j, s_i, s_j)
    arguments
      sys (1,1) struct
        body_i (1,1) string
        body_j (1,1) string
        s_i (2,1) double = [0; 0]
        s_j (2,1) double = [0; 0]
    end
    check_body_exists(sys, body_i)
    check_body_exists(sys, body_j)

    %c_id = coordinate_name_to_id(coord_name);

   joint = struct();
    joint.body_i_qidx = body_name_to_qidx(sys, body_i);
    joint.body_j_qidx = body_name_to_qidx(sys, body_j);
    joint.s_i = s_i;
    joint.s_j = s_j;

    sys.joints.translational= [sys.joints.translational, joint];
end