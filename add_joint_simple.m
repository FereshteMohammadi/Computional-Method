function sys = add_joint_simple(sys, body_name, coord_name, value)
    arguments
        sys (1,1) struct
        body_name (1,1) string
        coord_name (1,1) string
        value (1,1) double = 0
    end
    check_body_exists(sys, body_name)

    c_id = coordinate_name_to_id(coord_name);

    j = struct();
    j.body_qidx = body_name_to_qidx(sys, body_name);
    j.coord_id = c_id;
    j.coord_value = value;

    sys.joints.simple = [sys.joints.simple, j];
end