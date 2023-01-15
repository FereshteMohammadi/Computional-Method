function b_idx = body_name_to_qidx(sys, body_name)

b_id = body_name_to_id(sys, body_name);
b_idx = body_idx(b_id);
end