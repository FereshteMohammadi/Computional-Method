function b_id = body_name_to_id(sys, body_name)

b_id = 1;
while sys.bodies(b_id).name ~= body_name
    b_id = b_id + 1;
end