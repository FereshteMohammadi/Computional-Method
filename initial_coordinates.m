function q = initial_coordinates(sys)
%INITIAL_COORDINATES Get initial coordinates from the system
q = zeros(3 * length(sys.bodies), 1);
for i = 1:length(sys.bodies)
    q(body_idx(i)) = [sys.bodies(i).location; 
        sys.bodies(i).orientation];
end