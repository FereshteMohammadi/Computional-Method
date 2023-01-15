function M = mass_matrix(sys)
%MASS_MATRIX Get mass matrix from the system
m = zeros(3 * length(sys.bodies(1:4)), 1); % mass matrix diagonal

for i = 1: length(sys.bodies)
    m(body_idx(i)) = [sys.bodies(i).mass; 
        sys.bodies(i).mass;
        sys.bodies(i).moment_of_inertia];
end
M = diag(m);