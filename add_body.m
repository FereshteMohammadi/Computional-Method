function sys = add_body(sys, name, location, orientation, mass, moment_of_inertia)
%ADD_BODY Add single body to the multibody system
    arguments
        sys (1,1) struct
        name (1,1) string
        location (2,1) double = [0; 0]
        orientation (1,1) double = 0 
        mass (1,1) double = 1 
        moment_of_inertia (1,1) double = 1 
    end
    body = struct("name", name);
    body.location = location;
    body.orientation = orientation;
    body.mass = mass;
    body.moment_of_inertia = moment_of_inertia;

    sys.bodies = [sys.bodies, body];
end