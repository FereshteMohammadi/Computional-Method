function check_body_exists(sys, name)
    arguments
        sys (1,1) struct
        name (1,1) string
    end

    for b = sys.bodies
        if strcmp(b.name, name)
            return
        end
    end
    error("Body %s not found in the system.", name)
end