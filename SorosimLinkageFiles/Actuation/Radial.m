classdef Radial
    properties(Access=private)
        r_pos    % position of the radial actuator, interval \belongs [0, 1]
    end

    methods
        function R = Radial(r_start, r_end)
            if r_start < 0 || r_start > 1
                error("r_start must be in [0, 1]");
            end

            if r_end < 0 || r_start > 1
                error("r_end must be in [0, 1]");
            end
            
            if r_end <= r_start
                error("r_end must be greater than r_start");
            end

            R.r_pos = [r_star r_end];
        end
    end

    methods
        function r_pos = get_r_pos(R)
            r_pos = R.r_pos;
        end
    end
end
