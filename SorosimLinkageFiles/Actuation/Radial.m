classdef Radial
    properties(Access=private)
        rc    % X coordinates of the radial actuation
        local % true if the pressure is applied to local frame
    end

    methods
        function R = Radial(cx, local)
            if cx<0 || cx > 1
                error('cx must be in [0, 1]');
            end
            R.rc = cx;
            R.local = local;
        end
        
        %% Getter methods
        function cx = get_rc(R)
            cx = R.rc;
        end
        
        function local = get_local(R)
            local = R.local;
        end
    end
end
