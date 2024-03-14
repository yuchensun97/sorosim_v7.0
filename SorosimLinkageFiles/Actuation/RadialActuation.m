classdef RadialActuation
    properties(Access=private)
        n_ract    %number of radial actuations
        rc        %(n_ract x 2) array of radial actuator position
    end

    methods % constructor
        function RA = RadialActuation(varargin)
            %% initializing radial actuators
            if nargin == 0
                RA.n_ract = 0;
                RA.rc = zeros(0, 2);
            else
                n_ract = nargin;
                rc = zeros(n_ract, 2);
                
                for i=1:n_ract
                    R = varargin{i};
                    if ~isa(R, 'Radial')
                        error('Input must be of Radial class');
                    else
                        rc(i, :) = R.get_r_pos;
                    end
                end
                if RadialActuation.checkOverlay(rc)
                    error('Overlay of radial actuations is not allowed');
                end
                RA.n_ract = n_ract;
                RA.rc = rc;
            end
        end
    end

    methods(Access=private, Static=true)
        function overlay = checkOverlay(intervals)
            sorted_intervals = sortrows(intervals, 1);
            overlay = false;
            n = size(sorted_intervals, 1)-1;
            for i=1:n
                curr_end = sorted_intervals(i, 2);
                next_start = sorted_intervals(i+1, 1);
                if curr_end > next_start
                    overlay = true;
                    break;
                end
            end
        end
    end

    methods % getter
        function n = get_n_sact(RA)
            n = RA.n_ract;
        end
        
        function rc = get_rc(RA)
            rc = RA.rc;
        end
    end
end
