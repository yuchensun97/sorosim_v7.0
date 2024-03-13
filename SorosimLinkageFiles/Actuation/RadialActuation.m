classdef RadialActuation
    properties(Access=private)
        n_ract    %number of radial actuations
        rc        %(n_ract x 1) array of radial position
        local     %(n_ract x 1) array of if the actuator is local or not
    end

    methods % constructor
        function RA = RadialActuation(varargin)
            %% initializing radial actuators
            if nargin == 0
                RA.n_ract = 0;
                RA.rc = zeros(0, 0);
                RA.local = zeros(0, 0);
            else
                n_ract = nargin;
                rc = zeros(n_ract, 0);
                local = zeros(n_ract, 0);
                
                for i=1:n_ract
                    R = varargin{i};
                    if ~isa(R, 'Radial')
                        error('Input must be of Radial class');
                    else
                        rc(i, 1) = R.get_rc();
                        local(i, 1) = R.get_local();
                    end
                end

                RA.n_ract = n_ract;
                RA.rc = rc;
                RA.local = local;
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
        
        function local = get_local(RA)
            local = RA.local;
        end
    end
end
