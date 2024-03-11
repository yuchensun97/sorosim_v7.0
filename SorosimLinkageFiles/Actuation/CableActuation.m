classdef CableActuation
    properties(Access = private)
        n_sact  % number of cable actuations
        dc_fn   % (n_sact x 1) array of local cable position function handler
        dcp_fn  % (n_sact x 1) array of space derivative of the local cable position
    end

    methods % constructor
        function CA = CableActuation(varargin)
            %% initializing Cable vector
            if nargin == 0
                CA.n_sact = 0;
                CA.dc_fn = cell(0,0);
                CA.dcp_fn = cell(0, 0);
            else
                n_sact = nargin;
                dc_fn = cell(n_sact, 0);
                dcp_fn = cell(n_sact, 0);

                for i=1:n_sact
                    C = varargin{i};
                    if ~isa(C, 'Cable')
                        error('Input must be of Cable class');
                    else
                        dc_fn{i,1} = C.get_dc_fn();
                        dcp_fn{i,1} = C.get_dcp_fn();
                    end
                end

                CA.n_sact = n_sact;
                CA.dc_fn = dc_fn;
                CA.dcp_fn = dcp_fn;
            end
        end
    end

    methods
        % Getter for n_sact
        function n = get_n_sact(obj)
            n = obj.n_sact;
        end
        
        % Getter for dc_fn
        function dc = get_dc_fn(obj)
            dc = obj.dc_fn;
        end
        
        % Getter for dcp_fn
        function dcp = get_dcp_fn(obj)
            dcp = obj.dcp_fn;
        end
    end

end
