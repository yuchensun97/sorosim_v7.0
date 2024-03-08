classdef Cable
    properties(Access=private)
        dc_fn  % 3d coordinates of the cable as a function of X
        dcp_fn % 3d coordinates of the cable as a function of X
    end

    methods
        function C = Cable(cy_fn, cz_fn)
            if ~isa(cy_fn, 'function_handle')
                error('cy_fn should be a function handle');
            end
            if ~isa(cz_fn, 'function_handle')
                error('cz_fn should be a function handle');
            end
            cy_fn_s = func2str(cy_fn);
            if ~startsWith(cy_fn_s, '@(X)')
                error('cy_fn should be the form of @(X)...');
            end
            cz_fn_s = func2str(cz_fn);
            if ~startsWith(cz_fn_s, '@(X)')
                error('cz_fn should be the form of @(X)...');
            end
            cyp_fn = @(X)derivative(cy_fn, X);
            czp_fn = @(X)derivative(cz_fn, X);
            dc_fn = @(X)[0;cy_fn(X);cz_fn(X)];
            dcp_fn = @(X)[0;cyp_fn(X);czp_fn(X)];
            C.dc_fn = dc_fn;
            C.dcp_fn = dcp_fn;
        end
    end

    methods(Access = private, Static=true) % statics method
        function dydx = derivative(f, X)
            % f, function handler
            % X, point at which to compute the derivative
            h = 1e-5;
            dydx = (f(X+h)-f(X))/h;
        end
    end

    methods
        function dc = get_dc_fn(C)
            dc = C.dc_fn;
        end
        
        function dcp = get_dcp_fn(C)
            dcp = C.dcp_fn;
        end
    end
end
