function ux = TMact(t, Xs, Pmax, Tp, bp_s, bp_e)
    % Inputs:
        %   t    -- scalar, time
        %   Xs   -- (nip, 1) vector, integration points
        %   Fmax -- scalar, maximum force applied to the cable
        %   Tp -- scalar, propagation time
        %   bp_s -- scalar, initial bend point
        %   bp_e -- scalar, final bend point
    % Outputs:
        %   ux -- (nip,1) vector, LM tension at integration points at t
    
        if t < Tp
            mu = (bp_e-bp_s)/Tp * t + bp_s; % propangation function
            alpha = Pmax;
        else
            mu = bp_e;
            alpha = Pmax;
        end
    
        ux = alpha * (1- 1./(1+exp(-200*(Xs-mu))));
        ux = -ux;
    end
    