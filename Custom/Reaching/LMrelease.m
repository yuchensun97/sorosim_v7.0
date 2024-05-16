function ux = LMrelease(t, Xs, Fmax, Tp, bp_s, bp_e)
    % Inputs:
    %   t    -- scalar, time
    %   Xs   -- (nip, 1) vector, integration points
    %   Fmax -- scalar, maximum force applied to the cable
    %   Tp -- scalar, propagation time
    %   bp_s -- scalar, initial bend point
    %   bp_e -- scalar, final bend point
    % Outputs:
    %   ux -- (nip,1) vector, LM tension at integration points at t

    if t<Tp
        beta = Fmax/Tp*(Tp-t);
        alpha = bp_s + (bp_e-bp_s)/Tp*t;
    else
        beta = 0;
        alpha = bp_s;
    end

    ux = beta * (1 - 1./(1+exp(-50*(Xs - alpha))));
    ux = -ux';
end
