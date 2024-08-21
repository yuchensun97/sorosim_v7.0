function ux = LM(t, Xs, Fmax, Fmin, Tp, bp, be)
% Inputs:
    %   t    -- scalar, time
    %   Xs   -- (nip, 1) vector, integration points
    %   Fmax -- scalar, maximum force applied to the cable
    %   Fmin -- scalar, minimum force applied to the cable
    %   Tp -- scalar, propagation time
    %   bp -- scalar, initial bend point
    %   be -- scalar, final bend point
    % Outputs:
    %   ux -- (nip,1) vector, LM tension at integration points at t

    if t < Tp
        mu = (bp-be) / Tp * t + be;
        alpha = (Fmax - Fmin)/Tp * t + Fmin;
    else
        mu = bp;
        alpha = Fmax;
    end

    ux = alpha * (1- 1./(1+exp(-200*(Xs-mu))));
    ux = -ux;
end
