function ux = OM(t, Xs, Fmax, Ts, Tp, bp)
% Inputs:
    %   t    -- scalar, time
    %   Xs   -- (nip, 1) vector, integration points
    %   Fmax -- scalar, maximum force applied to the cable
    %   Ts -- scalar, start time
    %   Tp -- scalar, propagation end time
    %   bp -- scalar, initial bend point
    % Outputs:
    %   ux -- (nip,1) vector, LM tension at integration points at t

    if t < Ts
        mu = 0;
        alpha = 0;
    elseif t < Tp
        mu = bp/(Tp-Ts) * (t-Ts);
        alpha = Fmax;
    else
        mu = bp;
        alpha = Fmax;
    end

    ux = alpha * (1- 1./(1+exp(-200*(Xs-mu))));
    ux = -ux;
end
