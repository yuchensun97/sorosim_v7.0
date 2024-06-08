function ux = LMrelease(t, Xs, Fmax, Fmin, Tp, bp_s, bp_e)
    % Inputs:
    %   t    -- scalar, time
    %   Xs   -- (nip, 1) vector, integration points
    %   Fmax -- scalar, maximum force applied to the cable
    %   Fmin -- scalar, minimum force applied to the cable
    %   Tp -- scalar, propagation time
    %   bp_s -- scalar, initial bend point
    %   bp_e -- scalar, final bend point
    % Outputs:
    %   ux -- (nip,1) vector, LM tension at integration points at t

    if t < Tp
        mu = (bp_e-bp_s)/Tp * t + bp_s; % propangation function
        alpha = (Fmax-Fmin)/Tp^2 * t^2 - 2 * (Fmax-Fmin)/Tp * t + Fmax;
        % alpha = -(Fmax-Fmin)/Tp*t+Fmax;
    else
        mu = bp_e;
        alpha = Fmin;
    end

    % ux = alpha * exp(-(Xs-mu).^2/(2*0.17^2));
    ux = alpha * (1- 1./(1+exp(-40*(Xs-mu))));
    ux = -ux;
end
