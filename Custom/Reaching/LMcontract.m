function ux = LMcontract(t, Xs, Fmax, Tp, bp_s, bp_e)
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
        alpha = Fmax;
    else
        mu = bp_e;
        alpha = Fmax;
    end

    % ux = alpha * exp(-(Xs-mu).^2/(2*0.17^2));
    ux = alpha * (1- 1./(1+exp(-200*(Xs-mu))));
    ux = -ux;
end
