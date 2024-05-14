function ux = LMrelease(t, Xs, Fmax, fend)
    % Inputs:
    %   t    -- scalar, time
    %   Xs   -- (nip, 1) vector, integration points
    %   Fmax -- scalar, maximum force applied to the cable
    %   fend -- scalar, force end at fend (equivalent to 'a' in the initial description)

    % Constants
    T = 1.5;  % Ramping time (time it takes for the force to drop to 0)
    Tp = 5;   % Propagation time (time it takes for the effect to propagate along the cable)

    % Calculate start and end time of force decay at each point
    startTime = (Xs / fend) * Tp;
    endTime = startTime + (1 - Xs / fend) * T;

    % Initialize cable tension at integration points to zero
    ux = zeros(length(Xs), 1);

    % Apply force calculations only where applicable
    for i = 1:length(Xs)
        if Xs(i) > fend
            ux(i) = 0;  % Ensure that for x > fend, u(x, t) is always zero
        elseif t < startTime(i)
            ux(i) = (1 - Xs(i) / fend) * Fmax;  % Before decay starts
        elseif t >= startTime(i) && t <= endTime(i)
            ux(i) = (1 - (t - startTime(i)) / (endTime(i) - startTime(i))) * (1 - Xs(i) / fend) * Fmax;
        else
            ux(i) = 0;  % After decay has finished
        end
    end
    ux = -ux;
end
