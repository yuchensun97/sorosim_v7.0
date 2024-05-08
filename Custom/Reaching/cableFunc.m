function ux = cableFunc(t, Xs)
% input:
%   t  -- scalar, time
%   Xs -- (nip, 1) vector, integration points
% returns:
%   ux -- (nip, 1) vector, cable tension at integration points at t

    Fmax = 15; % maximum force 15 N
    T = 0.1; % ramping time
    Tp = 2; % propangation time

    pe = (t-T)/(Tp-T); % propangation end
    ps = t/Tp; % propangation start

    nip = length(Xs); % number of integration points
    ux = zeros(nip, 1); % cable tension at integration points at t

    ux(Xs < pe) = 0;
    ux(Xs >= pe & Xs < ps) = (1-ps)*Fmax*(Xs(Xs >= pe & Xs < ps)-pe)/(ps-pe);
    ux(Xs >= ps & Xs < 1) = (1-Xs(Xs >= ps & Xs < 1))*Fmax;
    ux(Xs >= 1) = 0;
end
