function ux = LMrelease(t, Xs, Fmax, fend)
% input:
%   t  -- scalar, time
%   Xs -- (nip, 1) vector, integration points
%   Fmax -- scalar, maximum force applied to the cable
%   fend -- scalar, force end at fend
% returns:
%   ux -- (nip, 1) vector, cable tension at integration points at t

    T = 0.2; % ramping time
    Tp = 2; % propangation time

    pe = (t-T)/(Tp-T); % propangation end
    ps = t/Tp; % propangation start

    nip = length(Xs); % number of integration points
    ux = zeros(nip, 1); % cable tension at integration points at t

    ux(Xs < pe) = 0;
    ux(Xs >= pe & Xs < ps) = -(1-ps/fend)*Fmax*(Xs(Xs >= pe & Xs < ps)-pe)/(ps-pe);
    ux(Xs >= ps & Xs < fend) = -(1-Xs(Xs >= ps & Xs < fend)/fend)*Fmax;
    ux(Xs >= fend) = 0;
end
