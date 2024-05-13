function ux = TMcontract(t, Xs, Pmax, fend)
% input:
%   t  -- scalar, time
%   Xs -- (nip, 1) vector, integration points
%   Fmax -- scalar, maximum force applied to the cable
%   fend -- scalar, force end at fend
% returns:
%   ux -- (nip, 1) vector, cable tension at integration points at t

    hit = sqrt(-log(0.02)); % sigmoid function reach 98%

    T = 0.2; % ramping time
    Tp = 2; % propangation time

    nip = length(Xs); % number of integration points
    ux = zeros(nip, 1); % boundary stress at integration points at t.

    c = hit/T; % compressed factor
    X = t- Xs(Xs < fend)/fend*Tp;
    ux(Xs < fend) = Pmax*heaviside(X).*(1-exp(-(c*X).^2));
end