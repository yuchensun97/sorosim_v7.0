function Bq_rho = ComputeRadialActuation(Tr, rc, u_rho)
    % Compute the actuation loads for the cable
    % rc: (n_ract, 2) vector, radial actuation location
    % u_rho: (nip, 1), the actuation loads at intergration points
    % returns:
    % Bq_rho: ndof_rho x 1, the actuation loads in q-space
    Bq_rho = zeros(Tr.ndof_rho, Tr.n_ract);
    nsig = Tr.nsig;
    B_rho = Tr.Twists(2).B_rho;
    Xs = Tr.Twists(2).Xs;
    Ws = Tr.Twists(2).Ws;
    r_fn = Tr.Link.r_fn;
    ld = Tr.Link.L;

    j = 1; % index of rc intervals
    for ii = 1:nsig
        if Xs(ii)>=rc(j,1) && Xs(ii)<=rc(j,2)
            if Ws(ii) > 0
                % do the integration
                z = r_fn(Xs(ii));
                B_rho_here = B_rho(ii, :);
                Bq_rho(:,j) = Bq_rho(:,j)+2*pi*(ld*Ws(ii)*B_rho_here'*z^2*u_rho(ii));
            end
        end
        if Xs(ii) > rc(j, 2)
            if j == Tr.n_ract
                break;
            end
            j=j+1;
        end
    end
end
