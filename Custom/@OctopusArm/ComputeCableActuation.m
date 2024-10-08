function Bq_xi = ComputeCableActuation(Tr, dc, dcp, q_xi, q_rho, u_xi)
    % Compute the actuation loads for the cable
    % dc: 3 x nip, the direction of the cable
    % dcp: 3 x nip, the derivative of the direction of the cable
    % q_xi: ndof_xi x 1, the generalized coordinates of the xi
    % q_rho: ndof_rho x 1, the generalized coordinates of the rho
    % u_xi: nip x 1, the actuation loads at intergration points
    % returns:
    % Bq_xi: ndof_xi x 1, the actuation loads in q-space
    Bq_xi = zeros(Tr.ndof_xi, 1);
    Twists = Tr.Twists(2);

    xi_star = Twists.xi_star;
    rho_star = Twists.rho_star;
    ld = Tr.Link.L;
    Ws = Twists.Ws;
    nip = Twists.nip;
    ndof_xi = Twists.dof_xi;
    ndof_rho = Twists.dof_rho;
    dcp = dcp/ld;% convert back to the [0,L] domain

    for ii = 1:nip
        if Ws(ii)>0 && u_xi(ii)
            dc_here = dc(:, ii);
            dcp_here = dcp(:, ii); 
            xi_here = xi_star(6*(ii-1)+1:6*ii, 1);
            B_xi_here = Twists.B_xi(6*(ii-1)+1:6*ii, :);
            if ndof_xi>0
                xi_here = B_xi_here*q_xi+xi_here;
            end

            rho_here = rho_star(ii, 1);
            rho_prime_here = 0;
            B_rho_here = Twists.B_rho(ii, :);
            B_rho_prime_here = Twists.B_rho_prime(ii, :);

            if ndof_rho>0
                rho_here = B_rho_here*q_rho+rho_here;
                rho_prime_here = B_rho_prime_here*q_rho + rho_prime_here;
            end

            xihat_here123 = [rho_prime_here -xi_here(3)*rho_here xi_here(2)*rho_here xi_here(4);
                             xi_here(3)*rho_here rho_prime_here -xi_here(1)*rho_here xi_here(5);
                             -xi_here(2)*rho_here xi_here(1)*rho_here rho_prime_here xi_here(6)];
            tang = xihat_here123*[dc_here;1] + rho_here*dcp_here;
            tang = tang/norm(tang);
            Btau = [rho_here*[0 -dc_here(3) dc_here(2);
                     dc_here(3) 0 -dc_here(1);
                     -dc_here(2) dc_here(1) 0]*tang;
                    tang];
            Bq_xi = Bq_xi + ld*Ws(ii)*B_xi_here'*Btau*u_xi(ii);
        end
    end
end
