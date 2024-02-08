% The first two terms of K_rho in the second equation
function K_rho_part = findK_rho_part(Tr)
    ndof = Tr.ndof_rho;
    K_rho_part = zeros(ndof, ndof);

    % joint
    dof_joint = Tr.Twists(1).dof_rho;
    K_rho_part(1:dof_joint, 1:dof_joint) = zeros(dof_joint);

    % soft body
    ndof = ndof - dof_joint;
    ld = Tr.Link.L;
    Es = Tr.Twists(2).Es;
    Ws = Tr.Twists(2).Ws;
    sigma_rho = Tr.Twists(2).sigma_rho;
    nip = Tr.Twists(2).nip;

    Ktemp = zeros(ndof, ndof);
    phi = Tr.Twists(2).B_rho;
    phi_p = Tr.Twists(2).B_rho_prime;

    for ii=1:nip
        if Ws(ii)>0
            Es_here = Es((ii-1)*6+1:ii*6,:);
            sigma_here = Es_here(1,1);
            sigma_rho_here = sigma_rho(ii);
            Ktemp = Ktemp+ld*Ws(ii)*phi_p(ii,:)'*sigma_here*phi_p(ii,:)-...
                    ld*Ws(ii)*phi(ii,:)'*sigma_rho_here*phi(ii,:);
        end
    end
    K_rho_part(1+dof_joint:end, 1+dof_joint:end) = Ktemp;
end
