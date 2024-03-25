% damping matrix for the lateral equation
function D_rho = findD_rho(Tr)
    ndof = Tr.ndof_rho;
    D_rho = zeros(ndof, ndof);

    % joint
    dof_joint = Tr.Twists(1).dof_rho;
    D_rho(1:dof_joint, 1:dof_joint) = zeros(dof_joint);

    % soft body
    ndof = ndof - dof_joint;
    ld = Tr.Link.L;
    Gs = Tr.Twists(2).Gs;
    Ws = Tr.Twists(2).Ws;
    gamma = Tr.Twists(2).gamma_rho;
    nip = Tr.Twists(2).nip;

    Dtemp = zeros(ndof, ndof);
    phi = Tr.Twists(2).B_rho;
    phi_p = Tr.Twists(2).B_rho_prime;

    for ii=1:nip
        if Ws(ii)>0
            Gs_here = Gs((ii-1)*6+1:ii*6,:);
            sigma_here = Gs_here(1,1);
            gamma_here = gamma(ii);
            Dtemp = Dtemp+ld*Ws(ii)*phi_p(ii,:)'*sigma_here*phi_p(ii,:)+...
                    ld*Ws(ii)*phi(ii,:)'*gamma_here*phi(ii,:);
        end
    end
    D_rho(1+dof_joint:end, 1+dof_joint:end) = Dtemp;
end
