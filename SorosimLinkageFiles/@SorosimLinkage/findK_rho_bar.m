% elastic matrix for the strain term in the lateral equation
function K_rho_bar = findK_rho_bar(Tr)
    ndof_xi = Tr.ndof_xi;
    ndof_rho = Tr.ndof_rho;
    K_rho_bar = zeros(ndof_rho, ndof_xi);

    % joint
    dof_joint_xi = Tr.Twists(1).dof_xi;
    dof_joint_rho = Tr.Twists(1).dof_rho;
    K_rho_bar(1:dof_joint_rho, 1:dof_joint_xi) = zeros(dof_joint_rho, dof_joint_xi);

    % soft body
    ndof_xi = ndof_xi - dof_joint_xi;
    ndof_rho = ndof_rho - dof_joint_rho;
    ld = Tr.Link.L;
    sigma = Tr.Twists(2).sigma_xi;
    Ws = Tr.Twists(2).Ws;
    nip = Tr.Twists(2).nip;

    Ktemp = zeros(ndof_rho, ndof_xi);
    phi_xi = Tr.Twists(2).B_xi;
    phi_rho = Tr.Twists(2).B_rho;

    for ii=1:nip
        if Ws(ii)>0
            sigma_here = sigma(ii);
            phi_xi_here = phi_xi((ii-1)*6+1:ii*6,:);
            psi = phi_xi_here(4,:);
            Ktemp = Ktemp + ld*Ws(ii)*phi_rho(ii,:)'*sigma_here*psi;
        end
    end
    K_rho_bar(1+dof_joint_xi:end, 1+dof_joint_rho:end) = Ktemp;
end
