% damping matrix for the lateral term in the strain equation
%elastic matrix for the \rho term in strain equation
function D_xi_bar = findD_xi_bar(Tr)
    ndof_xi = Tr.ndof_xi;
    ndof_rho = Tr.ndof_rho;
    D_xi_bar = zeros(ndof_xi, ndof_rho);

    % joint
    dof_joint_xi = Tr.Twists(1).dof_xi;
    dof_joint_rho = Tr.Twists(1).dof_rho;
    D_xi_bar(1:dof_joint_xi, 1:dof_joint_rho) = zeros(dof_joint_xi, dof_joint_rho);

    % soft body
    ndof_xi = ndof_xi - dof_joint_xi;
    ndof_rho = ndof_rho - dof_joint_rho;
    ld = Tr.Link.L;
    gamma = Tr.Twists(2).gamma_xi;
    Ws = Tr.Twists(2).Ws;
    nip = Tr.Twists(2).nip;

    Dtemp = zeros(ndof_xi, ndof_rho);
    phi_xi = Tr.Twists(2).B_xi;
    phi_rho = Tr.Twists(2).B_rho;

    for ii=1:nip
        if Ws(ii)>0
            gamma_here = gamma(ii);
            psi = zeros(6, ndof_rho);
            psi(4, :) = phi_rho(ii,:);
            Dtemp = Dtemp + ld*Ws(ii)*phi_xi((ii-1)*6+1:ii*6,:)'*gamma_here*psi;
        end
    end
    D_xi_bar(1+dof_joint_xi:end, 1+dof_joint_rho:end) = Dtemp;
end
