% elastic matrix in strain equations
function M_rho = findM_rho(Tr)
    ndof = Tr.ndof_rho;
    M_rho = zeros(ndof, ndof);

    % joint
    dof_joint = Tr.Twists(1).dof_xi;
    M_rho(1:dof_joint, 1:dof_joint) = zeros(dof_joint);

    % soft body
    ndof = ndof - dof_joint;
    ld = Tr.Link.L;
    Ms = Tr.Twists(2).Ms;
    Ws = Tr.Twists(2).Ws;
    nip = Tr.Twists(2).nip;

    Mtemp = zeros(ndof, ndof);
    phi = Tr.Twists(2).B_rho;

    for ii=1:nip
        if Ws(ii)>0
            Ms_here = Ms((ii-1)*6+1:ii*6,:);
            II_here = Ms_here(2, 2) + Ms_here(3, 3);
            Mtemp = Mtemp + ld*Ws(ii)*phi(ii,:)'*II_here*phi(ii,:);
        end
    end
    M_rho(1+dof_joint:end, 1+dof_joint:end) = Mtemp;
end
