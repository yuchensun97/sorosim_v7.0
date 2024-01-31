% elastic matrix in strain equations
function K_xi = findK_xi(Tr)
    ndof = Tr.ndof_xi;
    K_xi = zeros(ndof_xi, ndof_xi);

    % joint
    dof_joint = Tr.Twists(1).dof_xi;
    K_xi(1:dof_joint, 1:dof_joint) = zeros(dof_joint);

    % soft body
    ndof = ndof - dof_joint;
    ld = Tr.Link.L;
    Es = Tr.Twists(2).Es;
    Ws = Tr.Twists(2).Ws;
    nip = Tr.Twists(2).nip;

    Ktemp = zeros(ndof, ndof);
    phi = Tr.Twists(2).B_xi;

    for ii=1:nip
        if Ws(ii)>0
            Es_here = Es((ii-1)*6+1:ii*6,:);
            Ktemp = Ktemp + ld*Ws(ii)*phi((ii-1)*6+1:ii*6,:)'*Es_here*phi((ii-1)*6+1:ii*6,:);
        end
    end
    K_xi(1+dof_joint:end, 1+dof_joint:end) = Ktemp;
end
