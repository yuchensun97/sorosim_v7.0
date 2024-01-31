% damping matrix for the strain equation
function D_xi = findD_xi(Tr)
    ndof = Tr.ndof_xi;
    D_xi = zeros(ndof_xi, ndof_xi);

    % joint
    dof_joint = Tr.Twists(1).dof_xi;
    D_xi(1:dof_joint, 1:dof_joint) = zeros(dof_joint);

    % soft body
    ndof = ndof - dof_joint;
    ld = Tr.Link.L;
    Gs = Tr.Twists(2).Gs;
    Ws = Tr.Twists(2).Ws;
    nip = Tr.Twists(2).nip;

    Dtemp = zeros(ndof, ndof);
    phi = Tr.Twists(2).B_xi;

    for ii=1:nip
        if Ws(ii)>0
            Gs_here = Gs((ii-1)*6+1:ii*6,:);
            Dtemp = Dtemp + ld*Ws(ii)*phi((ii-1)*6+1:ii*6,:)'*Gs_here*phi((ii-1)*6+1:ii*6,:);
        end
    end
    D_xi(1+dof_joint:end, 1+dof_joint:end) = Dtemp;
end
