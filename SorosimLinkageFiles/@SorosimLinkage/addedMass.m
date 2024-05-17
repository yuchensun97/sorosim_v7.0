% compute the added mass matrix in referenced configuration

function M_added = addedMass(Tr)

    rho_w = Tr.rho_w;
    r_fn = Tr.Link.r_fn;
    nip = Tr.Twists(2).nip;
    Xs = Tr.Twists(2).Xs;

    M_added = zeros(6*nip, 6);
    Bly = 0.6;
    Blz = 0.6;

    for ii=1:nip
        r_here = r_fn(Xs(ii));
        M_added((ii-1)*6+1:ii*6,:) = pi * r_here^2 * rho_w * diag([0 0 0 0 Bly Blz]);
    end
end
