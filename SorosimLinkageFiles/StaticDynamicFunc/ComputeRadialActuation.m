function Bq_rho = ComputeRadialActuation(Tr, rc, r_local, q_rho)
    Bq_rho = zeros(Tr.ndof_rho, Tr.n_ract);
    nsig = Tr.nsig;
    n_ract = Tr.n_ract;

    r_fn = Tr.Link.r_fn;
    Xs = Tr.Twist(2).Xs;

    %TODO: compute \Phi_a
    for i_ract = 1:n_ract
        rc_here = rc(i_ract); %X coordinates of the radial actuator
        r_here = r_fn(rc_here); %radius of the current cross-section
        local_here = r_local(i_ract);
        
    end

    %TODO: multiply by rho basis
    %TODO: then integrate over the rod length
end
