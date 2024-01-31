function J_rho_prime = Jacobianprime(Tr)

    ndof_rho = Tr.ndof_rho;
    nsig = Tr.nsig;
    J_rho_prime = zeros(nsig, ndof_rho);

    % joint
    dof_rho_joint = Tr.Twists(1).dof_rho;
    B_rho_prime = Tr.Twists(1).B_rho_prime;
    
    if dof_rho_joint > 0
        J_rho_prime(:, 1:dof_rho_joint) = B_rho_prime;
    end
    J_rho_prime(:,dof_rho_joint+1:end) = Tr.Twists(2).B_rho_prime;
end
