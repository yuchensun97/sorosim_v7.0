function Phi_Prime_rho = Phi_Prime_Rho_LegendrePolynomial(X, Bdof, Bodr)
    %PHI_Prime_Rho_LegendrePolynomial generate the derivatives of basis function
    %wrt arclength of Legendre
    %   X -- varites from 0 to 1
    %   Bdof -- tells if the inflation ratio is on or off
    %   eg. Bdof = 1 implies that the inflation ratio is on
    %   Bodr -- defines the order of inflation ratio
    %   eg. Bodr = 1 implies linear inflation ratio.
    
    dof = sum(Bdof.*(Bodr+1));
    Phi_rho = zeros(1, dof);
    Phi_Prime_rho = zeros(1, dof);
    
    X = 2*X - 1; % transform from [0, 1] to [-1, 1]	
    k = 1;
    P0 = 1;
    P1 = X;
    P0_prime = 0;
    P1_prime = 1;
    for j = 1:dof
        if j==1
            Phi_rho(k) = P0;
            Phi_Prime_rho(k) = P0_prime;
        end
        if j==2
            Phi_rho(k) = P1;
            Phi_Prime_rho(k) = P1_prime;
        end
        if j>2
            n = j - 2;
            P1t = P1;
            P1 = ((2*n+1)*X*P1 - n*P0)/(n+1);
            P0 = P1t;
            m = k - 1; % the current order of the polynomial
            P1_prime = m * Phi_rho(k-1) + X * Phi_Prime_rho(k-1);
            % recurrence update
            Phi_rho(k) = P1;
            Phi_Prime_rho(k) = P1_prime;
        end
        k = k + 1;
    end
end
    