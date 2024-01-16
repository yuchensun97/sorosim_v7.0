function Phi_rho = Phi_Rho_LegendrePolynomial(X, Bdof, Bodr)
%PHI_Rho_LegendrePolynomial generate the basis function of Legendre
%   X -- varites from 0 to 1
%   Bdof -- tells if the inflation ratio is on or off
%   eg. Bdof = 1 implies that the inflation ratio is on
%   Bodr -- defines the order of inflation ratio
%   eg. Bodr = 1 implies linear inflation ratio.

dof = sum(Bdof.*(Bodr+1));
Phi_rho = zeros(1, dof);

X = 2*X - 1; % transform from [0, 1] to [-1, 1]	
k = 1;
P0 = 1;
P1 = X;
for j = 1:dof
    if j==1
        Phi_rho(k) = P0;
    end
    if j==2
        Phi_rho(k) = P1;
    end
    if j>2
        n = j - 2;
        P1t = P1;
        P1 = ((2*n+1)*X*P1 - n*P0)/(n+1);
        P0 = P1t;
        Phi_rho(k) = P1;
    end
    k = k + 1;
end

end
