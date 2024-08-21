function Phi_Prime_Rho = Phi_Prime_Rho_Hermitian_robin(X, Bdof, Bodr)
    %Phi_Prime_Rho_Hermitian generate the basis function of robin boundary condition
%   X -- varites from 0 to 1
%   Bdof -- tells if the inflation ratio is on or off
%   eg. Bdof = 1 implies that the inflation ratio is on
%   Bodr -- defines the order of inflation ratio
%   eg. Bodr = 1 implies linear inflation ratio.
    nele = Bodr; %gives the number of elements
    dof = Bdof*(2*nele+2);
    Phi_Prime_Rho = zeros(1,dof);

    if nele == 0
        error("Error: number of elements must be greater than 0")
    end

    if dof == 0
        return
    end

    k = 1;
    w = 1/nele;
    a = 0;
    for j=1:nele
        b = a + w;
        if X>=a && X<=b
            Xc = (X-a)/(b-a);
            dXcdX = 1/(b-a);
            Phi_Prime_Rho(k)   = (-6*Xc+6*Xc^2)*dXcdX;
            Phi_Prime_Rho(k+1) = (1-4*Xc+3*Xc^2)*dXcdX;
            Phi_Prime_Rho(k+2) = (6*Xc-6*Xc^2)*dXcdX;
            Phi_Prime_Rho(k+3) = (-2*Xc+3*Xc^2)*dXcdX;
        end
        k=k+2;
        a=a+w;
    end
end
