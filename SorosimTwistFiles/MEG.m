%This function computes Ms, Es, and Gs for a given Xs for the jth division
%of a link. Modified from Anup Mathew (15.12.2022)

function [Ms,Es,Gs]= MEG(Link, Xs)
    np = length(Xs);
    r_fn     = Link.r_fn;
    %updating:
    r_nGauss = zeros(np,1);
    for ii=1:np
        r_nGauss(ii) = r_fn(Xs(ii));
    end
    Iy_p = (pi/4)*r_nGauss.^4;
    Iz_p = Iy_p;
    Ix_p = Iy_p+Iz_p;
    A_p  = pi*r_nGauss.^2;
    
    Ms = zeros(6*np,6); %inertia
    Es = zeros(6*np,6); %stiffness
    Gs = zeros(6*np,6); %damping
    
    Rho0  = Link.Rho0;
    G    = Link.G;
    E    = Link.E;
    Eta  = Link.Eta;
    for ii=1:np
        Ms((ii-1)*6+1:ii*6,:) = Rho0*diag([Ix_p(ii),Iy_p(ii),Iz_p(ii),A_p(ii),A_p(ii),A_p(ii)]);
        Es((ii-1)*6+1:ii*6,:) = diag([G*Ix_p(ii),E*Iy_p(ii),E*Iz_p(ii),E*A_p(ii),G*A_p(ii),G*A_p(ii)]);
        Gs((ii-1)*6+1:ii*6,:) = Eta*diag([Ix_p(ii),3*Iy_p(ii),3*Iz_p(ii),3*A_p(ii),A_p(ii),A_p(ii)]);
    end
    
end %eof
