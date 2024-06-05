%This function computes Ms, Es, and Gs for a given Xs for the jth division
%of a link. Modified from Anup Mathew (15.12.2022)

function [sigma,gamma]= EG_rho_bar(Link, Xs)
    np = length(Xs);
    r_fn     = Link.r_fn;
    %updating:
    r_nGauss = zeros(np,1);
    for ii=1:np
        r_nGauss(ii) = r_fn(Xs(ii));
    end
    A_p  = pi*r_nGauss.^2;
    
    sigma = zeros(np,1); %stiffness
    gamma = zeros(np,1); %damping
    
    G    = Link.G; % shear
    Lam = Link.Lamda; 
    Eta  = Link.Eta;
    for ii=1:np
        sigma(ii) = 4*(G + Lam)*A_p(ii);
        gamma(ii) = 8*Eta*A_p(ii);
    end
    
end %eof
