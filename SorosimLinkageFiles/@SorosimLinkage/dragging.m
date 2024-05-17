% compute the dragging/lifting matrix in referenced configuration

function DL = dragging(Tr)
    Cdx = 0;
    Cdy = 1.1;
    Cdz = Cdy;
    Clx = -0.1;
    Cly = 0;
    Clz = 0;
    rho_w = Tr.rho_w;

    r_fn = Tr.Link.r_fn;
    nip = Tr.Twists(2).nip;
    Xs = Tr.Twists(2).Xs;

    DL = zeros(6*nip, 6);

    for ii=1:nip
        r_here = r_fn(Xs(ii));
        DL((ii-1)*6+1:ii*6, :) = r_here * rho_w * [zeros(3, 3) zeros(3, 3);
                                                   zeros(3, 3) [0.5*pi*Cdx -Clz Cly; 
                                                                Clz Cdy -Clx; 
                                                                -Cly Clx Cdz]];                                            
    end
end
