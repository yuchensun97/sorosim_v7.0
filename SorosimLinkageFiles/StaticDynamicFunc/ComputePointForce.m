function Fp = ComputePointForce(Tr, J_xi, g, t)
    np = Tr.np;
    Fp = zeros(Tr.ndof_xi, 1);
    Xs = Tr.Twists(2).Xs;
    nsig = Tr.nsig;

    for ip = 1:np
        Fp_loc = Tr.Fp_loc(ip);
        Fp_vec = Tr.Fp_vec(ip)(t);
        LocalForce = Tr.LocalForce(ip);
        [rows, cols] = size(Fp_vec);
        if rows ~= 6 || cols ~= 1
            error('Input Point Force must be 6x1')
        end
        
        for ii=1:nsig
            if Xs(ii)==Fp_loc
                if ~LocalForce
                    g_here = g((ii-1)*4+1:ii*4, :);
                    g_here(1:3, 4) = zeros(3,1);
                    Ad_g_here_inv = dinamico_Adjoint(ginv(g_here));
                    Fp_vec = Ad_g_here_inv*Fp_vec;
                end
                Fp=Fp+J_xi((ii-1)*4+1:ii*4,:);
                break;
            end
        end
    end
end
