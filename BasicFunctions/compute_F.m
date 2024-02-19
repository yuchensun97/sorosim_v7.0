function F = compute_F(twist, rho, rho_s, X)
%COMPUTE_F compute the linearized deformation tensor 
% author: Yuchen SUN
%   twist -- 6X1, strain twist
%   rho   -- scalar, inflation ratio
%   rho_s -- scalar, derivate of inflation ratio wrt arclength s
%   X     -- 2X1, material point position vector at s
%% calculation starts here
a = twist(4) + X(1) * rho_s - X(2) * twist(3);
b = twist(5) + X(1) * twist(3) - X(2) * rho_s;
c = twist(6) - X(1) * twist(2) + X(2) * twist(1);

F = [rho, 0, a;
     0, rho, b;
     0, 0, c]

end

