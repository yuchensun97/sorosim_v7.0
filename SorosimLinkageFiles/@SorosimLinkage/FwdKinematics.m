function [g, rho] = FwdKinematics(Tr, q_xi, q_rho)
%FwdKinematics: numericallly map general coordinates to Euclidean space
%   Detailed explanation goes here
%   Tr      : SorosimTwist class
%   q_xi    : general coordinates for strains
%   q_rho   : general coordinates for inflation ratio
%   returns:
%   g       : 4x4 homogeneous transformation matrix
%   rho     : inflation ratio

    % initialization
    nsig = Tr.nsig;
    ndof_xi = Tr.ndof;
    ndof_rho = Tr.ndof_rho;
    g_here = Tr.g_base; % g at base
    rho_here = Tr.rho_base; % rho at base
    J_xi_here = zeros(6, ndof_xi); % Jacobian of twist at base
    J_rho_here = zeros(1, ndof_xi); % Jacobian of rho at base, the base for rho itself is 0
    g = zeros(4*nsig, 4);
    rho = zeros(nsig, 1);
    g(1:4, :) = g_here;
    rho(1) = rho_here;
    Xs = Tr.Twist(2).Xs;

    if Tr.Z_order==4
        B_Z1 = Tr.Twist(2).B_Z1_xi;
        B_Z2 = Tr.Twist(2).B_Z2_xi;
    else    % Z_order==2
        B_Z = Tr.Twist(2).B_Z_xi;
    end

    for ii = 2:nsig
        %TODO: kinematics and Jacobian updates here
    end

end %eof

