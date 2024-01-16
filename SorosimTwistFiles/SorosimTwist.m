classdef SorosimTwist
    %SSorosimTwist class that assigns DoFs and Bases to the soft link and
    %   can calculate the twist and inflation ratio given joint angles
    
    % get of the the multiple links loop and base type
    properties
        B_xi_dof    %(6x1) array specifying the allowable DoFs of a soft piece. 1 if allowed, 0 if not
        B_xi_odr    %(6x1) array specifying the order of allowed DoFs.(0: constant, 1: linear, 2: quadratic, ...)
        dof_xi      %DoFs of the each twist base
        B_rho_dof   %scalr specifying the allowable inflation ratio. 1 if allowed, 0 if not
        B_rho_odr   %scalar specifying the order of allowed inflation ratio.(0: constant, 1: linear, 2: quadratic, ...)
        dof_rho     %DoFs of the each inflation ratio base

        % for twist
        Bh_xi       %Function handler for base
        B_xi        %(6xdof) Base matrix calculated at lumped joints or ((6xnGauss)xdof) base matrices computed at every significant points of a soft division
        B_Z1_xi     %Base calculated at 4th order first Zanna point (Xs+Z1*(delta(Xs)))
        B_Z2_xi     %Base calculated at 4th order second Zanna point (Xs+Z2*(delta(Xs)))
        B_Z_xi      %Base calculated at 2nd order Zanna point
        
        % for inflation ratio
        Bh_rho      %Function handler for base
        B_rho       %(1xdof) Base matrix calculated at lumped joints or ((1xnGauss)xdof) base matrices computed at every significant points of a soft division

        % initial position
        xi_star     %(6x1) Reference strain vector at the initial position
        rho_star    %Reference inflation ratio at the initial position

        Link        %Link associated with this twist only for soft link
        nip         %Number of integration points including boundary points
        Xs          %integration points
        Ws          %integration weights
        Ms          %(6nipx6) Inertia matrix of cross section.
        Es          %(6nipx6) Stiffness matrix.
        Gs          %(6nipx6) Damping matrix.
    end
    
    methods
        % TODO: handle multiple constructors
        function T = SorosimTwist(link, B_xi, B_rho)
            %SorosimTwist Constructor
            %   link: current soft link
            %   B_xi: (6x2) array specifying the allowable DoFs and oder of a soft piece.
            %         1st column: 1 if allowed, 0 if not
            %         2nd column: 0: constant, 1: linear, 2: quadratic, ...
            %   B_rho: (1x2) array specifying the allowable DoFs and oder of a soft piece.
            %         1st column: 1 if allowed, 0 if not
            %         2nd column: 0: constant, 1: linear, 2: quadratic, ...
            %   stars: (6x1) reference strain vector and (1x1) reference inflation ratio
            T.Link = link;

            % Zanna allocation is applied to twist
            Z1 = 1/2-sqrt(3)/6;     %Zanna quadrature coefficients 4th order
            Z2 = 1/2+sqrt(3)/6;     %Zanna quadrature coefficients 4th order
            Z = 1/2;                %Zanna quadrature coefficients 2nd order   
            
            % setup integration points and weights
            nGauss = 10;            %number of Gauss points
            [Xs, Ws, nip] = GaussQuadrature(nGauss);

            % load B_xi_dof, B_rho_dof, B_xi_odr, B_rho_odr 
            T.B_xi_dof = B_xi(:,1);
            T.B_rho_dof = B_rho(:,1);
            T.B_xi_odr = B_xi(:,2);
            T.B_rho_odr = B_rho(:,2);

            % for twist
            dof_xi = sum(B_xi_dof.*(B_xi_odr+1));

            B_xi = zeros(nip*6, dof_xi);
            B_xi_Z1 = zeros(nip*6, dof_xi);
            B_xi_Z2 = zeros(nip*6, dof_xi);
            B_xi_Z = zeros(nip*6, dof_xi);

            X = Xs(1);
            B_xi(1:6, :) = Phi_Xi_LegendrePolynomial(X, B_xi_dof, B_xi_odr);
            for ii=2:nip
                X = Xs(ii);
                B_xi((ii-1)*6+1:ii*6, :) = Phi_Xi_LegendrePolynomial(X, B_xi_dof, B_xi_odr);
                X = Xs(ii-1)+Z1*(Xs(ii)-Xs(ii-1));
                B_xi_Z1((ii-2)*6+1:ii*6, :) = Phi_Xi_LegendrePolynomial(X, B_xi_dof, B_xi_odr);
                X = Xs(ii-1)+Z2*(Xs(ii)-Xs(ii-1));
                B_xi_Z2((ii-2)*6+1:ii*6, :) = Phi_Xi_LegendrePolynomial(X, B_xi_dof, B_xi_odr);
                X = Xs(ii-1)+Z*(Xs(ii)-Xs(ii-1));
                B_xi_Z((ii-2)*6+1:ii*6, :) = Phi_Xi_LegendrePolynomial(X, B_xi_dof, B_xi_odr);
            end
            file = 'Phi_Xi_LegendrePolynomial';
            Bh_xi = str2func(['@(X, B_xi_dof, B_xi_odr)', file, '(X, B_xi_dof, B_xi_odr)']);	

            % for inflation ratio
            dof_rho = sum(B_rho_dof.*(B_rho_odr+1));

            B_rho = zeros(nip, dof_rho);
            X = Xs(1);
            B_rho(1, :) = Phi_Rho_LegendrePolynomial(X, B_rho_dof, B_rho_odr);
            for ii=2:nip
                X = Xs(ii);
                B_rho(ii, :) = Phi_Rho_LegendrePolynomial(X, B_rho_dof, B_rho_odr);
            end
            file = 'Phi_Rho_LegendrePolynomial';
            Bh_rho = str2func(['@(X, B_rho_dof, B_rho_odr)', file, '(X, B_rho_dof, B_rho_odr)']);

            % initial position, simpilified to the undeformed position
            xi_star = zeros(6*nip, 4); % precomputation at all gauss and zannah guess points
            xi_star(6:6:end, :) = ones(4, 4);
            rho_star = ones(nip, 1);

            % TODO: precompute Ms, Es, Gs

            T.xi_star = xi_star;
            T.rho_star = rho_star;
            T.nip = nip;
            T.Xs = Xs;
            T.Ws = Ws;
            T.B_xi = B_xi;
            T.B_Z1_xi = B_xi_Z1;
            T.B_Z2_xi = B_xi_Z2;
            T.B_Z_xi = B_xi_Z;
            T.Bh_xi = Bh_xi;
            T.B_rho = B_rho;
            T.Bh_rho = Bh_rho;
            T.dof_xi = dof_xi;
            T.dof_rho = dof_rho;
        end

    end
end

