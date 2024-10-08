classdef SorosimTwist < handle
    %SSorosimTwist class that assigns DoFs and Bases to the soft link and
    %   can calculate the twist and inflation ratio given joint angles
    
    % get of the the multiple links loop and base type
    properties
        basisType   % Base type for rho(dirichlet, neumann, robin, mixed)
        B_xi_dof    %(6x1) array specifying the allowable DoFs of a soft piece. 1 if allowed, 0 if not
        B_xi_odr    %(6x1) array specifying the order of allowed DoFs.(0: constant, 1: linear, 2: quadratic, ...)
        dof_xi      %DoFs of the each twist base
        B_rho_dof   %scalr specifying the allowable inflation ratio. 1 if allowed, 0 if not
        B_rho_odr   %scalar specifying the order of allowed inflation ratio.(0: constant, 1: linear, 2: quadratic, ...)
        dof_rho     %DoFs of the each inflation ratio base

        % for twist
        Bh_xi       %Function handler for base
        B_xi        %(6nxdof) Base matrix calculated at lumped joints or ((6xnGauss)xdof) base matrices computed at every significant points of a soft division
        B_Z1_xi     %Base calculated at 4th order first Zanna point (Xs+Z1*(delta(Xs)))
        B_Z2_xi     %Base calculated at 4th order second Zanna point (Xs+Z2*(delta(Xs)))
        B_Z_xi      %Base calculated at 2nd order Zanna point
        
        % for inflation ratio
        Bh_rho      %Function handler for base
        Bh_rho_prime%Function handler for base
        B_rho       %(nxdof) Base matrix calculated at lumped joints or ((1xnGauss)xdof) base matrices computed at every significant points of a soft division
        B_rho_prime %(nxdof) Derivative of B_rho

        % initial position
        xi_star     %(6nx1) Reference strain vector at the initial position
        rho_star    %Reference inflation ratio at the initial position
        xi_starfn   %Function handler for xi_star
        rho_starfn  %Function handler for rho_star

        Link        %Link associated with this twist only for soft link
        nip         %Number of integration points including boundary points
        nGauss      %Number of Gauss points
        Xs          %integration points
        Ws          %integration weights
        Ms          %(6nipx6) Inertia matrix of cross section.
        Es          %(6nipx6) Stiffness matrix.
        Gs          %(6nipx6) Damping matrix.
        sigma_xi    %(nipx1) stiffness term for rho in strain equation, or for nu in lateral equation
        gamma_xi    %(nipx1) damping term for rho in strain equation, or for nu in lateral equation
        sigma_rho   %(nipx1) stiffness term for rho in lateral equation
        gamma_rho   %(nipx1) damping term for rho in lateral equation

        Xadd        %additional integration points(nx1)
    end
    
    methods
        % TODO: handle multiple constructors
        function T = SorosimTwist(varargin)
            %SorosimTwist Constructor
            if nargin == 4
                link = varargin{1};
                B_xi_in = varargin{2};
                B_rho_in = varargin{3};
                basisType = varargin{4};
                T.basisType = basisType;
                %   link: current soft link
                %   B_xi_in: (6x2) array specifying the allowable DoFs and oder of a soft piece.
                %         1st column: 1 if allowed, 0 if not
                %         2nd column: 0: constant, 1: linear, 2: quadratic, ...
                %   B_rho_in: (1x2) array specifying the allowable DoFs and oder of a soft piece.
                %         1st column: 1 if allowed, 0 if not
                %         2nd column: 0: constant, 1: linear, 2: quadratic, ...
                %   stars: (6x1) reference strain vector and (1x1) reference inflation ratio
                T.Link = link;

                % Zanna allocation is applied to twist
                Z1 = 1/2-sqrt(3)/6;     %Zanna quadrature coefficients 4th order
                Z2 = 1/2+sqrt(3)/6;     %Zanna quadrature coefficients 4th order
                Z = 1/2;                %Zanna quadrature coefficients 2nd order   
                
                % setup integration points and weights
                nGauss = 20;            %number of Gauss points
                T.nGauss = nGauss;
                [Xs, Ws, nip] = GaussQuadrature(nGauss);

                % load B_xi_dof, B_rho_dof, B_xi_odr, B_rho_odr 
                B_xi_dof = B_xi_in(:, 1);
                B_rho_dof = B_rho_in(1,1);
                B_xi_odr = B_xi_in(:, 2);
                B_rho_odr = B_rho_in(1, 2);
                T.B_xi_dof = B_xi_in(:,1);
                T.B_rho_dof = B_rho_in(1,1);
                T.B_xi_odr = B_xi_in(:,2);
                T.B_rho_odr = B_rho_in(1,2);

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
                    B_xi_Z1((ii-2)*6+1:(ii-1)*6, :) = Phi_Xi_LegendrePolynomial(X, B_xi_dof, B_xi_odr);
                    X = Xs(ii-1)+Z2*(Xs(ii)-Xs(ii-1));
                    B_xi_Z2((ii-2)*6+1:(ii-1)*6, :) = Phi_Xi_LegendrePolynomial(X, B_xi_dof, B_xi_odr);
                    X = Xs(ii-1)+Z*(Xs(ii)-Xs(ii-1));
                    B_xi_Z((ii-2)*6+1:(ii-1)*6, :) = Phi_Xi_LegendrePolynomial(X, B_xi_dof, B_xi_odr);
                end
                file = 'Phi_Xi_LegendrePolynomial';
                Bh_xi = str2func(['@(X, Bdof, Bodr)', file, '(X, Bdof, Bodr)']);	

                % for inflation ratio
                switch basisType
                    case 'legendre'
                        dof_rho = sum(B_rho_dof.*(B_rho_odr+1));
                        B_rho = zeros(nip, dof_rho);
                        X = Xs(1);
                        B_rho(1, :) = Phi_Rho_LegendrePolynomial(X, B_rho_dof, B_rho_odr);
                        B_rho_prime(1, :) = Phi_Prime_Rho_LegendrePolynomial(X, B_rho_dof, B_rho_odr);
                        for ii=2:nip
                            X = Xs(ii);
                            B_rho(ii, :) = Phi_Rho_LegendrePolynomial(X, B_rho_dof, B_rho_odr);
                            B_rho_prime(ii, :) = Phi_Prime_Rho_LegendrePolynomial(X, B_rho_dof, B_rho_odr);
                        end
                        file = 'Phi_Rho_LegendrePolynomial';
                        file_prime = 'Phi_Prime_Rho_LegendrePolynomial';
                        Bh_rho = str2func(['@(X, Bdof, Bodr)', file, '(X, Bdof, Bodr)']);
                        Bh_rho_prime = str2func(['@(X, Bdof, Bodr)', file_prime, '(X, Bdof, Bodr)']);
                    case 'mixed'
                        dof_rho = sum(B_rho_dof.*(2*B_rho_odr));
                        B_rho = zeros(nip, dof_rho);
                        X = Xs(1);
                        B_rho(1, :) = Phi_Rho_Hermitian(X, B_rho_dof, B_rho_odr);
                        B_rho_prime(1, :) = Phi_Prime_Rho_Hermitian(X, B_rho_dof, B_rho_odr);
                        for ii=2:nip
                            X = Xs(ii);
                            B_rho(ii, :) = Phi_Rho_Hermitian(X, B_rho_dof, B_rho_odr);
                            B_rho_prime(ii, :) = Phi_Prime_Rho_Hermitian(X, B_rho_dof, B_rho_odr);
                        end
                        file = 'Phi_Rho_Hermitian';
                        file_prime = 'Phi_Prime_Rho_Hermitian';
                        Bh_rho = str2func(['@(X, Bdof, Bodr)', file, '(X, Bdof, Bodr)']);
                        Bh_rho_prime = str2func(['@(X, Bdof, Bodr)', file_prime, '(X, Bdof, Bodr)']);
                    case 'dirichlet'
                        dof_rho = sum(B_rho_dof.*(2*B_rho_odr));
                        B_rho = zeros(nip, dof_rho);
                        X = Xs(1);
                        B_rho(1, :) = Phi_Rho_Hermitian_dirichlet(X, B_rho_dof, B_rho_odr);
                        B_rho_prime(1, :) = Phi_Prime_Rho_Hermitian_dirichlet(X, B_rho_dof, B_rho_odr);
                        for ii=2:nip
                            X = Xs(ii);
                            B_rho(ii, :) = Phi_Rho_Hermitian_dirichlet(X, B_rho_dof, B_rho_odr);
                            B_rho_prime(ii, :) = Phi_Prime_Rho_Hermitian_dirichlet(X, B_rho_dof, B_rho_odr);
                        end
                        file = 'Phi_Rho_Hermitian_dirichlet';
                        file_prime = 'Phi_Prime_Rho_Hermitian_dirichlet';
                        Bh_rho = str2func(['@(X, Bdof, Bodr)', file, '(X, Bdof, Bodr)']);
                        Bh_rho_prime = str2func(['@(X, Bdof, Bodr)', file_prime, '(X, Bdof, Bodr)']);
                    case 'robin'
                        dof_rho = sum(B_rho_dof.*(2*B_rho_odr+2));
                        B_rho = zeros(nip, dof_rho);
                        X = Xs(1);
                        B_rho(1, :) = Phi_Rho_Hermitian_robin(X, B_rho_dof, B_rho_odr);
                        B_rho_prime(1, :) = Phi_Prime_Rho_Hermitian_robin(X, B_rho_dof, B_rho_odr);
                        for ii=2:nip
                            X = Xs(ii);
                            B_rho(ii, :) = Phi_Rho_Hermitian_robin(X, B_rho_dof, B_rho_odr);
                            B_rho_prime(ii, :) = Phi_Prime_Rho_Hermitian_robin(X, B_rho_dof, B_rho_odr);
                        end
                        file = 'Phi_Rho_Hermitian_robin';
                        file_prime = 'Phi_Prime_Rho_Hermitian_robin';
                        Bh_rho = str2func(['@(X, Bdof, Bodr)', file, '(X, Bdof, Bodr)']);
                        Bh_rho_prime = str2func(['@(X, Bdof, Bodr)', file_prime, '(X, Bdof, Bodr)']);
                    case 'neumann'
                        dof_rho = sum(B_rho_dof.*(2*B_rho_odr));
                        B_rho = zeros(nip, dof_rho);
                        X = Xs(1);
                        B_rho(1, :) = Phi_Rho_Hermitian_neumann(X, B_rho_dof, B_rho_odr);
                        B_rho_prime(1, :) = Phi_Prime_Rho_Hermitian_neumann(X, B_rho_dof, B_rho_odr);
                        for ii=2:nip
                            X = Xs(ii);
                            B_rho(ii, :) = Phi_Rho_Hermitian_neumann(X, B_rho_dof, B_rho_odr);
                            B_rho_prime(ii, :) = Phi_Prime_Rho_Hermitian_neumann(X, B_rho_dof, B_rho_odr);
                        end
                        file = 'Phi_Rho_Hermitian_neumann';
                        file_prime = 'Phi_Prime_Rho_Hermitian_neumann';
                        Bh_rho = str2func(['@(X, Bdof, Bodr)', file, '(X, Bdof, Bodr)']);
                        Bh_rho_prime = str2func(['@(X, Bdof, Bodr)', file_prime, '(X, Bdof, Bodr)']);
                    otherwise
                        error("Invalid basis type. Must be 'legendre','mixed', 'dirichlet', 'robin', or 'neumann'.");
                end

                % initial position, simpilified to the undeformed position
                xi_star = zeros(6*nip, 4); % precomputation at all gauss and zannah guess points
                xi_star(4:6:end, :) = ones(nip, 4);
                rho_star = ones(nip, 1);

                xi_starfn = @(X)[0 0 0 1 0 0]';
                rho_starfn = @(X)1;

                % TODO: precompute Ms, Es, Gs
                [Ms, Es, Gs] = MEG(link, Xs);
                [sigma_xi, gamma_xi] = EG_xi_bar(link, Xs);
                [sigma_rho, gamma_rho] = EG_rho_bar(link, Xs);

                T.xi_starfn = xi_starfn;
                T.rho_starfn = rho_starfn;
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
                T.B_rho_prime = B_rho_prime;
                T.Bh_rho = Bh_rho;
                T.Bh_rho_prime = Bh_rho_prime;
                T.dof_xi = dof_xi;
                T.dof_rho = dof_rho;
                T.Ms = Ms;
                T.Es = Es;
                T.Gs = Gs;
                T.sigma_xi = sigma_xi;
                T.gamma_xi = gamma_xi;
                T.sigma_rho = sigma_rho;
                T.gamma_rho = gamma_rho;
            elseif nargin == 2
                T.B_xi = varargin{1};
                T.B_rho = varargin{2};
                T.B_rho_prime = [];
            elseif nargin == 0
                nGauss = 10;
                [Xs, Ws, nip] = GaussQuadrature(nGauss);
                T.nip = nip;
                T.Xs = Xs;
                T.Ws = Ws;
                xi_starfn = @(X)[0 0 0 1 0 0]';
                rho_starfn = @(X)1;
                T.xi_starfn = xi_starfn;
                T.rho_starfn = rho_starfn;
                T.B_xi = [];
                T.B_rho = [];
                T.B_rho_prime = [];
                T.B_xi_dof = zeros(6, 1);
                T.B_rho_dof = 0;
                T.B_xi_odr = zeros(6, 1);
                T.B_rho_odr = 0;
                T.dof_xi = 0;
                T.dof_rho = 0;
            end

        end
    end

    methods
        %% set and property update fuctions
        function T = UpdateBh(T)
            if isempty(T.dof_xi) || isempty(T.dof_rho)
                return
            end
            file_xi = 'Phi_Xi_LegendrePolynomial';
            T.Bh_xi = str2func(['@(X, Bdof, Bodr)', file_xi, '(X, Bdof, Bodr)']);
            switch T.basisType
                case 'legendre'
                    file_rho = 'Phi_Rho_LegendrePolynomial';
                    file_prime = 'Phi_Prime_Rho_LegendrePolynomial';
                    T.Bh_rho = str2func(['@(X, Bdof, Bodr)', file_rho, '(X, Bdof, Bodr)']);
                    T.Bh_rho_prime = str2func(['@(X, Bdof, Bodr)', file_prime, '(X, Bdof, Bodr)']);
                case 'dirichlet'
                    file_rho = 'Phi_Rho_Hermitian_dirichlet';
                    file_prime = 'Phi_Prime_Rho_Hermitian_dirichlet';
                    T.Bh_rho = str2func(['@(X, Bdof, Bodr)', file_rho, '(X, Bdof, Bodr)']);
                    T.Bh_rho_prime = str2func(['@(X, Bdof, Bodr)', file_prime, '(X, Bdof, Bodr)']);
                case 'neumann'
                    file_rho = 'Phi_Rho_Hermitian_neumann';
                    file_prime = 'Phi_Prime_Rho_Hermitian_neumann';
                    T.Bh_rho = str2func(['@(X, Bdof, Bodr)', file_rho, '(X, Bdof, Bodr)']);
                    T.Bh_rho_prime = str2func(['@(X, Bdof, Bodr)', file_prime, '(X, Bdof, Bodr)']);
                case 'robin'
                    file_rho = 'Phi_Rho_Hermitian_robin';
                    file_prime = 'Phi_Prime_Rho_Hermitian_robin';
                    T.Bh_rho = str2func(['@(X, Bdof, Bodr)', file_rho, '(X, Bdof, Bodr)']);
                    T.Bh_rho_prime = str2func(['@(X, Bdof, Bodr)', file_prime, '(X, Bdof, Bodr)']);
                case 'mixed'
                    file_rho = 'Phi_Rho_Hermitian';
                    file_prime = 'Phi_Prime_Rho_Hermitian';
                    T.Bh_rho = str2func(['@(X, Bdof, Bodr)', file_rho, '(X, Bdof, Bodr)']);
                    T.Bh_rho_prime = str2func(['@(X, Bdof, Bodr)', file_prime, '(X, Bdof, Bodr)']);
                otherwise
                    error("Invalid basis type. Must be 'legendre','mixed', 'dirichlet', 'robin', or 'neumann'.")
            end
        end

        %% updates np, Xs, Ws, and dof
        function T = UpdateIntegration(T)
            if isempty(T.dof_xi) || isempty(T.dof_rho)
                return
            end
            if T.nGauss < 5
                T.nGauss = 5;
            end
            [T.Xs, T.Ws, T.nip] = GaussQuadrature(T.nGauss);
        end

        %% updates dof
        function T = Updatedof(T)
            if isempty(T.dof_xi) || isempty(T.dof_rho)
                return
            end
            T.dof_xi = sum(T.B_xi_dof.*(T.B_xi_odr+1));
            switch T.basisType
                case 'legendre'
                    T.dof_rho = sum(T.B_rho_dof.*(T.B_rho_odr+1));
                case 'mixed'
                    T.dof_rho = sum(T.B_rho_dof.*(2*T.B_rho_dof));
                case 'dirichlet'
                    T.dof_rho = sum(T.B_rho_dof.*(2*T.B_rho_dof));
                case 'robin'
                    T.dof_rho = sum(T.B_rho_dof.*(2*T.B_rho_dof+2));
                case 'neumann'
                    T.dof_rho = sum(T.B_rho_dof.*(2*T.B_rho_dof));
                otherwise
                    error("Invalid basis type. Must be 'legendre','mixed', 'dirichlet', 'robin', or 'neumann'.");
            end
        end

        %% updates all the bases properties
        function UpdateBs(T)
            if isempty(T.dof_xi) || isempty(T.dof_rho)
                return
            end

            Z1     = 1/2-sqrt(3)/6;          %Zanna quadrature coefficient
            Z2     = 1/2+sqrt(3)/6;          %Zanna quadrature coefficient
            Z      = 1/2;                    %Zanna quadrature coefficient

            T.B_xi    = zeros(T.nip*6, T.dof_xi);
            T.B_Z1_xi = zeros(T.nip*6, T.dof_xi);
            T.B_Z2_xi = zeros(T.nip*6, T.dof_xi);
            T.B_Z_xi  = zeros(T.nip*6, T.dof_xi);
            T.B_rho   = zeros(T.nip, T.dof_rho);

            ii = 1;
            X = T.Xs(ii);
            T.B_xi(1:6, :) = T.Bh_xi(X, T.B_xi_dof, T.B_xi_odr);
            T.B_rho(ii, :) = T.Bh_rho(X, T.B_rho_dof, T.B_rho_odr);
            T.B_rho_prime(ii, :) = T.Bh_rho_prime(X, T.B_rho_dof, T.B_rho_odr);
            for ii=2:T.nip
                X = T.Xs(ii);
                T.B_xi((ii-1)*6+1:ii*6, :) = T.Bh_xi(X, T.B_xi_dof, T.B_xi_odr);
                X = T.Xs(ii-1)+Z1*(T.Xs(ii)-T.Xs(ii-1));
                T.B_Z1_xi((ii-2)*6+1:(ii-1)*6, :) = T.Bh_xi(X, T.B_xi_dof, T.B_xi_odr);
                X = T.Xs(ii-1)+Z2*(T.Xs(ii)-T.Xs(ii-1));
                T.B_Z2_xi((ii-2)*6+1:(ii-1)*6, :) = T.Bh_xi(X, T.B_xi_dof, T.B_xi_odr);
                X = T.Xs(ii-1)+Z*(T.Xs(ii)-T.Xs(ii-1));
                T.B_Z_xi((ii-2)*6+1:(ii-1)*6, :) = T.Bh_xi(X, T.B_xi_dof, T.B_xi_odr);
                T.B_rho(ii, :) = T.Bh_rho(X, T.B_rho_dof, T.B_rho_odr);
                T.B_rho_prime(ii, :) = T.Bh_rho_prime(X, T.B_rho_dof, T.B_rho_odr);
            end
        end

        %% updates viscoleastic matrix
        function T = UpdateMEG(T)
            if isempty(T.dof_xi) || isempty(T.dof_rho)
                return
            end
            [T.Ms, T.Es, T.Gs] = MEG(T.Link, T.Xs);
            [T.sigma_xi, T.gamma_xi] = EG_xi_bar(T.Link, T.Xs);
            [T.sigma_rho, T.gamma_rho] = EG_rho_bar(T.Link, T.Xs);
        end

        %% updates initial position
        % TODO: update CUSTOMIZED initial position
        function T= UpdateXRStar(T)
            if isempty(T.dof_xi) || isempty(T.dof_rho)
                return
            end
            T.xi_star = zeros(6*T.nip, 4); % precomputation at all gauss and zannah guess points
            T.xi_star(4:6:end, :) = ones(T.nip, 4);
            T.rho_star = ones(T.nip, 1);
        end

        function T = Add_more_X(T)
            for k = 1:size(T.Xadd)
                X1 = T.Xadd(k);
                if ~any(T.Xs == X1) %if already present, don't add
                    iv = find(T.Xs > X1);
                    iX1 = iv(1);
                    T.Xs = [T.Xs(1:iX1-1); X1; T.Xs(iX1:end)];
                    T.Ws = [T.Ws(1:iX1-1); 0; T.Ws(iX1:end)];
                    T.nip = T.nip + 1;
                end
            end
        end

        function UpdateAll(T)
            if isempty(T.dof_xi) || isempty(T.dof_rho)
                return
            end
            T.Updatedof();
            T.UpdateIntegration();
            T.Add_more_X();
            T.UpdateBs();
            T.UpdateXRStar();
            T.UpdateMEG();
        end

        function UpdatePreCompute(T)
            T.UpdateBs();
            T.UpdateXRStar();
            T.UpdateMEG();
        end

        function T = set.B_xi_dof(T, val)
            T.B_xi_dof = val;
            T.Updatedof();
            T.UpdateBs();
        end

        function T = set.B_rho_dof(T, val)
            T.B_rho_dof = val;
            T.Updatedof();
            T.UpdateBs();
        end

        function T = set.B_xi_odr(T, val)
            T.B_xi_odr = val;
            T.UpdateBs();
        end

        function T = set.B_rho_odr(T, val)
            T.B_rho_odr = val;
            T.UpdateBs();
        end

        function T = set.Bh_xi(T, val)
            T.Bh_xi = val;
            T.UpdateAll();
        end

        function T = set.Bh_rho(T, val)
            T.Bh_rho = val;
            T.UpdateAll();
        end

        function T = set.Xadd(T, val)
            T.Xadd = val;
%             T.Add_more_X();
            T.UpdateAll();
        end

        function s = saveobj(T)
            % save properties into a struct to avoid set method when loading
            % Pass that struct onto the SAVE command.
            s.basisType = T.basisType;
            s.B_xi_dof = T.B_xi_dof;
            s.B_rho_dof = T.B_rho_dof;
            s.B_xi_odr = T.B_xi_odr;
            s.B_rho_odr = T.B_rho_odr;
            s.dof_xi = T.dof_xi;
            s.dof_rho = T.dof_rho;
            s.B_xi = T.B_xi;
            s.B_Z1_xi = T.B_Z1_xi;
            s.B_Z2_xi = T.B_Z2_xi;
            s.B_Z_xi = T.B_Z_xi;
            s.B_rho = T.B_rho;
            s.B_rho_prime = T.B_rho_prime;
            s.Bh_xi = T.Bh_xi;
            s.Bh_rho = T.Bh_rho;
            s.Bh_rho_prime = T.Bh_rho_prime;
            s.xi_star = T.xi_star;
            s.xi_starfn = T.xi_starfn;
            s.rho_star = T.rho_star;
            s.rho_starfn = T.rho_starfn;
            s.Link = T.Link;
            s.nip = T.nip;
            s.Xs = T.Xs;
            s.Ws = T.Ws;
            s.Ms = T.Ms;
            s.Es = T.Es;
            s.Gs = T.Gs;
            s.sigma_xi = T.sigma_xi;
            s.gamma_xi = T.gamma_xi;
            s.sigma_rho = T.sigma_rho;
            s.gamma_rho = T.gamma_rho;
            s.Xadd = T.Xadd;
        end
    end

    methods(Static)
        %% TODO: loadobj
        function T = loadobj(s)
            % Construct a new object T
            T = SorosimTwist;
            % Assign the properties saved in s to the new object T
            T.basisType = s.basisType;
            T.B_xi_dof = s.B_xi_dof;
            T.B_rho_dof = s.B_rho_dof;
            T.B_xi_odr = s.B_xi_odr;
            T.B_rho_odr = s.B_rho_odr;
            T.dof_xi = s.dof_xi;
            T.dof_rho = s.dof_rho;
            T.B_xi = s.B_xi;
            T.B_Z1_xi = s.B_Z1_xi;
            T.B_Z2_xi = s.B_Z2_xi;
            T.B_Z_xi = s.B_Z_xi;
            T.B_rho = s.B_rho;
            T.B_rho_prime = s.B_rho_prime;
            T.Bh_xi = s.Bh_xi;
            T.Bh_rho = s.Bh_rho;
            T.Bh_rho_prime = s.Bh_rho_prime;
            T.xi_star = s.xi_star;
            T.xi_starfn = s.xi_starfn;
            T.rho_star = s.rho_star;
            T.rho_starfn = s.rho_starfn;
            T.Link = s.Link;
            T.nip = s.nip;
            T.Xs = s.Xs;
            T.Ws = s.Ws;
            T.Ms = s.Ms;
            T.Es = s.Es;
            T.Gs = s.Gs;
            T.sigma_xi = s.sigma_xi;
            T.gamma_xi = s.gamma_xi;
            T.sigma_rho = s.sigma_rho;
            T.gamma_rho = s.gamma_rho;
            T.Xadd = s.Xadd;
        end
    end
end

