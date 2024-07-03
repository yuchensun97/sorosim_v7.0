classdef SorosimLinkage
    %SOROSIMLINKAGE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %%
        %General Properties
        ndof_xi     %number of DOF for twist
        ndof_rho    %number of DOF for inflation ratio
        nsig        %number of configuration computed
        Link        % SorosimLink Object
        Twists      % SorosimTwist vector

        g_base      %(4x4)initial pose of the link wrt to the inertial frame
        rho_base    %initial inflation ratio
        Z_order     %order of Zannah collocation (2, 4 or 6) defualt value is 4

        %External Force Properties
        Gravity     %logical 1 if gravity is present, 0 if not. (default value: 1)
        G           %gravity vector (default value: [0 0 0 0 0 -9.81])

        %Pre-computed elastic Properties
        Damped      %1 if the soft links are elastically damped 0 if not
        K_xi        %Generalized stiffness matrix for strains
        D_xi        %Generalized damping matrix for strains
        K_xi_bar    %Generalized stiffness matrix for the \rho term in strain equation
        D_xi_bar    %Generalized damping matrix for the \rho term in strain equation
        K_rho_part  %Generalized stiffness matrix for lateral equation
        D_rho       %Generalized damping matrix for lateral equation
        K_rho_bar   %Generalized stiffness matrix for the strain term in the lateral equation
        D_rho_bar   %Generalized damping matrix for the strain term in the lateral equation
        M_rho       %Generalized mass matrix for the lateral equation.

        % actuation
        ActuatedL    %1 if the soft links are actuated in longitudinal direction, 0 if not
        ActuatedR    %1 if the soft links are actuated in radial direction, 0 if not

        % environment
        Water       %1 if the soft link is immersed in water
        rho_w       %density of water
        DL          %drag/lifting matrix in referenced configuration
        M_added     %added mass matrix in referenced configuration

        % cable actuator for soft links
        n_sact       %number of soft link actuators
        dc           %(n_sactx1) cells of local cable position (0, yp, zp) at Gauss quadrature points of all active soft divisions
        dcp          %(n_sactx1) cells of space derivative of the local cable position (0, yp',zp')
        CableActuator     %CableActuation class of parameterized functions corresponding to the y and z coodinates of the cable

        % radial actuator for soft links
        n_ract      %number of radial actuators
        rc          %(n_ractx1) array of radial actuator position
        RadialActuator %RadialActuation classs 

        % custom actuation
        CAP          %true if custom actuation is preseent, false is not
        CAS          %true to apply a custom actuator strength

        % point force
        PointForce  %1 if the soft links are actuated, 0 if not
        LocalForce  %1 if point force/moment is a local force, and 0 if it is wrt global frame
        np          %number of point forces
        Fp_loc      %cell, location of the point force/moment
        Fp_vec      %cell, value of the point force/moment

        %% Plotting Properties
        PlotParameters    %struct containing plotting parameters
        %%%PlotParameters is a struct with following elements%%%
        %Lscale
        %CameraPosition     CameraPosition with respect to origin
        %CameraUpVector     Orientation of normal
        %CameraTarget       Target location
        %Light              logical 1 if the light is on, 0 if not. (default value: 1)
        %Az_light           Light Azimuth wrt camera
        %El_light           Light Elevation wrt camera
        %X_lim              X limit [X_lower_limt X_upper_limit]
        %Y_lim              Y limit [Y_lower_limt Y_upper_limit]
        %Z_lim              Z limit [Z_lower_limt Z_upper_limit]
        %FrameRateValue     FrameRate for dyanmic plot
        %ClosePrevious      logical 0 to not close previous image, 1 to close. (default value: 1)
    end
    
    methods
        function Tr = SorosimLinkage(Link, varargin)
            %SorosimLinkage Constructor

            %% initialization starts here
            Tr.Link = Link;
            Tr.g_base = Link.gi;
            Tr.rho_base = Link.rhoi;
            Lscale = Link.L;
            Tr.Z_order = 4;

            % Twists for base and the body
            % VTwists(1): fix base Twist
            % VTwists(2): Soft body Twist
            VTwists = SorosimTwist.empty(Tr.Link.npie, 0);
            VTwists(1) = SorosimTwist([], []);
            VTwists(1).dof_xi = 0;
            VTwists(1).dof_rho = 0;

            VTwists(2) = SorosimTwist(Link, Link.B_xi, Link.B_rho, Link.basisType);
            Tr.Twists = VTwists;
            Tr.ndof_xi = VTwists(1).dof_xi + VTwists(2).dof_xi;
            Tr.ndof_rho =VTwists(1).dof_rho + VTwists(2).dof_rho;
            Tr.nsig = VTwists(2).nip;

            %% input parser
            p = inputParser;
            checkLink = @(x)isa(x, 'SorosimLink');
            defaultGravity = false;
            defaultDamping = false;
            defaultWater = false;
            defaultActuationL = false;
            defaultActuationR = false;
            defaultPointForce = false;
            defaultLocalForce = false(0, 0);
            defaultFp_loc = zeros(0, 0); % integration point location
            defaultFp_vec = cell(0, 0); % force should be function handler
            defaultCableActuator = CableActuation();
            checkCableActuation = @(x)isa(x, 'CableActuation');
            defaultRadialActuator = RadialActuation();
            checkRadialActuator = @(x)isa(x, 'RadialActuation');

            addRequired(p, 'Link', checkLink);
            addOptional(p, 'Damped', defaultDamping, @islogical);
            addOptional(p, 'Gravity', defaultGravity, @islogical);
            addOptional(p, 'PointForce', defaultPointForce, @islogical);
            addOptional(p, 'Water', defaultWater, @islogical);
            addOptional(p, 'ActuationL', defaultActuationL, @islogical);
            addOptional(p, 'ActuationR', defaultActuationR, @islogical);
            addParameter(p, 'Fp_loc', defaultFp_loc, @isvector);
            addParameter(p, 'LocalForce', defaultLocalForce, @islogical);
            addParameter(p, 'Fp_vec', defaultFp_vec, @iscell);
            addParameter(p, 'CableActuator', defaultCableActuator, checkCableActuation);
            addParameter(p, 'RadialActuator', defaultRadialActuator, checkRadialActuator);

            parse(p, Link, varargin{:});

            %% External Force Properties
            Tr.Gravity = p.Results.Gravity;
            if Tr.Gravity
                Tr.G = [0 0 0 0 0 -9.81]';
            end

            Tr.Water = p.Results.Water;
            if Tr.Water
                Tr.rho_w = 1000;
                Tr.DL = dragging(Tr);
                Tr.M_added = addedMass(Tr);
            end

            %% Constant coefficients
            % stiffness
            K_xi = findK_xi(Tr);
            K_xi_bar = findK_xi_bar(Tr);
            K_rho_part = findK_rho_part(Tr);
            K_rho_bar = findK_rho_bar(Tr);
            Tr.K_xi = K_xi;
            Tr.K_xi_bar = K_xi_bar;
            Tr.K_rho_part = K_rho_part;
            Tr.K_rho_bar = K_rho_bar;

            %damping
            D_xi = findD_xi(Tr);
            D_xi_bar = findD_xi_bar(Tr);
            D_rho = findD_rho(Tr);
            D_rho_bar = findD_rho_bar(Tr);
            Tr.D_xi = D_xi;
            Tr.D_xi_bar = D_xi_bar;
            Tr.D_rho_bar = D_rho_bar;
            Tr.D_rho = D_rho;
            Tr.Damped = p.Results.Damped;

            % mass
            M_rho = findM_rho(Tr);
            Tr.M_rho = M_rho;

            %% point force
            Tr.PointForce = p.Results.PointForce;
            if Tr.PointForce
                Fp_loc = p.Results.Fp_loc;
                LocalForce = p.Results.LocalForce;
                Fp_vec = p.Results.Fp_vec;
                nFp = length(Fp_loc);
                Tr.np = nFp;

                sz_local = size(LocalForce);
                if sz_local(1)~=nFp
                    error('Mismatch size of LocalForce and Fp_loc');
                end
                if sz_local(2)~=1
                    error('dimension exceed');
                end

                sz_vec = size(Fp_vec);
                if sz_vec(1)~=nFp
                    error('Mismatch size of Fp_vec and Fp_loc');
                end
                if sz_vec(2)~=1
                    error('dimension exceed');
                end
                for ip = 1:nFp
                    Fp_vec_h = Fp_vec{ip};
                    if ~isa(Fp_vec_h, 'function_handle')
                        error('Fp_vec should be a function handle');
                    end
                    Fp_vec_s = func2str(Fp_vec_h);
                    if ~startsWith(Fp_vec_s, '@(t)')
                        error('The function handler should be the form of @(t)...');
                    end
                end

                Tr.Fp_loc = Fp_loc;
                Tr.Twists(2).Xadd = Fp_loc;
                Tr.nsig = Tr.Twists(2).nip;
                Tr.Fp_vec = Fp_vec;
                Tr.LocalForce = LocalForce;
            end

            %% Actuation
            Tr.ActuatedL = p.Results.ActuationL;
            Tr.ActuatedR = p.Results.ActuationR;
            Tr.n_sact = 0;
            Tr.n_ract = 0;
            if Tr.ActuatedR
                RadialActuator = p.Results.RadialActuator;
                n_ract = RadialActuator.get_n_ract();
                if n_ract ==0
                    error('You must have at least 1 radial actuator');
                end
                rc = RadialActuator.get_rc();
                Tr.rc = rc;
                Tr.n_ract = n_ract;
            end

            if Tr.ActuatedL
                CableActuator = p.Results.CableActuator;
                n_sact = CableActuator.get_n_sact();
                if n_sact == 0
                    error('You must have at least 1 cable actuator');
                end
                dc_fn = CableActuator.get_dc_fn();
                dcp_fn = CableActuator.get_dcp_fn();
                dc = cell(n_sact, 1);
                dcp = cell(n_sact, 1);
                Xs = Tr.Twists(2).Xs;
                for i=1:n_sact
                    dc_fn_here = dc_fn{i};
                    dcp_fn_here = dcp_fn{i};
                    dc_sact = arrayfun(dc_fn_here, Xs, 'UniformOutput', false);
                    dcp_sact = arrayfun(dcp_fn_here, Xs, 'UniformOutput', false);
                    dc{i} = cell2mat(dc_sact');
                    dcp{i} = cell2mat(dcp_sact');
                end
                Tr.dc = dc;
                Tr.dcp = dcp;
                Tr.n_sact = n_sact;
            end

            %% Plot parameters
            PlotParameters.Lscale         = Lscale;
            PlotParameters.CameraPosition = [Lscale*2 -Lscale/2 Lscale/2];
            PlotParameters.CameraTarget   = [0 0 0];
            PlotParameters.CameraUpVector = [0 0 1];
            PlotParameters.Light          = true;
            PlotParameters.Az_light       = 0;
            PlotParameters.El_light       = 0;
            PlotParameters.X_lim          = [-0.8*Lscale 1.2*Lscale];
            PlotParameters.Y_lim          = [-0.5*Lscale 1.2*Lscale];
            PlotParameters.Z_lim          = [-0.5*Lscale 0.5*Lscale];
            PlotParameters.FrameRateValue = 50;
            PlotParameters.ClosePrevious  = false;
            PlotParameters.Position       = [0.1300 0.1100 0.7750 0.8150]; %default value (normalized)

            Tr.PlotParameters = PlotParameters;

        end %end of class constructor
    end

    methods
        %% Methods
        [g, rho] = FwdKinematics(Tr, q_xi, q_rho);
        fh = plotq(Tr, q_xi, q_rho);
        [J_xi, J_rho] = Jacobian(Tr, q_xi);
        Jd_xi = Jacobiandot(Tr, q_xi, qd_xi);
        Jp_rho = Jacobianprime(Tr);
        K_xi = findK_xi(Tr);
        K_xi_bar = findK_xi_bar(Tr);
        K_rho_part = findK_rho_part(Tr);
        K_rho_bar = findK_rho_bar(Tr);
        D_xi = findD_xi(Tr);
        D_xi_bar = findD_xi_bar(Tr);
        D_rho = findD_rho(Tr);
        D_rho_bar = findD_rho_bar(Tr);
        M_rho = findM_rho(Tr);
        DL = dragging(Tr);
        M_added = addedMass(Tr);
    end

    methods
        %% setter
        function Tr = Update(Tr)
            if isempty(Tr.Z_order)
                return
            end
            Tr.ndof_xi = Tr.Twists(1).dof_xi + ...
                        Tr.Twists(2).dof_xi;
            Tr.ndof_rho = Tr.Twists(1).dof_rho + ...
                        Tr.Twists(2).dof_rho;
            Tr.nsig = Tr.Twists(2).nip;
        end

        function Tr = UpdateTwist(Tr)
            Tr.Twists(2) = SorosimTwist(Tr.Link, Tr.Link.B_xi, Tr.Link.B_rho);
            Tr.Update();
        end

        function Tr = set.Twists(Tr, val)
            if ~isa(val,'SorosimTwist')
                error('Input must be a SorosimTwist object')
            end
            Tr.Twists = val;
            Tr.Update();
        end

    end
end
