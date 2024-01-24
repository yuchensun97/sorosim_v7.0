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
        function Tr = SorosimLinkage(Link)
            %SorosimLinkage Constructor
            
            % check input
            if nargin < 1
                error('Not enough input arguments')
            end
            if ~isa(Link,'SorosimLink')
                error('Input must be a SorosimLink object')
            end

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

            VTwists(2) = SorosimTwist(Link, Link.B_xi, Link.B_rho);
            Tr.Twists = VTwists;
            Tr.ndof_xi = VTwists(1).dof_xi + VTwists(2).dof_xi;
            Tr.ndof_rho =VTwists(1).dof_rho + VTwists(2).dof_rho;
            Tr.nsig = VTwists(2).nip;

            %% Plot parameters
            PlotParameters.Lscale         = Lscale;
            PlotParameters.CameraPosition = [Lscale*2 -Lscale/2 Lscale/2];
            PlotParameters.CameraTarget   = [0 0 0];
            PlotParameters.CameraUpVector = [0 0 1];
            PlotParameters.Light          = true;
            PlotParameters.Az_light       = 0;
            PlotParameters.El_light       = 0;
            PlotParameters.X_lim          = [-2*Lscale 2*Lscale];
            PlotParameters.Y_lim          = [-2*Lscale 2*Lscale];
            PlotParameters.Z_lim          = [-2*Lscale 2*Lscale];
            PlotParameters.FrameRateValue = 50;
            PlotParameters.ClosePrevious  = false;
            PlotParameters.Position       = [0.1300 0.1100 0.7750 0.8150]; %default value (normalized)

            Tr.PlotParameters = PlotParameters;

            % ignore other properties for Kinematics

        end %end of class constructor
    end

    methods
        %% Methods
        [g, rho] = FwdKinematics(Tr, q_xi, q_rho);
        fh = plotq(Tr, q_xi, q_rho);
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
