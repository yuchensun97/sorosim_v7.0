classdef SorosimLink
    %SorosimLink rigid joint and a body
    %   author: Yuchen Sun 
    %   modified from Anup Mathew 02.03.2022
    
    properties
        %General Properties
        npie=2;    %number of pieces. 2 for soft link

        %Geometric Properties
        L          %length of the link
        r_fn          %radius as a function of X1 (X1=X/L, X1 varies from 0 to 1)[m]
        r_base     %radius at the base [m]
        r_tip      %radius at the tip [m]
        gi=eye(4); %Transformation from joint to center of area
        gf=eye(4); %Transformation to joint from center of area
        rhoi=1;    %inflation ratio at the joint, scalar

        %Material Properties
        E          %Young's modulus [Pa]
        Poi        %Poisson's ratio, scalar
        G          %Shear modulus [Pa]
        Lamda      %first Lame constant [Pa]
        Rho0       %density of the material [kg/m^3]
        Eta        %viscosity [Pa.s]

        %Motion Properties
        B_xi       %motions allowed and their order, (6x2)
        B_rho      %allowable inflation and its order, (1x2)

        %basis Properties
        basisType  %basis type for rho

        %Plot Properties
        color      %color of the link (random by default)
        n_l        %number of cross sections per link
        n_r        %number of points per cross section

        Lscale     %scaling factor for plotting symbols or axes
    end
    
    methods
        function Li = SorosimLink(varargin)
            %Sorosim Link constructor
            if nargin == 1
                filename = varargin{1};
                %read from json file
                fid = fopen(filename);
                raw = fread(fid,inf);
                str = char(raw');
                fclose(fid);
                data = jsondecode(str);

                % check if input values are valid
                if data.length <= 0
                    error('length must be positive')
                end
                if data.base_radius <= 0
                    error('base_radius must be positive')
                end
                if data.tip_radius < 0
                    error('tip_radius must be positive')
                end
                if data.young <= 0
                    error('young must be positive')
                end
                if data.poisson < 0 || data.poisson > 0.5
                    error('poisson must be between 0 and 0.5')
                end
                if data.density <= 0
                    error('density must be positive')
                end
                if data.viscousity < 0
                    error('viscousity must be positive')
                end
                if data.cs <= 0
                    error('cross sections per link must be positive')
                end
                if data.points <= 0
                    error('points per cross section must be positive')
                end
                if size(data.motion, 1)~=6 || size(data.motion, 2)~=1
                    error('motion must be a (6x1) vector')
                end
                if ~all(data.motion == 0 | data.motion == 1)
                    error('motion must be a binary vector')
                end
                if size(data.motion_order, 1)~=6 || size(data.motion_order, 2)~=1
                    error('motion_order must be a (6x1) vector')
                end
                if ~all(data.motion_order >= 0)
                    error('motion_order must be a positive vector')
                end
                if ~(data.inflation==0 || data.inflation==1)
                    error('inflation must be 0 or 1')
                end
                if data.inflation_order <= 0
                    error('inflation_order must be positive')
                end

                %assign values
                Li.basisType = data.basisType;
                Li.L = data.length;
                r_base = data.base_radius;
                r_tip = data.tip_radius;
                Li.r_base = r_base;
                Li.r_tip = r_tip;
                Li.r_fn = @(X1) X1.*(r_tip-r_base) + r_base;
                Li.B_xi = [data.motion data.motion_order];
                Li.B_rho = [data.inflation data.inflation_order];
                Li.E = data.young;
                Li.Poi = data.poisson;
                Li.G = Li.E/(2*(1+Li.Poi));
                Li.Lamda = Li.E*Li.Poi/((1+Li.Poi)*(1-2*Li.Poi));
                Li.Rho0 = data.density;
                Li.Eta = data.viscousity;
                Li.color = rand(1,3);
                Li.n_l = data.cs;
                Li.n_r = data.points;
                A0 = pi*r_base^2;

            elseif nargin == 0
                Li.L = 0.5;
                Li.r_base = 0.02;
                Li.r_tip = 0.02;
                Li.r_fn = @(X1) X1.*(Li.r_tip - Li.r_base) + Li.r_base;
                Li.E = 1e6;
                Li.Poi = 0.49;
                Li.G = Li.E/(2*(1+Li.Poi));
                Li.Rho0 = 1000;
                Li.Eta = 1e4;
                Li.Lamda = Li.E*Li.Poi/((1+Li.Poi)*(1-2*Li.Poi));
                Li.color = rand(1,3);
                Li.n_l = 25;
                Li.n_r = 18;
                r_base = Li.r_fn(0);
                A0 = pi*r_base^2;
                Li.B_xi = [1 1 1 1 1 1;
                           0 0 0 0 0 0]';
                Li.B_rho = [0, 0];
                Li.basisType = 'mixed';
            else
                error('Wrong number of input arguments')
            end

            Lscale = (A0 * Li.L)^(1/3);
            Li.Lscale = Lscale;     
        end
    end

    methods
        function Li = Update(Li)
            r_fn_new = @(X1) X1.*(Li.r_tip-Li.r_base) + Li.r_base;
            Li.r_fn = r_fn_new;

            A0 = pi*Li.r_base^2;
            Lscale_now = (A0 * Li.L)^(1/3);
            Li.Lscale = Lscale_now;
        end

        function Li = UpdateG(Li)
            if isempty(Li.color)
                return
            end
            Li.G = Li.E/(2*(1+Li.Poi));
        end

        %% set properties
        function Li = set.r_base(Li, val)
            Li.r_base = val;
            Li.Update();
        end

        function Li = set.r_tip(Li, val)
            Li.r_tip = val;
            Li.Update();
        end

        function Li = set.L(Li, val)
            Li.L = val;
            Li.Update();
        end

        function Li = set.Rho0(Li, val)
            Li.Rho0 = val;
            Li.Update();
        end

        function Li = set.E(Li, val)
            Li.E = val;
            Li.UpdateG();
        end

        function Li = set.Poi(Li, val)
            Li.Poi = val;
            Li.UpdateG();
        end

        function Li = set.B_xi(Li, val)
            Li.B_xi = val;
        end

        function Li = set.B_rho(Li, val)
            Li.B_rho = val;
        end

    end
end
