classdef SorosimLink
    %SorosimLink rigid joint and a body
    %   modified from Anup Mathew 02.03.2022
    
    properties
        %General Properties
        npie=2;    %number of pieces. 2 for soft link

        %Geometric Properties
        L          %length of the link
        r          %radius as a function of X1 (X1=X/L, X1 varies from 0 to 1)[m]
        gi=eye(4); %Transformation from joint to center of area
        gf=eye(4); %Transformation to joint from center of area
        rho=1;     %inflation ratio at the joint, scalar

        %Material Properties
        E          %Young's modulus [Pa]
        Poi        %Poisson's ratio, scalar
        G          %Shear modulus [Pa]
        Rho0       %density of the material [kg/m^3]
        Eta        %viscosity [Pa.s]

        %Plot Properties
        color      %color of the link (random by default)
        n_l        %number of cross sections per link
        n_p        %number of points per cross section

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
                %assign values
                Li.L = data.length;
                Li.r = str2func(['@(X1)' data.radius]);
                Li.E = data.young;
                Li.Poi = data.poisson;
                Li.G = Li.E/(2*(1+Li.Poi));
                Li.Rho0 = data.density;
                Li.Eta = data.viscousity;
                Li.color = data.color;
                Li.n_l = data.cs;
                Li.n_p = data.points;
                Li.Lscale = data.scale;

            elseif nargin == 0
                %default values, same as SimpleLinkage.L1
                Li.L = 0.5;
                Li.r = @(X1) X1.*(-1.0/1.0e+2) + 3.0/1.0e+2;
                Li.E = 1e6;
                Li.Poi = 0.5;
                Li.G = Li.E/(2*(1+Li.Poi));
                Li.Rho0 = 1000;
                Li.Eta = 1e-3;
                Li.color = rand(1,3);
                Li.n_l = 25;
                Li.n_p = 10;
                Li.Lscale = 0.1122;
            else
                error('Wrong number of input arguments')
            end
                
        end
    end
end

