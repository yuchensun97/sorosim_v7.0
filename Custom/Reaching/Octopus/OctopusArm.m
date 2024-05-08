classdef OctopusArm < SorosimLinkage
    % OctopusArm inherited from SorosimLinkage
    % Deal with the case when the actuation loads are functions of both time and space
    methods
        function Tr = OctopusArm(Link, varargin)
            Tr@SorosimLinkage(Link, varargin);
        end
    end
    methods
        % @Override
        ydot = derivatives(Tr, t, qqd, uqt_xi, uqt_rho);
        err = equilibrium(Tr, qu, u_xi, u_rho);
        Bq_xi = ComputeCableActuation(Tr, dc, dcp, q_xi, q_rho, uqt_xi);
        [t, qqd] = dynamics(Tr, qqd0, uqt_xi, uqt_rho, odetype, dt, tmax);
        q = statics(Tr, qu0, u_xi, u_rho)
    end
end
