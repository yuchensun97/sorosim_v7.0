clc;
clear;
close all;
%% analytical solution starts here
% initialization
E = 1e6; %young
G = E/(2*(1+0.5)); %shear, incompressible
r = 0.02; % radius
L = 0.5;  % length
A = pi*r^2;
I33 = 2*(pi/4)*r^4;
Ftip = -500;
bc = 'neumann';

% hand-on solution
nu3 = Ftip/((E-(E-2*G)^2/(4*G))*A)+1;
rho = -(E-2*G)/(4*G)*(nu3-1)+1;

% numerical solution
% ds=0.01;
% x0 = [1;0]; % intial guess of rho and rho_prime
% options = optimoptions('fsolve','Algorithm','trust-region-dogleg','Display','iter','MaxFunctionEvaluations',2e10);
% x0_true = fsolve(@(x0)shooting(x0, E, G, A, I33, L, Ftip, bc), x0, options);
% [s, x] = ode45(@(s,x)compress(s, x, E, G, A, I33, Ftip), 0:ds:L, x0_true);

% function dxds = compress(s, x, E, G, A, I33, Ftip)
%     rho = x(1);
%     rhod = x(2);
%     nu_hat = Ftip/(E*A)-(E-2*G)/E*(rho-1); %nu3-1
%     rhodd = (4*G*A*(rho-1)+(E-2*G)*A*nu_hat)/(G*I33);
%     dxds = [rhod;rhodd];
% end

% function err = shooting(x0, E, G, A, I33, L, Ftip, bc)
%     options = odeset('RelTol',1e-3);
%     ds = 0.01;
%     [s, x] = ode45(@(s,x)compress(s, x, E, G, A, I33, Ftip), 0:ds:L, x0, options);
%     if bc == 'neumann' %#ok<BDSCA>
%         rhod0 = x(1, 2);
%         rhodl = x(end,2);
%         err = [rhod0; rhodl];
%     elseif bc == 'dirichlet' %#ok<BDSCA>
%         rho0 = x(1, 1);
%         rhol = x(1, end);
%         err = [rho0; rhol];
%     else
%         error('incorrect bounddary condition');
%     end
% end
