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

% isotropy
nu3 = -15/((E-(E-2*G)^2/(4*G))*A)+1;
rho = -(E-2*G)/(4*G)*(nu3-1)+1;

res = 4*G*A*(rho-1)+(E-2*G)*A*(nu3-1);

A = [0 0 -(E-2*G)/E;
     0 0 1;
     (E-2*G)*A/(G*I33) 4*A/I33 0];
detA = det(A);
