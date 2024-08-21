%This toolbox allows the modelling and simulation of static and dynamic
%soft robot open chains
%Run this file to initialize the Rigid-Soft Robotics Toolbox
%Last modified by Anup Teejo Mathew - 18/01/2022
% Modified by Yuchen SUN <sun.yuchen22@u.nus.edu>

clc
clear variables

addpath(genpath('Basicfunctions'))
addpath(genpath('Custom'))
addpath('SorosimLinkFiles')
addpath('SorosimTwistFiles')
addpath(genpath('SorosimLinkageFiles')) %include subfolders
addpath(genpath('tests'))


if exist('.\Basis_properties.mat','file')
    delete('Basis_properties.mat')
end
if exist('.\cableactuation.mat','file')
    delete('cableactuation.mat')
end
if exist('.\CablePoints.mat','file')
    delete('CablePoints.mat')
end

disp('Welcome to SoRoSim Toolbox')
disp('Provide the soft link properties file. See ./Custom/Octopus.json for an example')
disp('Type LinkName=SorosimLink(filename) to create the link (joint and body)')
disp('Type LinkageName=SorosimLinkage(LinkName) to create linkages')
disp('For static equilibrium problem type [q,u]=LinkageName.statics')
disp('For dynamics problem type [t,qqd] = LinkageName.dynamics')
