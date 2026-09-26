%%%%% Accelerated Gradient method for EIT
%%%%% Complete Electrode Model (CEM) Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% May 2025, Sanwar Uddin Ahmad, sahmad@cvsu.edu, VSU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Creating mesh and true sigma model
% Initialize
clear, close all
format long;
%%
body.rc = 0.15; %0.15; % radius of the circle (m)
body.NumSrc = 32; % number of the electrodes
body.zl = 1.2; %*ones(1,body.NumSrc); % contact impedance (Ohm)
body.el = 0.022; % length of the electrode (m)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  %
% Note: el needs to be determined manually after generating the mesh   %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Specifying data
% Create geometry of boundary
cr = 1; %rc;
gd=[4; 0.0; 0.0; cr; cr; 0.0];
dl=decsg(gd);

%%%%%%%%%% MAKE SURE (number of bdy nodes) >= 2*number of electrodes %%%%%%%%%%
hval = 0.15;
% Create mesh
[p,e,t]=initmesh(dl,'hmax',hval);
[p,e,t]=refinemesh(dl,p,e,t);
nodes1=length(p);
[p1,e1,t1]=refinemesh(dl,p,e,t);
p1 = p1*body.rc;
p = p*body.rc;
nodes=length(p1);
tris=length(t1);
sigTrue=ones(nodes,1); % conductivity sigma at refined mesh 
sigTrue1=ones(nodes1,1); % sigma at coarser mesh

inArea= 0.3; 
outArea= 0.07; 
%% Finding electrode positions
electrode = elec_ind(body,p,e);
body.elec = electrode;
nl = body.NumSrc;

%% True conductivity
% [sigTrue,sigTrue1] = ex1(p,p1,inArea,outArea,body); % single
[sigTrue,sigTrue1] = ex2(p,p1,inArea,outArea,body); % double
% [sigTrue,sigTrue1] = ex33(p,p1,inArea,outArea,body); % one circle one rectangle one tringle
figure;pdeplot(p,e,t,'xydata',sigTrue1,'mesh','off');colormap(jet)

%% Current type
amp = 0.2; %5 mA
Ic = trig_current(body,amp);
npat = size(Ic,2);
body.npat = npat;
body.current = Ic;   

%% MIRGN Method for coarser mesh

sig_guess = outArea*ones(size(p,2),1);
sig_g = pdeintrp(p,t,sig_guess);
sig_T = pdeintrp(p,t,sigTrue1');
load('ex1_firgn_0_1n_new.mat');

%% Computing the jacobian

[~,u,Kn] = fwd_solver_eit2D(p,e,t,sig_g',body);
J = cem_jac_vauk(p,t,u,inv(Kn),body);
J1 = J;

%% Reconstruction
Vmeas = ex1_firgn_0_1n_new.data;
lambda = ex1_firgn_0_1n_new.lambda;
numbit = 100;
Tol = ex1_firgn_0_1n_new.Tol;
l2err_sig = zeros(numbit,1); % Reconstruction error

% [sig_all,l2err_sig,st,obj,iter] = IRGN(Vmeas,sig_g,sig_T,J1,lambda,numbit,p,e,t,Tol,body); % running IRGN method
% [sig_all,l2err_sig,st,obj,iter] = block_kaczmarz(Vmeas,sig_g,sig_T,J1,lambda,numbit,p,e,t,Tol,body); % running bART method

%% Displaying last sigma model and predicted voltage in comparison to experimental data
k = iter-1;
figure; pdeplot(p,e,t,'xydata',sig_all(:,k),'mesh','off'); colormap(jet); 
