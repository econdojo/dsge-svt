% Function USER_MOD
%
% Purpose:    User input - model & measurement equations
%
% Format:     user_mod
%
% Input:      P         parameter structure
%             V         variable structure
%             para      model parameters
%             ssp       steady state/implied parameters
%
% Output:     G0,G1,Psi,Pi,CC - state equation (gensys.m input)
%
%                       G0*x(t) = G1*x(t-1) + Psi*epsilon(t) + Pi*eta(t)
%
%             Z,D,Sigma_u     - measurement equation
%
%                       y(t) = D + Z*x(t) + u(t), E[uu'] = Sigma_u (diagonal)
%
% Written by Fei Tan, Saint Louis University
% Updated: December 15, 2016

%% User Input Start Here

%---------------------------------------------
%              GENSYS Equations
%---------------------------------------------

%---------------------------------------------
j = j+1;   % (1) DIS
%---------------------------------------------
G0(j,V.y) = 1;
G0(j,V.E_y) = -1;
G0(j,V.g) = -(1-para(P.rho_G));
G0(j,V.R) = 1/para(P.tau);
G0(j,V.E_pi) = -1/para(P.tau);
G0(j,V.z) = -para(P.rho_Z)/para(P.tau);
%---------------------------------------------
j = j+1;   % (2) NKPC
%---------------------------------------------
G0(j,V.pi) = 1;
G0(j,V.E_pi) = -ssp(P.beta);
G0(j,V.y) = -para(P.kappa);
G0(j,V.g) = para(P.kappa);
%---------------------------------------------
j = j+1;   % (3) MP
%---------------------------------------------
% Output gap rule; unimodal
G0(j,V.R) = 1;
G0(j,V.pi) = -(1-para(P.rho_R))*para(P.psi_1);
G0(j,V.y) = -(1-para(P.rho_R))*para(P.psi_2);
G0(j,V.g) = (1-para(P.rho_R))*para(P.psi_2);
G1(j,V.R) = para(P.rho_R);
Psi(j,V.eps_R) = 1;
%---------------------------------------------
j = j+1;   % (4) G shock
%---------------------------------------------
G0(j,V.g) = 1;
G1(j,V.g) = para(P.rho_G);
Psi(j,V.eps_G) = 1;
%---------------------------------------------
j = j+1;   % (5) Z shock
%---------------------------------------------
G0(j,V.z) = 1;
G1(j,V.z) = para(P.rho_Z);
Psi(j,V.eps_Z) = 1;
%---------------------------------------------
j = j+1;   % (6) E_y error
%---------------------------------------------
G0(j,V.y) = 1;
G1(j,V.E_y) = 1;
Pi(j,V.fe_y) = 1;
%---------------------------------------------
j = j+1;   % (7) E_pi error
%---------------------------------------------
G0(j,V.pi) = 1;
G1(j,V.E_pi) = 1;
Pi(j,V.fe_pi) = 1;
%---------------------------------------------
j = j+1;   % (8) Lagged y
%---------------------------------------------
G0(j,V.d_y) = 1;
G1(j,V.y) = 1;

%---------------------------------------------
%           Measurement Equations
%---------------------------------------------
Z(V.YGR,V.y) = 1;
Z(V.YGR,V.z) = 1;
Z(V.YGR,V.d_y) = -1;
Z(V.INF,V.pi) = 4;
Z(V.INT,V.R) = 4;
D(V.YGR) = para(P.gamma_Q);
D(V.INF) = para(P.pi_A);
D(V.INT) = para(P.pi_A)+para(P.r_A)+4*para(P.gamma_Q);
% Measurement covariance (comment out if no measurement error)
Sigma_u(V.YGR) = 0.1207^2;
Sigma_u(V.INF) = 0.2919^2;
Sigma_u(V.INT) = 0.4476^2;
%Sigma_u = (0.2*std(Y)').^2;

%% User Input End Here