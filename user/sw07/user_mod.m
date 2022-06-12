% Function USER_MOD
%
% Purpose:    User input - model & measurement equations; see
%             Smets & Wouters (2007) for model details
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
j = j+1;   % (1)
%---------------------------------------------
G0(j,V.y) = 1;
G0(j,V.c) = -ssp(P.c_y);
G0(j,V.i) = -ssp(P.i_y);
G0(j,V.z) = -ssp(P.z_y);
G0(j,V.eps_g) = -1;
%---------------------------------------------
j = j+1;   % (2)
%---------------------------------------------
G0(j,V.y_star) = 1;
G0(j,V.c_star) = -ssp(P.c_y);
G0(j,V.i_star) = -ssp(P.i_y);
G0(j,V.z_star) = -ssp(P.z_y);
G0(j,V.eps_g) = -1;
%---------------------------------------------
j = j+1;   % (3)
%---------------------------------------------
G0(j,V.c) = 1;
G0(j,V.E_c) = -1/(1+para(P.h)/ssp(P.gamma));
G0(j,V.l) = -ssp(P.wl_c)*(para(P.sigma_c)-1)/para(P.sigma_c)/(1+para(P.h)/ssp(P.gamma));
G0(j,V.E_l) = -G0(j,V.l);
G0(j,V.r) = (1-para(P.h)/ssp(P.gamma))/para(P.sigma_c)/(1+para(P.h)/ssp(P.gamma));
G0(j,V.E_pi) = -G0(j,V.r);
G0(j,V.eps_b) = G0(j,V.r);
G1(j,V.c) = 1+G0(j,V.E_c);
%---------------------------------------------
j = j+1;   % (4)
%---------------------------------------------
G0(j,V.c_star) = 1;
G0(j,V.E_c_star) = -1/(1+para(P.h)/ssp(P.gamma));
G0(j,V.l_star) = -ssp(P.wl_c)*(para(P.sigma_c)-1)/para(P.sigma_c)/(1+para(P.h)/ssp(P.gamma));
G0(j,V.E_l_star) = -G0(j,V.l_star);
G0(j,V.r_star) = (1-para(P.h)/ssp(P.gamma))/para(P.sigma_c)/(1+para(P.h)/ssp(P.gamma));
G0(j,V.eps_b) = G0(j,V.r_star);
G1(j,V.c_star) = 1+G0(j,V.E_c_star);
%---------------------------------------------
j = j+1;   % (5)
%---------------------------------------------
G0(j,V.i) = 1;
G0(j,V.E_i) = -ssp(P.beta)*ssp(P.gamma)^(1-para(P.sigma_c))/(1+ssp(P.beta)*ssp(P.gamma)^(1-para(P.sigma_c)));
G0(j,V.q) = -1/para(P.varphi)/ssp(P.gamma)^2/(1+ssp(P.beta)*ssp(P.gamma)^(1-para(P.sigma_c)));
G0(j,V.eps_i) = -1;
G1(j,V.i) = 1+G0(j,V.E_i);
%---------------------------------------------
j = j+1;   % (6)
%---------------------------------------------
G0(j,V.i_star) = 1;
G0(j,V.E_i_star) = -ssp(P.beta)*ssp(P.gamma)^(1-para(P.sigma_c))/(1+ssp(P.beta)*ssp(P.gamma)^(1-para(P.sigma_c)));
G0(j,V.q_star) = -1/para(P.varphi)/ssp(P.gamma)^2/(1+ssp(P.beta)*ssp(P.gamma)^(1-para(P.sigma_c)));
G0(j,V.eps_i) = -1;
G1(j,V.i_star) = 1+G0(j,V.E_i_star);
%---------------------------------------------
j = j+1;   % (7)
%---------------------------------------------
G0(j,V.q) = 1;
G0(j,V.E_q) = -ssp(P.beta)*(1-ssp(P.delta))/ssp(P.gamma)^para(P.sigma_c);
G0(j,V.E_r_k) = -1-G0(j,V.E_q);
G0(j,V.r) = 1;
G0(j,V.E_pi) = -1;
G0(j,V.eps_b) = 1;
%---------------------------------------------
j = j+1;   % (8)
%---------------------------------------------
G0(j,V.q_star) = 1;
G0(j,V.E_q_star) = -ssp(P.beta)*(1-ssp(P.delta))/ssp(P.gamma)^para(P.sigma_c);
G0(j,V.E_r_k_star) = -1-G0(j,V.E_q_star);
G0(j,V.r_star) = 1;
G0(j,V.eps_b) = 1;
%---------------------------------------------
j = j+1;   % (9)
%---------------------------------------------
G0(j,V.y) = 1;
G0(j,V.k_s) = -para(P.Phi)*para(P.alpha);
G0(j,V.l) = -para(P.Phi)*(1-para(P.alpha));
G0(j,V.eps_a) = -para(P.Phi);
%---------------------------------------------
j = j+1;   % (10)
%---------------------------------------------
G0(j,V.y_star) = 1;
G0(j,V.k_s_star) = -para(P.Phi)*para(P.alpha);
G0(j,V.l_star) = -para(P.Phi)*(1-para(P.alpha));
G0(j,V.eps_a) = -para(P.Phi);
%---------------------------------------------
j = j+1;   % (11)
%---------------------------------------------
G0(j,V.k_s) = 1;
G0(j,V.z) = -1;
G1(j,V.k) = 1;
%---------------------------------------------
j = j+1;   % (12)
%---------------------------------------------
G0(j,V.k_s_star) = 1;
G0(j,V.z_star) = -1;
G1(j,V.k_star) = 1;
%---------------------------------------------
j = j+1;   % (13)
%---------------------------------------------
G0(j,V.z) = 1;
G0(j,V.r_k) = -(1-para(P.psi))/para(P.psi);
%---------------------------------------------
j = j+1;   % (14)
%---------------------------------------------
G0(j,V.z_star) = 1;
G0(j,V.r_k_star) = -(1-para(P.psi))/para(P.psi);
%---------------------------------------------
j = j+1;   % (15) 
%---------------------------------------------
G0(j,V.k) = 1;
G0(j,V.i) = -(1-(1-ssp(P.delta))/ssp(P.gamma));
G0(j,V.eps_i) = -(1-(1-ssp(P.delta))/ssp(P.gamma))*para(P.varphi)*ssp(P.gamma)^2*(1+ssp(P.beta)*ssp(P.gamma)^(1-para(P.sigma_c)));
G1(j,V.k) = 1+G0(j,V.i);
%---------------------------------------------
j = j+1;   % (16)
%---------------------------------------------
G0(j,V.k_star) = 1;
G0(j,V.i_star) = -(1-(1-ssp(P.delta))/ssp(P.gamma));
G0(j,V.eps_i) = -(1-(1-ssp(P.delta))/ssp(P.gamma))*para(P.varphi)*ssp(P.gamma)^2*(1+ssp(P.beta)*ssp(P.gamma)^(1-para(P.sigma_c)));
G1(j,V.k_star) = 1+G0(j,V.i_star);
%---------------------------------------------
j = j+1;   % (17) 
%---------------------------------------------
G0(j,V.mu_p) = 1;
G0(j,V.k_s) = -para(P.alpha);
G0(j,V.l) = para(P.alpha);
G0(j,V.w) = 1;
G0(j,V.eps_a) = -1;
%---------------------------------------------
j = j+1;   % (18)
%---------------------------------------------
G0(j,V.w_star) = 1;
G0(j,V.k_s_star) = -para(P.alpha);
G0(j,V.l_star) = para(P.alpha);
G0(j,V.eps_a) = -1;
%---------------------------------------------
j = j+1;   % (19) 
%---------------------------------------------
G0(j,V.pi) = 1;
G0(j,V.E_pi) = -ssp(P.beta)*ssp(P.gamma)^(1-para(P.sigma_c))/(1+para(P.iota_p)*ssp(P.beta)*ssp(P.gamma)^(1-para(P.sigma_c)));
G0(j,V.mu_p) = (1-ssp(P.beta)*ssp(P.gamma)^(1-para(P.sigma_c))*para(P.xi_p))*(1-para(P.xi_p))/(1+para(P.iota_p)*ssp(P.beta)*ssp(P.gamma)^(1-para(P.sigma_c)))/(1+(para(P.Phi)-1)*ssp(P.eps_p))/para(P.xi_p);
G0(j,V.eps_p) = -1;
G1(j,V.pi) = para(P.iota_p)/(1+para(P.iota_p)*ssp(P.beta)*ssp(P.gamma)^(1-para(P.sigma_c)));
%---------------------------------------------
j = j+1;   % (20) 
%---------------------------------------------
G0(j,V.mu_w) = 1;
G0(j,V.w) = -1;
G0(j,V.l) = para(P.sigma_l);
G0(j,V.c) = 1/(1-para(P.h)/ssp(P.gamma));
G1(j,V.c) = para(P.h)/ssp(P.gamma)/(1-para(P.h)/ssp(P.gamma));
%---------------------------------------------
j = j+1;   % (21)
%---------------------------------------------
G0(j,V.r_k) = 1;
G0(j,V.w) = -1;
G0(j,V.k) = 1;
G0(j,V.l) = -1;
%---------------------------------------------
j = j+1;   % (22)
%---------------------------------------------
G0(j,V.r_k_star) = 1;
G0(j,V.w_star) = -1;
G0(j,V.k_star) = 1;
G0(j,V.l_star) = -1;
%---------------------------------------------
j = j+1;   % (23)
%---------------------------------------------
G0(j,V.w) = 1;
G0(j,V.E_w) = -ssp(P.beta)*ssp(P.gamma)^(1-para(P.sigma_c))/(1+ssp(P.beta)*ssp(P.gamma)^(1-para(P.sigma_c)));
G0(j,V.E_pi) = G0(j,V.E_w);
G0(j,V.pi) = (1+ssp(P.beta)*ssp(P.gamma)^(1-para(P.sigma_c))*para(P.iota_w))/(1+ssp(P.beta)*ssp(P.gamma)^(1-para(P.sigma_c)));
G0(j,V.mu_w) = (1-ssp(P.beta)*ssp(P.gamma)^(1-para(P.sigma_c))*para(P.xi_w))*(1-para(P.xi_w))/(1+ssp(P.beta)*ssp(P.gamma)^(1-para(P.sigma_c)))/(1+(ssp(P.lambda_w)-1)*ssp(P.eps_w))/para(P.xi_w);
G0(j,V.eps_w) = -1;
G1(j,V.w) = 1+G0(j,V.E_w);
G1(j,V.pi) = para(P.iota_w)*G1(j,V.w);
%---------------------------------------------
j = j+1;   % (24)
%---------------------------------------------
G0(j,V.w_star) = 1;
G0(j,V.l_star) = -para(P.sigma_l);
G0(j,V.c_star) = -1/(1-para(P.h)/ssp(P.gamma));
G1(j,V.c_star) = -para(P.h)/ssp(P.gamma)/(1-para(P.h)/ssp(P.gamma));
%---------------------------------------------
j = j+1;   % (25)
%---------------------------------------------
G0(j,V.r) = 1;
G0(j,V.pi) = -(1-para(P.rho))*para(P.r_pi);
G0(j,V.y) = -(1-para(P.rho))*para(P.r_y)-para(P.r_delta_y);
G0(j,V.y_star) = (1-para(P.rho))*para(P.r_y)+para(P.r_delta_y);
G0(j,V.eps_r) = -1;
G1(j,V.r) = para(P.rho);
G1(j,V.y) = -para(P.r_delta_y);
G1(j,V.y_star) = para(P.r_delta_y);
%---------------------------------------------
j = j+1;   % (26)
%---------------------------------------------
G0(j,V.eps_a) = 1;
G1(j,V.eps_a) = para(P.rho_a);
Psi(j,V.eta_a) = 1;
%---------------------------------------------
j = j+1;   % (27)
%---------------------------------------------
G0(j,V.eps_b) = 1;
G1(j,V.eps_b) = para(P.rho_b);
Psi(j,V.eta_b) = 1;
%---------------------------------------------
j = j+1;   % (28)
%---------------------------------------------
G0(j,V.eps_g) = 1;
G1(j,V.eps_g) = para(P.rho_g);
Psi(j,V.eta_g) = 1;
Psi(j,V.eta_a) = para(P.rho_ga);
%---------------------------------------------
j = j+1;   % (29)
%---------------------------------------------
G0(j,V.eps_i) = 1;
G1(j,V.eps_i) = para(P.rho_i);
Psi(j,V.eta_i) = 1;
%---------------------------------------------
j = j+1;   % (30)
%---------------------------------------------
G0(j,V.eps_r) = 1;
G1(j,V.eps_r) = para(P.rho_r);
Psi(j,V.eta_r) = 1;
%---------------------------------------------
j = j+1;   % (31)
%---------------------------------------------
G0(j,V.eps_p) = 1;
G0(j,V.d_eta_p) = -1;
G1(j,V.eps_p) = para(P.rho_p);
G1(j,V.d_eta_p) = -para(P.mu_p);
%---------------------------------------------
j = j+1;   % (32)
%---------------------------------------------
G0(j,V.d_eta_p) = 1;
Psi(j,V.eta_p) = 1;
%---------------------------------------------
j = j+1;   % (33)
%---------------------------------------------
G0(j,V.eps_w) = 1;
G0(j,V.d_eta_w) = -1;
G1(j,V.eps_w) = para(P.rho_w);
G1(j,V.d_eta_w) = -para(P.mu_w);
%---------------------------------------------
j = j+1;   % (34)
%---------------------------------------------
G0(j,V.d_eta_w) = 1;
Psi(j,V.eta_w) = 1;
%---------------------------------------------
j = j+1;   % (35)
%---------------------------------------------
G0(j,V.pi) = 1;
G1(j,V.E_pi) = 1;
Pi(j,V.fe_pi) = 1;
%---------------------------------------------
j = j+1;   % (36)
%---------------------------------------------
G0(j,V.c) = 1;
G1(j,V.E_c) = 1;
Pi(j,V.fe_c) = 1;
%---------------------------------------------
j = j+1;   % (37)
%---------------------------------------------
G0(j,V.l) = 1;
G1(j,V.E_l) = 1;
Pi(j,V.fe_l) = 1;
%---------------------------------------------
j = j+1;   % (38)
%---------------------------------------------
G0(j,V.q) = 1;
G1(j,V.E_q) = 1;
Pi(j,V.fe_q) = 1;
%---------------------------------------------
j = j+1;   % (39)
%---------------------------------------------
G0(j,V.r_k) = 1;
G1(j,V.E_r_k) = 1;
Pi(j,V.fe_r_k) = 1;
%---------------------------------------------
j = j+1;   % (40)
%---------------------------------------------
G0(j,V.i) = 1;
G1(j,V.E_i) = 1;
Pi(j,V.fe_i) = 1;
%---------------------------------------------
j = j+1;   % (41)
%---------------------------------------------
G0(j,V.w) = 1;
G1(j,V.E_w) = 1;
Pi(j,V.fe_w) = 1;
%---------------------------------------------
j = j+1;   % (42)
%---------------------------------------------
G0(j,V.c_star) = 1;
G1(j,V.E_c_star) = 1;
Pi(j,V.fe_c_star) = 1;
%---------------------------------------------
j = j+1;   % (43)
%---------------------------------------------
G0(j,V.l_star) = 1;
G1(j,V.E_l_star) = 1;
Pi(j,V.fe_l_star) = 1;
%---------------------------------------------
j = j+1;   % (44)
%---------------------------------------------
G0(j,V.q_star) = 1;
G1(j,V.E_q_star) = 1;
Pi(j,V.fe_q_star) = 1;
%---------------------------------------------
j = j+1;   % (45)
%---------------------------------------------
G0(j,V.r_k_star) = 1;
G1(j,V.E_r_k_star) = 1;
Pi(j,V.fe_r_k_star) = 1;
%---------------------------------------------
j = j+1;   % (46)
%---------------------------------------------
G0(j,V.i_star) = 1;
G1(j,V.E_i_star) = 1;
Pi(j,V.fe_i_star) = 1;
%---------------------------------------------
j = j+1;   % (47)
%---------------------------------------------
G0(j,V.d_y) = 1;
G1(j,V.y) = 1;
%---------------------------------------------
j = j+1;   % (48)
%---------------------------------------------
G0(j,V.d_c) = 1;
G1(j,V.c) = 1;
%---------------------------------------------
j = j+1;   % (49)
%---------------------------------------------
G0(j,V.d_i) = 1;
G1(j,V.i) = 1;
%---------------------------------------------
j = j+1;   % (50)
%---------------------------------------------
G0(j,V.d_w) = 1;
G1(j,V.w) = 1;

%---------------------------------------------
%           Measurement Equations
%---------------------------------------------
Z(V.YGR,V.y) = 1;
Z(V.YGR,V.d_y) = -1;
Z(V.CGR,V.c) = 1;
Z(V.CGR,V.d_c) = -1;
Z(V.IGR,V.i) = 1;
Z(V.IGR,V.d_i) = -1;
Z(V.WGR,V.w) = 1;
Z(V.WGR,V.d_w) = -1;
Z(V.HOUR,V.l) = 1;
Z(V.INF,V.pi) = 1;
Z(V.FFR,V.r) = 1;
D(V.YGR) = para(P.gamma_bar);
D(V.CGR) = para(P.gamma_bar);
D(V.IGR) = para(P.gamma_bar);
D(V.WGR) = para(P.gamma_bar);
D(V.HOUR) = para(P.l_bar);
D(V.INF) = para(P.pi_bar);
D(V.FFR) = ssp(P.r_bar);
% Measurement covariance (comment out if no measurement error)
% Sigma_u(V.YGR) = 0.1712^2;
% Sigma_u(V.CGR) = 0.1371^2;
% Sigma_u(V.IGR) = 0.4513^2;
% Sigma_u(V.WGR) = 0.1128^2;
% Sigma_u(V.HOUR) = 0.5816^2;
% Sigma_u(V.INF) = 0.1230^2;
% Sigma_u(V.FFR) = 0.1661^2;

%% User Input End Here