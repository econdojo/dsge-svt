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
% Output:     G0,G1,Psi,Pi,CC - GENSYS input
%             Sigma_e   system covariance
%             Sigma_u   measurement covariance
%             Z,D       measurement equation
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
G0(j,V.model.y) = 1;
G0(j,V.model.c) = -ssp(P.c_y);
G0(j,V.model.i) = -ssp(P.i_y);
G0(j,V.model.z) = -ssp(P.z_y);
G0(j,V.model.eps_g) = -1;
%---------------------------------------------
j = j+1;   % (2)
%---------------------------------------------
G0(j,V.model.y_star) = 1;
G0(j,V.model.c_star) = -ssp(P.c_y);
G0(j,V.model.i_star) = -ssp(P.i_y);
G0(j,V.model.z_star) = -ssp(P.z_y);
G0(j,V.model.eps_g) = -1;
%---------------------------------------------
j = j+1;   % (3)
%---------------------------------------------
G0(j,V.model.c) = 1;
G0(j,V.model.E_c) = -1/(1+para(P.h)/ssp(P.gamma));
G0(j,V.model.l) = -ssp(P.wl_c)*(para(P.sigma_c)-1)/para(P.sigma_c)/(1+para(P.h)/ssp(P.gamma));
G0(j,V.model.E_l) = -G0(j,V.model.l);
G0(j,V.model.r) = (1-para(P.h)/ssp(P.gamma))/para(P.sigma_c)/(1+para(P.h)/ssp(P.gamma));
G0(j,V.model.E_pi) = -G0(j,V.model.r);
G0(j,V.model.eps_b) = G0(j,V.model.r);
G1(j,V.model.c) = 1+G0(j,V.model.E_c);
%---------------------------------------------
j = j+1;   % (4)
%---------------------------------------------
G0(j,V.model.c_star) = 1;
G0(j,V.model.E_c_star) = -1/(1+para(P.h)/ssp(P.gamma));
G0(j,V.model.l_star) = -ssp(P.wl_c)*(para(P.sigma_c)-1)/para(P.sigma_c)/(1+para(P.h)/ssp(P.gamma));
G0(j,V.model.E_l_star) = -G0(j,V.model.l_star);
G0(j,V.model.r_star) = (1-para(P.h)/ssp(P.gamma))/para(P.sigma_c)/(1+para(P.h)/ssp(P.gamma));
G0(j,V.model.eps_b) = G0(j,V.model.r_star);
G1(j,V.model.c_star) = 1+G0(j,V.model.E_c_star);
%---------------------------------------------
j = j+1;   % (5)
%---------------------------------------------
G0(j,V.model.i) = 1;
G0(j,V.model.E_i) = -ssp(P.beta)*ssp(P.gamma)^(1-para(P.sigma_c))/(1+ssp(P.beta)*ssp(P.gamma)^(1-para(P.sigma_c)));
G0(j,V.model.q) = -1/para(P.varphi)/ssp(P.gamma)^2/(1+ssp(P.beta)*ssp(P.gamma)^(1-para(P.sigma_c)));
G0(j,V.model.eps_i) = -1;
G1(j,V.model.i) = 1+G0(j,V.model.E_i);
%---------------------------------------------
j = j+1;   % (6)
%---------------------------------------------
G0(j,V.model.i_star) = 1;
G0(j,V.model.E_i_star) = -ssp(P.beta)*ssp(P.gamma)^(1-para(P.sigma_c))/(1+ssp(P.beta)*ssp(P.gamma)^(1-para(P.sigma_c)));
G0(j,V.model.q_star) = -1/para(P.varphi)/ssp(P.gamma)^2/(1+ssp(P.beta)*ssp(P.gamma)^(1-para(P.sigma_c)));
G0(j,V.model.eps_i) = -1;
G1(j,V.model.i_star) = 1+G0(j,V.model.E_i_star);
%---------------------------------------------
j = j+1;   % (7)
%---------------------------------------------
G0(j,V.model.q) = 1;
G0(j,V.model.E_q) = -ssp(P.beta)*(1-ssp(P.delta))/ssp(P.gamma)^para(P.sigma_c);
G0(j,V.model.E_r_k) = -1-G0(j,V.model.E_q);
G0(j,V.model.r) = 1;
G0(j,V.model.E_pi) = -1;
G0(j,V.model.eps_b) = 1;
%---------------------------------------------
j = j+1;   % (8)
%---------------------------------------------
G0(j,V.model.q_star) = 1;
G0(j,V.model.E_q_star) = -ssp(P.beta)*(1-ssp(P.delta))/ssp(P.gamma)^para(P.sigma_c);
G0(j,V.model.E_r_k_star) = -1-G0(j,V.model.E_q_star);
G0(j,V.model.r_star) = 1;
G0(j,V.model.eps_b) = 1;
%---------------------------------------------
j = j+1;   % (9)
%---------------------------------------------
G0(j,V.model.y) = 1;
G0(j,V.model.k_s) = -para(P.Phi)*para(P.alpha);
G0(j,V.model.l) = -para(P.Phi)*(1-para(P.alpha));
G0(j,V.model.eps_a) = -para(P.Phi);
%---------------------------------------------
j = j+1;   % (10)
%---------------------------------------------
G0(j,V.model.y_star) = 1;
G0(j,V.model.k_s_star) = -para(P.Phi)*para(P.alpha);
G0(j,V.model.l_star) = -para(P.Phi)*(1-para(P.alpha));
G0(j,V.model.eps_a) = -para(P.Phi);
%---------------------------------------------
j = j+1;   % (11)
%---------------------------------------------
G0(j,V.model.k_s) = 1;
G0(j,V.model.z) = -1;
G1(j,V.model.k) = 1;
%---------------------------------------------
j = j+1;   % (12)
%---------------------------------------------
G0(j,V.model.k_s_star) = 1;
G0(j,V.model.z_star) = -1;
G1(j,V.model.k_star) = 1;
%---------------------------------------------
j = j+1;   % (13)
%---------------------------------------------
G0(j,V.model.z) = 1;
G0(j,V.model.r_k) = -(1-para(P.psi))/para(P.psi);
%---------------------------------------------
j = j+1;   % (14)
%---------------------------------------------
G0(j,V.model.z_star) = 1;
G0(j,V.model.r_k_star) = -(1-para(P.psi))/para(P.psi);
%---------------------------------------------
j = j+1;   % (15) 
%---------------------------------------------
G0(j,V.model.k) = 1;
G0(j,V.model.i) = -(1-(1-ssp(P.delta))/ssp(P.gamma));
G0(j,V.model.eps_i) = -(1-(1-ssp(P.delta))/ssp(P.gamma))*para(P.varphi)*ssp(P.gamma)^2*(1+ssp(P.beta)*ssp(P.gamma)^(1-para(P.sigma_c)));
G1(j,V.model.k) = 1+G0(j,V.model.i);
%---------------------------------------------
j = j+1;   % (16)
%---------------------------------------------
G0(j,V.model.k_star) = 1;
G0(j,V.model.i_star) = -(1-(1-ssp(P.delta))/ssp(P.gamma));
G0(j,V.model.eps_i) = -(1-(1-ssp(P.delta))/ssp(P.gamma))*para(P.varphi)*ssp(P.gamma)^2*(1+ssp(P.beta)*ssp(P.gamma)^(1-para(P.sigma_c)));
G1(j,V.model.k_star) = 1+G0(j,V.model.i_star);
%---------------------------------------------
j = j+1;   % (17) 
%---------------------------------------------
G0(j,V.model.mu_p) = 1;
G0(j,V.model.k_s) = -para(P.alpha);
G0(j,V.model.l) = para(P.alpha);
G0(j,V.model.w) = 1;
G0(j,V.model.eps_a) = -1;
%---------------------------------------------
j = j+1;   % (18)
%---------------------------------------------
G0(j,V.model.w_star) = 1;
G0(j,V.model.k_s_star) = -para(P.alpha);
G0(j,V.model.l_star) = para(P.alpha);
G0(j,V.model.eps_a) = -1;
%---------------------------------------------
j = j+1;   % (19) 
%---------------------------------------------
G0(j,V.model.pi) = 1;
G0(j,V.model.E_pi) = -ssp(P.beta)*ssp(P.gamma)^(1-para(P.sigma_c))/(1+para(P.iota_p)*ssp(P.beta)*ssp(P.gamma)^(1-para(P.sigma_c)));
G0(j,V.model.mu_p) = (1-ssp(P.beta)*ssp(P.gamma)^(1-para(P.sigma_c))*para(P.xi_p))*(1-para(P.xi_p))/(1+para(P.iota_p)*ssp(P.beta)*ssp(P.gamma)^(1-para(P.sigma_c)))/(1+(para(P.Phi)-1)*ssp(P.eps_p))/para(P.xi_p);
G0(j,V.model.eps_p) = -1;
G1(j,V.model.pi) = para(P.iota_p)/(1+para(P.iota_p)*ssp(P.beta)*ssp(P.gamma)^(1-para(P.sigma_c)));
%---------------------------------------------
j = j+1;   % (20) 
%---------------------------------------------
G0(j,V.model.mu_w) = 1;
G0(j,V.model.w) = -1;
G0(j,V.model.l) = para(P.sigma_l);
G0(j,V.model.c) = 1/(1-para(P.h)/ssp(P.gamma));
G1(j,V.model.c) = para(P.h)/ssp(P.gamma)/(1-para(P.h)/ssp(P.gamma));
%---------------------------------------------
j = j+1;   % (21)
%---------------------------------------------
G0(j,V.model.r_k) = 1;
G0(j,V.model.w) = -1;
G0(j,V.model.k) = 1;
G0(j,V.model.l) = -1;
%---------------------------------------------
j = j+1;   % (22)
%---------------------------------------------
G0(j,V.model.r_k_star) = 1;
G0(j,V.model.w_star) = -1;
G0(j,V.model.k_star) = 1;
G0(j,V.model.l_star) = -1;
%---------------------------------------------
j = j+1;   % (23)
%---------------------------------------------
G0(j,V.model.w) = 1;
G0(j,V.model.E_w) = -ssp(P.beta)*ssp(P.gamma)^(1-para(P.sigma_c))/(1+ssp(P.beta)*ssp(P.gamma)^(1-para(P.sigma_c)));
G0(j,V.model.E_pi) = G0(j,V.model.E_w);
G0(j,V.model.pi) = (1+ssp(P.beta)*ssp(P.gamma)^(1-para(P.sigma_c))*para(P.iota_w))/(1+ssp(P.beta)*ssp(P.gamma)^(1-para(P.sigma_c)));
G0(j,V.model.mu_w) = (1-ssp(P.beta)*ssp(P.gamma)^(1-para(P.sigma_c))*para(P.xi_w))*(1-para(P.xi_w))/(1+ssp(P.beta)*ssp(P.gamma)^(1-para(P.sigma_c)))/(1+(ssp(P.lambda_w)-1)*ssp(P.eps_w))/para(P.xi_w);
G0(j,V.model.eps_w) = -1;
G1(j,V.model.w) = 1+G0(j,V.model.E_w);
G1(j,V.model.pi) = para(P.iota_w)*G1(j,V.model.w);
%---------------------------------------------
j = j+1;   % (24)
%---------------------------------------------
G0(j,V.model.w_star) = 1;
G0(j,V.model.l_star) = -para(P.sigma_l);
G0(j,V.model.c_star) = -1/(1-para(P.h)/ssp(P.gamma));
G1(j,V.model.c_star) = -para(P.h)/ssp(P.gamma)/(1-para(P.h)/ssp(P.gamma));
%---------------------------------------------
j = j+1;   % (25)
%---------------------------------------------
G0(j,V.model.r) = 1;
G0(j,V.model.pi) = -(1-para(P.rho))*para(P.r_pi);
G0(j,V.model.y) = -(1-para(P.rho))*para(P.r_y)-para(P.r_delta_y);
G0(j,V.model.y_star) = (1-para(P.rho))*para(P.r_y)+para(P.r_delta_y);
G0(j,V.model.eps_r) = -1;
G1(j,V.model.r) = para(P.rho);
G1(j,V.model.y) = -para(P.r_delta_y);
G1(j,V.model.y_star) = para(P.r_delta_y);
%---------------------------------------------
j = j+1;   % (26)
%---------------------------------------------
G0(j,V.model.eps_a) = 1;
G1(j,V.model.eps_a) = para(P.rho_a);
Psi(j,V.shock.eta_a) = 1;
%---------------------------------------------
j = j+1;   % (27)
%---------------------------------------------
G0(j,V.model.eps_b) = 1;
G1(j,V.model.eps_b) = para(P.rho_b);
Psi(j,V.shock.eta_b) = 1;
%---------------------------------------------
j = j+1;   % (28)
%---------------------------------------------
G0(j,V.model.eps_g) = 1;
G1(j,V.model.eps_g) = para(P.rho_g);
Psi(j,V.shock.eta_g) = 1;
Psi(j,V.shock.eta_a) = para(P.rho_ga);
%---------------------------------------------
j = j+1;   % (29)
%---------------------------------------------
G0(j,V.model.eps_i) = 1;
G1(j,V.model.eps_i) = para(P.rho_i);
Psi(j,V.shock.eta_i) = 1;
%---------------------------------------------
j = j+1;   % (30)
%---------------------------------------------
G0(j,V.model.eps_r) = 1;
G1(j,V.model.eps_r) = para(P.rho_r);
Psi(j,V.shock.eta_r) = 1;
%---------------------------------------------
j = j+1;   % (31)
%---------------------------------------------
G0(j,V.model.eps_p) = 1;
G0(j,V.model.d_eta_p) = -1;
G1(j,V.model.eps_p) = para(P.rho_p);
G1(j,V.model.d_eta_p) = -para(P.mu_p);
%---------------------------------------------
j = j+1;   % (32)
%---------------------------------------------
G0(j,V.model.d_eta_p) = 1;
Psi(j,V.shock.eta_p) = 1;
%---------------------------------------------
j = j+1;   % (33)
%---------------------------------------------
G0(j,V.model.eps_w) = 1;
G0(j,V.model.d_eta_w) = -1;
G1(j,V.model.eps_w) = para(P.rho_w);
G1(j,V.model.d_eta_w) = -para(P.mu_w);
%---------------------------------------------
j = j+1;   % (34)
%---------------------------------------------
G0(j,V.model.d_eta_w) = 1;
Psi(j,V.shock.eta_w) = 1;
%---------------------------------------------
j = j+1;   % (35)
%---------------------------------------------
G0(j,V.model.pi) = 1;
G1(j,V.model.E_pi) = 1;
Pi(j,V.fore.pi) = 1;
%---------------------------------------------
j = j+1;   % (36)
%---------------------------------------------
G0(j,V.model.c) = 1;
G1(j,V.model.E_c) = 1;
Pi(j,V.fore.c) = 1;
%---------------------------------------------
j = j+1;   % (37)
%---------------------------------------------
G0(j,V.model.l) = 1;
G1(j,V.model.E_l) = 1;
Pi(j,V.fore.l) = 1;
%---------------------------------------------
j = j+1;   % (38)
%---------------------------------------------
G0(j,V.model.q) = 1;
G1(j,V.model.E_q) = 1;
Pi(j,V.fore.q) = 1;
%---------------------------------------------
j = j+1;   % (39)
%---------------------------------------------
G0(j,V.model.r_k) = 1;
G1(j,V.model.E_r_k) = 1;
Pi(j,V.fore.r_k) = 1;
%---------------------------------------------
j = j+1;   % (40)
%---------------------------------------------
G0(j,V.model.i) = 1;
G1(j,V.model.E_i) = 1;
Pi(j,V.fore.i) = 1;
%---------------------------------------------
j = j+1;   % (41)
%---------------------------------------------
G0(j,V.model.w) = 1;
G1(j,V.model.E_w) = 1;
Pi(j,V.fore.w) = 1;
%---------------------------------------------
j = j+1;   % (42)
%---------------------------------------------
G0(j,V.model.c_star) = 1;
G1(j,V.model.E_c_star) = 1;
Pi(j,V.fore.c_star) = 1;
%---------------------------------------------
j = j+1;   % (43)
%---------------------------------------------
G0(j,V.model.l_star) = 1;
G1(j,V.model.E_l_star) = 1;
Pi(j,V.fore.l_star) = 1;
%---------------------------------------------
j = j+1;   % (44)
%---------------------------------------------
G0(j,V.model.q_star) = 1;
G1(j,V.model.E_q_star) = 1;
Pi(j,V.fore.q_star) = 1;
%---------------------------------------------
j = j+1;   % (45)
%---------------------------------------------
G0(j,V.model.r_k_star) = 1;
G1(j,V.model.E_r_k_star) = 1;
Pi(j,V.fore.r_k_star) = 1;
%---------------------------------------------
j = j+1;   % (46)
%---------------------------------------------
G0(j,V.model.i_star) = 1;
G1(j,V.model.E_i_star) = 1;
Pi(j,V.fore.i_star) = 1;
%---------------------------------------------
j = j+1;   % (47)
%---------------------------------------------
G0(j,V.model.d_y) = 1;
G1(j,V.model.y) = 1;
%---------------------------------------------
j = j+1;   % (48)
%---------------------------------------------
G0(j,V.model.d_c) = 1;
G1(j,V.model.c) = 1;
%---------------------------------------------
j = j+1;   % (49)
%---------------------------------------------
G0(j,V.model.d_i) = 1;
G1(j,V.model.i) = 1;
%---------------------------------------------
j = j+1;   % (50)
%---------------------------------------------
G0(j,V.model.d_w) = 1;
G1(j,V.model.w) = 1;

%---------------------------------------------
%           Measurement Equations
%---------------------------------------------
SSR.Z(V.data.YGR,V.model.y) = 1;
SSR.Z(V.data.YGR,V.model.d_y) = -1;
SSR.Z(V.data.CGR,V.model.c) = 1;
SSR.Z(V.data.CGR,V.model.d_c) = -1;
SSR.Z(V.data.IGR,V.model.i) = 1;
SSR.Z(V.data.IGR,V.model.d_i) = -1;
SSR.Z(V.data.WGR,V.model.w) = 1;
SSR.Z(V.data.WGR,V.model.d_w) = -1;
SSR.Z(V.data.HOUR,V.model.l) = 1;
SSR.Z(V.data.INF,V.model.pi) = 1;
SSR.Z(V.data.FFR,V.model.r) = 1;
SSR.D(V.data.YGR) = para(P.gamma_bar);
SSR.D(V.data.CGR) = para(P.gamma_bar);
SSR.D(V.data.IGR) = para(P.gamma_bar);
SSR.D(V.data.WGR) = para(P.gamma_bar);
SSR.D(V.data.HOUR) = para(P.l_bar);
SSR.D(V.data.INF) = para(P.pi_bar);
SSR.D(V.data.FFR) = ssp(P.r_bar);

%---------------------------------------------
%              Covariance Matrix
%---------------------------------------------
% Shock covariance
SSR.Sigma_e(V.shock.eta_a) = para(P.sigma_a)^2;
SSR.Sigma_e(V.shock.eta_b) = para(P.sigma_b)^2;
SSR.Sigma_e(V.shock.eta_g) = para(P.sigma_g)^2;
SSR.Sigma_e(V.shock.eta_i) = para(P.sigma_i)^2;
SSR.Sigma_e(V.shock.eta_r) = para(P.sigma_r)^2;
SSR.Sigma_e(V.shock.eta_p) = para(P.sigma_p)^2;
SSR.Sigma_e(V.shock.eta_w) = para(P.sigma_w)^2;
% Measurement covariance (comment below out if no measurement error)
% SSR.Sigma_u(V.data.YGR) = 0.1712^2;
% SSR.Sigma_u(V.data.CGR) = 0.1371^2;
% SSR.Sigma_u(V.data.IGR) = 0.4513^2;
% SSR.Sigma_u(V.data.WGR) = 0.1128^2;
% SSR.Sigma_u(V.data.HOUR) = 0.5816^2;
% SSR.Sigma_u(V.data.INF) = 0.1230^2;
% SSR.Sigma_u(V.data.FFR) = 0.1661^2;
% SSR.Sigma_u = (0.2*std(data.Y)').^2;

%% User Input End Here