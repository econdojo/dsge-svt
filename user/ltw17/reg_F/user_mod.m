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
% Updated: August 24, 2017

%% User Input Start Here

%---------------------------------------------
%              GENSYS Equations
%---------------------------------------------

%---------------------------------------------
j = j+1;   % (1) Production function
%---------------------------------------------
G0(j,V.y) = 1;
G0(j,V.k) = -((ssp(P.y)+ssp(P.omeg))/ssp(P.y))*ssp(P.alpha);
G0(j,V.l) = -((ssp(P.y)+ssp(P.omeg))/ssp(P.y))*(1-ssp(P.alpha));
%---------------------------------------------
j = j+1;   % (2) Production factors
%---------------------------------------------
G0(j,V.r_k) = 1;
G0(j,V.w) = -1;
G0(j,V.k) = 1;
G0(j,V.l) = -1;
%---------------------------------------------
j = j+1;   % (3) Marginal cost
%---------------------------------------------
G0(j,V.mc) = 1;
G0(j,V.r_k) = -ssp(P.alpha);
G0(j,V.w) = ssp(P.alpha)-1;
%---------------------------------------------
j = j+1;   % (4) Phillips equation
%---------------------------------------------
G0(j,V.pi) = ssp(P.lamp);
G0(j,V.E_pi) = -ssp(P.lamp)*ssp(P.beta)/(1+ssp(P.beta)*para(P.chi_p));
G1(j,V.pi) = ssp(P.lamp)*para(P.chi_p)/(1+ssp(P.beta)*para(P.chi_p));
G0(j,V.mc) = -1;
G0(j,V.u_p) = -ssp(P.lamp);
%---------------------------------------------
j = j+1;   % (5) Savers' Lagrange multiplier
%---------------------------------------------
G0(j,V.lambda) = 1;
G0(j,V.u_a) = para(P.theta)/(ssp(P.expg)-para(P.theta));
G0(j,V.c_star) = ssp(P.expg)/(ssp(P.expg)-para(P.theta));
G1(j,V.c_star) = para(P.theta)/(ssp(P.expg)-para(P.theta));
G0(j,V.u_b) = -1;
G0(j,V.tau_c) = (ssp(P.tau_c)/(1+ssp(P.tau_c)));
%---------------------------------------------
j = j+1;   % (6) Consumption in utility
%---------------------------------------------
G0(j,V.c_star) = 1;
G0(j,V.c_s) = -ssp(P.c_s)/(ssp(P.c_s)+para(P.alpha_g)*ssp(P.g));
G0(j,V.g) = -para(P.alpha_g)*ssp(P.g)/(ssp(P.c_s)+para(P.alpha_g)*ssp(P.g));
%---------------------------------------------
j = j+1;   % (7) Euler equation
%---------------------------------------------
G0(j,V.lambda) = 1; 
G0(j,V.R) = -1;
G0(j,V.E_pi) = 1;
G0(j,V.E_lambda) = -1;
G0(j,V.u_a) = para(P.rhos_a);
%---------------------------------------------
j = j+1;   % (8) Capacity utilization
%---------------------------------------------
G0(j,V.r_k) = (1-para(P.psi))/para(P.psi);
G0(j,V.v) = -1;
G0(j,V.tau_k) = -((1-para(P.psi))/para(P.psi))*ssp(P.tau_k)/(1-ssp(P.tau_k));
%--------------------------------------------
j = j+1;   % (9) Capital FOC
%--------------------------------------------
G0(j,V.q) = 1;
G0(j,V.R) = 1;
G0(j,V.E_pi) = -1;
G0(j,V.E_q) = -ssp(P.beta)*exp(-ssp(P.gam))*(1-ssp(P.delta));
G0(j,V.E_r_k) = -ssp(P.beta)*exp(-ssp(P.gam))*ssp(P.R_k)*(1-ssp(P.tau_k));
G0(j,V.E_tau_k) = ssp(P.tau_k)*exp(-ssp(P.gam))*ssp(P.beta)*ssp(P.R_k);
%--------------------------------------------
j = j+1;   % (10) Investment FOC
%--------------------------------------------
G0(j,V.q) = -1/((1+ssp(P.beta))*para(P.s)*(exp(2*ssp(P.gam))));
G0(j,V.i) = 1;
G0(j,V.E_i) = -ssp(P.beta)/(1+ssp(P.beta));
G0(j,V.u_a)=(1-ssp(P.beta)*para(P.rhos_a))/(1+ssp(P.beta));
G1(j,V.i) = 1/(1+ssp(P.beta));
G0(j,V.u_i) = -1; % Note: shock normalized by 1/((1+beta)*s*(expg^2))
%--------------------------------------------
j = j+1;   % (11) Effective capital
%--------------------------------------------
G0(j,V.k) = 1;
G0(j,V.v) = -1;
G0(j,V.u_a) = 1;
G1(j,V.k_bar) = 1;
%--------------------------------------------
j = j+1;   % (12) Law of motion for capital
%--------------------------------------------
G0(j,V.k_bar) = 1;
G0(j,V.u_i) = -(1-(1-ssp(P.delta))*exp(-ssp(P.gam)))*(1+ssp(P.beta))*para(P.s)*(exp(2*ssp(P.gam))); % Note: extra terms due to normalization
G0(j,V.i) = -(1-(1-ssp(P.delta))*exp(-ssp(P.gam)));
G0(j,V.u_a) = (1-ssp(P.delta))*exp(-ssp(P.gam));
G1(j,V.k_bar) = (1-ssp(P.delta))*exp(-ssp(P.gam));
%--------------------------------------------
j = j+1;   % (13) Wage equation
%--------------------------------------------
G0(j,V.w) = 1+ssp(P.lamw);
G0(j,V.E_w) = -ssp(P.lamw)*(ssp(P.beta)/(1+ssp(P.beta)));
G0(j,V.pi) = (ssp(P.lamw)*(1+ssp(P.beta)*para(P.chi_w))/(1+ssp(P.beta)));
G0(j,V.E_pi) = -(ssp(P.lamw)*ssp(P.beta)/(1+ssp(P.beta)));
G0(j,V.l) = -para(P.xi);
G0(j,V.lambda) = 1;
G0(j,V.u_a) = ssp(P.lamw)*((1+ssp(P.beta)*para(P.chi_w)-para(P.rhos_a)*ssp(P.beta))/(1+ssp(P.beta)));
G0(j,V.tau_l) = -(ssp(P.tau_l)/(1-ssp(P.tau_l)));
G0(j,V.u_w) = -ssp(P.lamw);
G0(j,V.u_b) = -1;
G1(j,V.w) = (ssp(P.lamw)/(1+ssp(P.beta)));
G1(j,V.pi) = (ssp(P.lamw)*para(P.chi_w)/(1+ssp(P.beta)));
G1(j,V.u_a) = (ssp(P.lamw)*para(P.chi_w)/(1+ssp(P.beta)));
%--------------------------------------------
j = j+1;   % (14) Monetary policy rule
%--------------------------------------------
G0(j,V.R) = 1;
G1(j,V.R) = para(P.rho_r);
G0(j,V.pi) = -(1-para(P.rho_r))*para(P.phi_p);
G0(j,V.y) = -(1-para(P.rho_r))*para(P.phi_y);
G0(j,V.u_m) = -1;
%--------------------------------------------
j = j+1;   % (15) Aggregate resource constraint
%--------------------------------------------
G0(j,V.c) = ssp(P.c);
G0(j,V.i) = ssp(P.inv);
G0(j,V.y) = -ssp(P.y);
G0(j,V.g) = ssp(P.sgy)*ssp(P.y);
G0(j,V.v) = ssp(P.psi1)*ssp(P.k);
%--------------------------------------------
j = j+1;   % (16) Non-savers household's Budget
%--------------------------------------------
G0(j,V.c_n) = ssp(P.c_n)*(1+ssp(P.tau_c));
G0(j,V.tau_c) = ssp(P.tau_c)*ssp(P.c_n);
G0(j,V.w) = -ssp(P.w)*ssp(P.l)*(1-ssp(P.tau_l));
G0(j,V.l) = -ssp(P.w)*ssp(P.l)*(1-ssp(P.tau_l));
G0(j,V.tau_l) = ssp(P.w)*ssp(P.l)*ssp(P.tau_l); 
G0(j,V.z) = -ssp(P.z);
%--------------------------------------------
j = j+1;   % (17) Consumption aggregation
%--------------------------------------------
G0(j,V.c) = ssp(P.c);
G0(j,V.c_s) = -(1-ssp(P.mu_h))*ssp(P.c_s);
G0(j,V.c_n) = -ssp(P.mu_h)*ssp(P.c_n);
%--------------------------------------------
j = j+1;   % (18) Maturity structure
%--------------------------------------------
G0(j,V.R) = 1;
G0(j,V.E_p_b) = -ssp(P.rho)*ssp(P.p_b)/(1+ssp(P.rho)*ssp(P.p_b));
G0(j,V.p_b) = 1;
%--------------------------------------------
j = j+1;   % (19) Government budget constraint
%--------------------------------------------
G0(j,V.b) = ssp(P.sby);
G0(j,V.g) = -ssp(P.sgy);
G0(j,V.z) = -ssp(P.z)/ssp(P.y); 
G0(j,V.tau_k) = ssp(P.tau_k)*ssp(P.r_k)*ssp(P.ky);
G0(j,V.r_k) = ssp(P.tau_k)*ssp(P.r_k)*ssp(P.ky);
G0(j,V.k) = ssp(P.tau_k)*ssp(P.r_k)*ssp(P.ky);
G0(j,V.u_a) = ssp(P.sby)/ssp(P.beta);
G0(j,V.tau_l) = ssp(P.tau_l)*ssp(P.w)*ssp(P.ly);
G0(j,V.w) = ssp(P.tau_l)*ssp(P.w)*ssp(P.ly);
G0(j,V.l) = ssp(P.tau_l)*ssp(P.w)*ssp(P.ly);
G0(j,V.c) = ssp(P.tau_c)*ssp(P.cy);
G0(j,V.tau_c) = ssp(P.tau_c)*ssp(P.cy);
G1(j,V.b) = ssp(P.sby)/ssp(P.beta);
G0(j,V.p_b) = -ssp(P.sby)*ssp(P.rho)/ssp(P.pi_ss)*exp(-ssp(P.gam));
G1(j,V.p_b) = -ssp(P.sby)/ssp(P.beta); 
G0(j,V.pi) = ssp(P.sby)/ssp(P.beta); 
%--------------------------------------------
j = j+1;   % (20) G Rule
%--------------------------------------------
G0(j,V.g) = 1;
G1(j,V.sby) = -(1-para(P.rho_g))*ssp(P.gam_g);
G1(j,V.g) = para(P.rho_g);
G0(j,V.u_g) = -1;
%--------------------------------------------
j = j+1;   % (21) Capital tax rate rule
%--------------------------------------------
G0(j,V.tau_k) = 1;
G1(j,V.sby) = (1-ssp(P.rho_k))*ssp(P.gam_k);
G1(j,V.tau_k) = ssp(P.rho_k);
%--------------------------------------------
j = j+1;   % (22) Labor tax rate rule
%--------------------------------------------
G0(j,V.tau_l) = 1;
G1(j,V.sby) = (1-ssp(P.rho_l))*ssp(P.gam_l);
G1(j,V.tau_l) = ssp(P.rho_l);
%--------------------------------------------
j = j+1;   % (23) Consumption tax rate rule
%--------------------------------------------
G0(j,V.tau_c) = 1;
G1(j,V.tau_c) = ssp(P.rho_c);
%--------------------------------------------
j = j+1;   % (24) Z rule
%--------------------------------------------
G0(j,V.z) = 1;
G1(j,V.sby) = -(1-para(P.rho_z))*ssp(P.gam_z); 
G1(j,V.z) = para(P.rho_z);
G0(j,V.u_z) = -1;
%--------------------------------------------
j = j+1;   % (25) Fisher equation
%--------------------------------------------
G0(j,V.r) = 1;
G0(j,V.R) = -1;
G0(j,V.E_pi) = 1;
%--------------------------------------------
j = j+1;   % (26) sby defined
%--------------------------------------------
G0(j,V.sby) = 1;
G0(j,V.y) = 1;
G0(j,V.b) = -1;
%--------------------------------------------
j = j+1;   % (27) S defined
%--------------------------------------------
G0(j,V.S) = 1;
G0(j,V.tau_k) = -ssp(P.tau_k)*ssp(P.r_k)*ssp(P.k)/ssp(P.S);
G0(j,V.r_k) = -ssp(P.tau_k)*ssp(P.r_k)*ssp(P.k)/ssp(P.S);
G0(j,V.k) = -ssp(P.tau_k)*ssp(P.r_k)*ssp(P.k)/ssp(P.S);
G0(j,V.tau_l) = -ssp(P.tau_l)*ssp(P.w)*ssp(P.l)/ssp(P.S);
G0(j,V.w) = -ssp(P.tau_l)*ssp(P.w)*ssp(P.l)/ssp(P.S);
G0(j,V.l) = -ssp(P.tau_l)*ssp(P.w)*ssp(P.l)/ssp(P.S);
G0(j,V.tau_c) = -ssp(P.tau_c)*ssp(P.c)/ssp(P.S);
G0(j,V.c) = -ssp(P.tau_c)*ssp(P.c)/ssp(P.S);
G0(j,V.z) = ssp(P.z)/ssp(P.S);
G0(j,V.g) = ssp(P.g)/ssp(P.S);
%--------------------------------------------
j = j+1;   % (28) Define E_pi
%--------------------------------------------
G0(j,V.pi) = 1; 
G1(j,V.E_pi) = 1; 
Pi(j,V.f_pi) = 1;
%--------------------------------------------
j = j+1;   % (29) Define E_q
%--------------------------------------------
G0(j,V.q) = 1; 
G1(j,V.E_q) = 1; 
Pi(j,V.f_q) = 1;
%--------------------------------------------
j = j+1;   % (30) Define E_r_k
%--------------------------------------------
G0(j,V.r_k) = 1; 
G1(j,V.E_r_k) = 1; 
Pi(j,V.f_r_k) = 1;
%--------------------------------------------
j = j+1;   % (31) Define E_i
%--------------------------------------------
G0(j,V.i) = 1; 
G1(j,V.E_i) = 1; 
Pi(j,V.f_i) = 1;
%--------------------------------------------
j = j+1;   % (32) Define E_tau_k
%--------------------------------------------
G0(j,V.tau_k) = 1; 
G1(j,V.E_tau_k) = 1; 
Pi(j,V.f_tau_k) = 1;
%--------------------------------------------
j = j+1;   % (33) Define E_w
%--------------------------------------------
G0(j,V.w) = 1; 
G1(j,V.E_w) = 1; 
Pi(j,V.f_w) = 1;
%--------------------------------------------
j = j+1;   % (34) Define E_lambda
%--------------------------------------------
G0(j,V.lambda) = 1; 
G1(j,V.E_lambda) = 1; 
Pi(j,V.f_lambda) = 1;
%--------------------------------------------
j = j+1;   % (35) Define E_p_b
%--------------------------------------------
G0(j,V.p_b) = 1; 
G1(j,V.E_p_b) = 1; 
Pi(j,V.f_p_b) = 1;
%--------------------------------------------
j = j+1;   % (36) Growth rate of technology shock
%--------------------------------------------
G0(j,V.u_a) = 1; 
G1(j,V.u_a) = para(P.rhos_a); 
Psi(j,V.e_a) = 1;
%--------------------------------------------
j = j+1;   % (37) Preference Shock
%--------------------------------------------
G0(j,V.u_b) = 1; 
G1(j,V.u_b) = para(P.rhos_b); 
Psi(j,V.e_b) = 1;
%--------------------------------------------
j = j+1;   % (38) Investment shock
%--------------------------------------------
G0(j,V.u_i) = 1; 
G1(j,V.u_i) = para(P.rhos_i); 
Psi(j,V.e_i) = 1;
%--------------------------------------------
j = j+1;   % (39) Price markup shock
%--------------------------------------------
G0(j,V.u_p) = 1; 
G1(j,V.u_p) = para(P.rhos_p); 
Psi(j,V.e_p) = 1;
%--------------------------------------------
j = j+1;   % (40) Wage markup shock
%--------------------------------------------
G0(j,V.u_w) = 1; 
G1(j,V.u_w) = para(P.rhos_w); 
Psi(j,V.e_w) = 1;
%--------------------------------------------
j = j+1;   % (41) Monetary policy shock
%--------------------------------------------
G0(j,V.u_m) = 1; 
G1(j,V.u_m) = para(P.rhos_m); 
Psi(j,V.e_m) = 1;
%--------------------------------------------
j = j+1;   % (42) G shock
%--------------------------------------------
G0(j,V.u_g) = 1;
G1(j,V.u_g) = para(P.rhos_g); 
Psi(j,V.e_g) = 1;
%--------------------------------------------
j = j+1;   % (43) Z shock
%--------------------------------------------
G0(j,V.u_z)=1; 
G1(j,V.u_z) = para(P.rhos_z); 
Psi(j,V.e_z) = 1;
%---------------------------------------------
j = j+1;   % (44) Lagged c
%---------------------------------------------
G0(j,V.d_c) = 1;
G1(j,V.c) = 1;
%---------------------------------------------
j = j+1;   % (45) Lagged i
%---------------------------------------------
G0(j,V.d_i) = 1;
G1(j,V.i) = 1;
%---------------------------------------------
j = j+1;   % (46) Lagged w
%---------------------------------------------
G0(j,V.d_w) = 1;
G1(j,V.w) = 1;
%---------------------------------------------
j = j+1;   % (47) Lagged g
%---------------------------------------------
G0(j,V.d_g) = 1;
G1(j,V.g) = 1;
%---------------------------------------------
j = j+1;   % (48) Lagged b
%---------------------------------------------
G0(j,V.d_b) = 1;
G1(j,V.b) = 1;

%---------------------------------------------
%           Measurement Equations
%---------------------------------------------
Z(V.cobs,V.c) = 1;
Z(V.cobs,V.u_a) = 1;
Z(V.cobs,V.d_c) = -1;
Z(V.iobs,V.i) = 1;
Z(V.iobs,V.u_a) = 1;
Z(V.iobs,V.d_i) = -1;
Z(V.wobs,V.w) = 1;
Z(V.wobs,V.u_a) = 1;
Z(V.wobs,V.d_w) = -1;
Z(V.gobs,V.g) = 1;
Z(V.gobs,V.u_a) = 1;
Z(V.gobs,V.d_g) = -1;
Z(V.bobs,V.b) = 1;
Z(V.bobs,V.u_a) = 1;
Z(V.bobs,V.d_b) = -1;
Z(V.lobs,V.l) = 1;
Z(V.pobs,V.pi) = 1;
Z(V.robs,V.R) = 1;
D(V.cobs) = para(P.gam100);
D(V.iobs) = para(P.gam100);
D(V.wobs) = para(P.gam100);
D(V.gobs) = para(P.gam100);
D(V.bobs) = para(P.gam100);
D(V.lobs) = para(P.lbar);
D(V.pobs) = para(P.pbar);
D(V.robs) = para(P.pbar)+ssp(P.R_bar);
% Measurement covariance (comment out if no measurement error)
%Sigma_u = (0.2*std(Y)').^2;

%% User Input End Here