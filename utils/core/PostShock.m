function [shock,lamb] = PostShock(sdof,spec,P,V,Y)
% Function POSTSHOCK
%
% Purpose:    Draw posterior shocks & gamma precisions
%
% Format:     [shock,lamb] = PostShock(sdof,spec,P,V,Y)
%
% Input:      sdof      Student-t shock degrees of freedom
%             spec      model specification
%             P         structure of model parameters
%             V         structure of model variables
%             Y         data with ROW observations
%
% Output:     shock     structural shocks (next)
%             lamb      gamma precisions (next)
%
% Written by Fei Tan, Saint Louis University
% Updated: December 15, 2016

%% -------------------------------------------
%      Sample Shocks & Gamma Precisions
%---------------------------------------------

% Construct state space representation
para = P.mod.para;
G0 = zeros(V.mod.nvar);
G1 = zeros(V.mod.nvar);
Psi = zeros(V.mod.nvar,V.mod.nshock);
Pi = zeros(V.mod.nvar,V.mod.nfore);
CC = zeros(V.mod.nvar,1);
Sigma_u = zeros(V.mod.ndata,1);
Z = zeros(V.mod.ndata,V.mod.nvar);
D = zeros(V.mod.ndata,1);
ssp = zeros(P.mod.nssp,1); %#ok<NASGU>
j = 0; %#ok<NASGU>
user_ssp
user_mod
SSR.Sigma_u = Sigma_u;
SSR.Z = Z;
SSR.D = D;
[SSR.G,SSR.C,SSR.M,~,~,~,~,~] = gensys(G0,G1,CC,Psi,Pi);

% Stochastic volatility
T = length(spec.lamb);
SSR.sdof = Inf;
SSR.sv = false;
if isempty(P.mod.svp)
    Sigma_e = para(P.mod.svp_ss).^2;
    SSR.Sigma_e = repmat(Sigma_e,1,T)./repmat(spec.lamb,V.mod.nshock,1);
    SSR.Sigma_v = zeros(V.mod.nshock,1);
    SSR.B = zeros(V.mod.nshock);
    SSR.A = zeros(V.mod.nshock,1);
    spec.h = repmat(log(Sigma_e),1,T);
else
    Sigma_e = exp(para(P.mod.svp_ss));
    SSR.Sigma_e = exp(spec.h)./repmat(spec.lamb,V.mod.nshock,1);
    SSR.Sigma_v = para(P.mod.svp_vv);
    SSR.B = diag(para(P.mod.svp_ar));
    SSR.A = (eye(V.mod.nshock)-SSR.B)*para(P.mod.svp_ss);
end

% Durbin-Koopman simulation smoother
state = (eye(V.mod.nvar)-SSR.G)\SSR.C;
log_vol = (eye(V.mod.nshock)-SSR.B)\SSR.A;
[YY,shock] = SimuData_mex(SSR,T,state,log_vol);
SSR.C = zeros(V.mod.nvar,1);
SSR.D = zeros(V.mod.ndata,1);
fs = zeros(V.mod.nvar,1);
Omega_fs = dlyap_sym(SSR.G,SSR.M*diag(Sigma_e)*SSR.M',1);
se = DistSmoother_mex(Y-YY,SSR,fs,Omega_fs);
shock = se'+shock;

% Sample Gamma precisions
lamb = spec.lamb;
if ~isinf(sdof)
    for t = 1:T
        lamb(t) = gamrnd((sdof+V.mod.nshock)/2,2/(shock(t,:)/diag(exp(spec.h(:,t)))*shock(t,:)'+sdof));
    end
end

%-------------------- END --------------------
