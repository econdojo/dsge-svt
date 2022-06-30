function PostPred(h,sdof,Y,savepath)
% Function POSTPRED
%
% Purpose:    Posterior predictive analysis; see
%             Diebold, Schorfheide & Shin (2017), Real-time forecast evaluation
%             of DSGE models with stochastic volatility, Journal of Econometrics
%
% Format:     PostPred(h,sdof,Y,savepath)
%
% Input:      h         forecasting horizon
%             sdof      Student-t shock degrees of freedom
%             Y         data
%             savepath  save path
%
% Output:     none
%
% Written by Fei Tan, Saint Louis University
% Updated: June 3, 2020

%% -------------------------------------------
%           Sample Predictive Data
%---------------------------------------------

% Initialization
load([savepath filesep 'tarb_full.mat'])    %#ok<LOAD> load MCMC draws
[P,V] = ParVar;                             % load parameters & variables
M = size(chain_para,1);                     %#ok<USENS> % number of post draws
[T,n] = size(Y);                            % time span & number of observables
data_h = zeros(h,n,M);                      % predictive data

% Sample predictive data
progressbar('Real-Time Forecasting in Progress')

for iter = 1:M
    % Construct state space representation
    para = chain_para(iter,:)';
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
    if isinf(sdof)
        lamb = ones(1,T);
    else
        lamb = chain_lamb(iter,:); %#ok<USENS>
    end
    
    if isempty(P.mod.svp)
        Sigma_e = para(P.mod.svp_ss).^2;
        SSR.Sigma_e = repmat(Sigma_e,1,T)./repmat(lamb,V.mod.nshock,1);
        SSR.Sigma_v = zeros(V.mod.nshock,1);
        SSR.B = zeros(V.mod.nshock);
        SSR.A = zeros(V.mod.nshock,1);
    else
        log_vol = chain_h(:,:,iter); %#ok<USENS>
        Sigma_e = exp(para(P.mod.svp_ss));
        SSR.Sigma_e = exp(log_vol)./repmat(lamb,V.mod.nshock,1);
        SSR.Sigma_v = para(P.mod.svp_vv);
        SSR.B = diag(para(P.mod.svp_ar));
        SSR.A = (eye(V.mod.nshock)-SSR.B)*para(P.mod.svp_ss);
    end

    % Sample predictive data
    fs = (eye(V.mod.nvar)-SSR.G)\SSR.C;
    Omega_fs = dlyap_sym(SSR.G,SSR.M*diag(Sigma_e)*SSR.M',1);
    [fs,Omega_fs,~] = KalmanFilter(Y,SSR,fs,Omega_fs);
    state = mvt_rnd(fs(:,end)',Omega_fs(:,:,end),Inf,1)';
    SSR.Sigma_e = repmat(Sigma_e,1,h);
    SSR.sdof = sdof;
    SSR.sv = ~isempty(P.mod.svp);
    if isempty(P.mod.svp)
        [data_h(:,:,iter),~] = SimuData(SSR,h,state);
    else
        [data_h(:,:,iter),~] = SimuData(SSR,h,state,log_vol(:,end));
    end
    progressbar(iter/M)
end

save([savepath filesep 'tarb_full.mat'],'data_h','-append');

%-------------------- END --------------------
