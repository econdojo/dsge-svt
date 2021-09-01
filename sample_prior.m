function sample_prior(spec)
% Function SAMPLE_PRIOR
%
% Purpose:    Prior predictive analysis; see
%             Chib & Ergashev (2009), Analysis of multifactor affine yield
%             curve models, Journal of American Statistical Association
%
% Format:     sample_prior(spec)
%
% Input:      spec      TaRB-MH specification (structure)
%
% Output:     none
%
% Written by Fei Tan, Saint Louis University
% Updated: June 3, 2020

%% -------------------------------------------
%               Sampling Prior
%---------------------------------------------

% Initialization
[P,V] = ParVar;                        % load parameters & variables
prior_para = zeros(spec.npd,P.mod.npara);   % parameter
Y_m = zeros(spec.npd,V.mod.ndata);     % data mean
Y_sd = zeros(spec.npd,V.mod.ndata);    % data standard deviation
Y_ar = zeros(spec.npd,V.mod.ndata);    % data autocorrelation
iter1 = 0;
iter2 = 0;

% Sample prior
progressbar('Prior Sampling in Progress')

while iter1<spec.npd
    % Sample parameters
    para = zeros(P.mod.npara,1);
    for k = 1:P.mod.npara
        para(k) = prior_rnd(P.mod.para1(k),P.mod.para2(k),P.mod.type{k});
    end
    iter2 = iter2+1;
    
    % Construct state space representation
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
    if any(any(isinf([G0 G1]))) || any(any(isnan([G0 G1])))
        continue
    else
        [SSR.G,SSR.C,SSR.M,~,~,~,~,eu] = gensys(G0,G1,CC,Psi,Pi);
        if any(eu-1)
            continue
        else       % determinacy
            iter1 = iter1+1;
            prior_para(iter1,:) = para';
            progressbar(iter1/spec.npd)
        end
    end
    
    % Stochastic volatility
    SSR.sdof = spec.sdof;
    SSR.sv = ~isempty(P.mod.svp);
    if isempty(P.mod.svp)
        Sigma_e = para(P.mod.svp_ss).^2;
        SSR.Sigma_v = zeros(V.mod.nshock,1);
        SSR.B = zeros(V.mod.nshock);
        SSR.A = zeros(V.mod.nshock,1);
    else
        Sigma_e = exp(para(P.mod.svp_ss));
        SSR.Sigma_v = para(P.mod.svp_vv);
        SSR.B = diag(para(P.mod.svp_ar));
        SSR.A = (eye(V.mod.nshock)-SSR.B)*para(P.mod.svp_ss);
    end
    SSR.Sigma_e = repmat(Sigma_e,1,spec.T);

    % Simulate data
    [Y,~] = SimuData(SSR,spec.T);
    Y_m(iter1,:) = mean(Y);
    Y_sd(iter1,:) = std(Y);
    for k = 1:V.mod.ndata
        acf = autocorr(Y(:,k),'NumLags',1);
        Y_ar(iter1,k) = acf(2);
    end
end

c = spec.npd/iter2;     % normalization constant for truncated prior
save([spec.savepath filesep 'tarb_prior.mat'],'prior_para','Y_m','Y_sd','Y_ar','c');

%-------------------- END --------------------
