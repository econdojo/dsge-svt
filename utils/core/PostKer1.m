function postker = PostKer1(theta,block,chi,spec,P,V,data)
% Function POSTKER1
%
% Purpose:    Evaluate log posterior kernel density
%
% Format:     postker = PostKer1(theta,block,chi,spec,P,V,data)
%
% Input:      theta     block of model parameters
%             block     block of parameter indices
%             chi       optimization tuning parameter
%             spec      model specification
%             P         structure of model parameters
%             V         structure of model variables
%             data      data structure with ROW observations
%
% Output:     postker   (-)log posterior kernel density
%
% Written by Fei Tan, Saint Louis University
% Updated: July 18, 2017

%% -------------------------------------------
%          Evaluate Posterior Kernel
%---------------------------------------------

% Construct state space representation
if ~chi && (any(theta'<P.mod.lb(block)) || any(theta'>P.mod.ub(block)))
    postker = 1e10;
    return
else
    P.mod.para(block) = theta';
    para = ParTran(P.mod.para',P.mod.lb',P.mod.ub',chi,'con')';
    G0 = zeros(V.mod.nvar);
    G1 = zeros(V.mod.nvar);
    Psi = zeros(V.mod.nvar,V.mod.nshock);
    Pi = zeros(V.mod.nvar,V.mod.nfore);
    CC = zeros(V.mod.nvar,1);
    Sigma_u = zeros(V.mod.ndata,1);
    Z = zeros(V.mod.ndata,V.mod.nvar);
    D = zeros(V.mod.ndata,1);
    Y = data.Y; %#ok<NASGU>
    ssp = zeros(P.mod.nssp,1); %#ok<NASGU>
    j = 0; %#ok<NASGU>
    user_ssp
    user_mod
    SSR.Sigma_u = Sigma_u;
    SSR.Z = Z;
    SSR.D = D;
    if any(any(isinf([G0 G1]))) || any(any(isnan([G0 G1])))
        postker = 1e10;
        return
    end
    [SSR.G,SSR.C,SSR.M,~,~,~,~,eu] = gensys(G0,G1,CC,Psi,Pi);
end

% Evaluate posterior kernel
if any(eu-1)
    postker = 1e10;
    return
else
    % Choose model type
    if spec.var
        % DSGE log marginal likelihood
        loglik = LogLik(SSR,spec.prior,spec.nlag,data);
    else
        % Stochastic volatility
        T = length(spec.lamb);
        SSR.sdof = spec.sdof;
        SSR.sv = spec.sv;
        if isempty(P.mod.svp)
            Sigma_e = para(P.mod.svp_ss).^2;
            SSR.Sigma_e = repmat(Sigma_e,1,T)./repmat(spec.lamb,V.mod.nshock,1);
            SSR.Sigma_v = zeros(V.mod.nshock,1);
            SSR.B = zeros(V.mod.nshock);
            SSR.A = zeros(V.mod.nshock,1);
        else
            Sigma_e = exp(para(P.mod.svp_ss));
            if SSR.sv || spec.homo
                SSR.Sigma_e = repmat(Sigma_e,1,T);
            else
                SSR.Sigma_e = exp(spec.h)./repmat(spec.lamb,V.mod.nshock,1);
            end
            SSR.Sigma_v = para(P.mod.svp_vv);
            SSR.B = diag(para(P.mod.svp_ar));
            SSR.A = (eye(V.mod.nshock)-SSR.B)*para(P.mod.svp_ss);
        end
        
        % Conventional/mixture Kalman filter
        fs = (eye(V.mod.nvar)-SSR.G)\SSR.C;
        Omega_fs = dlyap_sym(SSR.G,SSR.M*diag(Sigma_e)*SSR.M',1);
        N = 1000*V.mod.nshock;
        [~,~,loglik] = MixKalman_mex(data.Y,SSR,fs,Omega_fs,N);
        loglik = loglik(spec.init+1:end);    % condition on initial observations
    end
    
    % Compute (-)log posterior kernel density
    logprior = 0;
    if length(block)<=P.mod.ncvp
        pid = P.mod.cvp;
    elseif length(block)<P.mod.npara
        pid = [P.mod.cvp;P.mod.svp_ss];
    else
        pid = 1:P.mod.npara;
    end
    for k = 1:length(pid)
        logprior = logprior+prior_pdf(para(pid(k)),P.mod.para1(pid(k)),P.mod.para2(pid(k)),P.mod.type{pid(k)});
    end
    postker = -(sum(loglik)+logprior);
end

%-------------------- END --------------------
