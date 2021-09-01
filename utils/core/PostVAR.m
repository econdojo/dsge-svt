function PostVAR(prior,nlag,data,savepath)
% Function POSTVAR
%
% Purpose:    Sample posterior VAR parameters; see
%             Del Negro & Schorfheide (2004), Priors from general equilibrium
%             models for VARs, International Economic Review
%
% Format:     PostVAR(prior,nlag,data,savepath)
%
% Input:      prior     DSGE prior weight
%             nlag      number of VAR lags
%             data      structure with VAR data matrices
%             savepath  save path
%
% Output:     none
%
% Written by Fei Tan, Saint Louis University
% Updated: July 18, 2017

%% -------------------------------------------
%            Sample VAR Posterior
%---------------------------------------------

% Initialization
load([savepath filesep 'tarb_full.mat'],'chain_para') % load MCMC draws
[P,V] = ParVar;                   % load parameters & variables
M = size(chain_para,1);           % number of post draws
[T,n] = size(data.Y);             % time span & number of data variables
k = 1+n*nlag;                     % number of regressors
Phi = zeros(k,n,M);               % VAR matrix coefficients
Sigma_u = zeros(n,n,M);           % VAR innovation matrix

% Sample posterior VAR parameters
progressbar('VAR Posterior Sampling in Progress')

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
    SSR.Sigma_e = para(P.mod.svp_ss).^2;
    SSR.Sigma_u = Sigma_u;
    SSR.Z = Z;
    SSR.D = D;
    [SSR.G,SSR.C,SSR.M,~,~,~,~,~] = gensys(G0,G1,CC,Psi,Pi);
    
    % Assemble DSGE population moments
    [E_y,E_yy_k] = Moment(SSR,nlag);
    E_yx = [E_y reshape(E_yy_k(:,:,2:end),n,n*nlag)];
    E_xx = zeros(n*nlag);
    for i = 1:nlag
        for j = i:nlag
            E_xx(((i-1)*n+1):(i*n),((j-1)*n+1):(j*n)) = E_yy_k(:,:,j-i+1);
        end
    end
    E_xx = triu(E_xx)+triu(E_xx,1)';
    E_xx = [1 repmat(E_y',1,nlag);repmat(E_y,nlag,1) E_xx]; %#ok<AGROW>
    
    % VAR posterior
    if isinf(prior)               % VAR subject to DSGE restrictions
        Phi(:,:,iter) = E_xx\E_yx';
        Sigma_u(:,:,iter) = E_yy_k(:,:,1)-E_yx*Phi(:,:,iter);
    else
        Phi_post = (prior*T*E_xx+data.X'*data.X)\(prior*T*E_yx'+data.X'*data.Y);
        E_uu_post = (prior*T*E_yy_k(:,:,1)+data.Y'*data.Y-(prior*T*E_yx+data.Y'*data.X)*Phi_post)/((prior+1)*T);
        R = cholmod(E_uu_post); E_uu_post = R'*R;
        Sigma_u(:,:,iter) = iwishrnd((prior+1)*T*E_uu_post,(prior+1)*T-k);
        Phi(:,:,iter) = reshape(mvt_rnd(reshape(Phi_post,k*n,1)',kron(Sigma_u(:,:,iter),inv(prior*T*E_xx+data.X'*data.X)),Inf,1)',k,n);
    end
    progressbar(iter/M)
end

save([savepath filesep 'tarb_full.mat'],'Phi','Sigma_u','-append');

%-------------------- END --------------------
