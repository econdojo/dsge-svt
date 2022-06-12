function sample_post(spec)
% Function SAMPLE_POST
%
% Purpose:    Implement full TaRB-MH MCMC algorithm; see
%             Chib & Ramamurthy (2010), TaRB MCMC methods with application
%             to DSGE models, Journal of Econometrics
%
% Format:     sample_post(spec)
%
% Input:      spec      TaRB-MH specification (structure)
%
% Output:     none
%
% Written by Fei Tan, Saint Louis University
% Updated: July 18, 2017

%% -------------------------------------------
%           Full TaRB-MH Algorithm
%---------------------------------------------

% Find posterior mode under constant volatility
find_mode(spec)
fprintf('Mode search done! Beginning estimation...');
pause(30)

% Check validity of inputs
if spec.M-spec.N<=500
    error('Number of draws after burn-in <= 500.')
elseif spec.pdof<=2
    error('Student-t degrees of freedom <= 2.')
elseif spec.prob<0 || spec.prob>1
    error('Blocking probability <0 or >1.')
elseif spec.freq<0 || spec.freq>1
    error('Tailoring frequency <0 or >1.')
end

% Initialization
load([spec.savepath filesep 'chain_init.mat'],'xopt','Hopt','pid')   % load post mode & (-)inverse Hessian
[P,V] = ParVar;                             % load parameters & variables
data.Y = importdata(spec.datpath);          % import data
data.Y = data.Y((spec.datrange(1)+1):(end-spec.datrange(2)),:);
T = size(data.Y,1);                         % time span
sdof = spec.sdof; spec.sdof = Inf;          % shock t degrees of freedom
sv = ~isempty(P.mod.svp); spec.sv = false;  % stochastic volatility
chain_para = zeros(spec.M,P.mod.npara);     % parameter markov chain
P.mod.para(pid) = xopt';                    % always most current draw
rej1 = 0;                                   % number of overall rejections
rej2 = 0;                                   % number of recycle rejections
nb1 = 0;                                    % number of overall blocks
nb2 = 1e-3;                                 % number of recycle blocks
next = 1;                                   % next round of randomization

% Model-specific setup
if strcmp(spec.mod,'dsge')
    % DSGE t shock precision
    if ~isinf(sdof)
        chain_lamb = zeros(spec.M,T);
    end
    if ~isfield(spec,'lamb') || isinf(sdof)
        spec.lamb = ones(1,T);
    end
    
    % DSGE stochastic volatility
    if sv
        chain_h = zeros(V.mod.nshock,T,spec.M);  % log volatilities
        spec.homo = false;
        [p10,m10,v10] = mix10;
        w = zeros(10,1);
        mixid = repmat(randsample(10,T,true,p10)',V.mod.nshock,1);
    end
    if ~isfield(spec,'h') && sv
        spec.h = repmat(P.mod.para(P.mod.svp_ss),1,T);
    end
else
    % Construct VAR/BMM data matrices
    data = DatMat(data.Y,spec.mod,spec.nlag);
end

% Implement full TaRB-MH MCMC algorithm
tic;
progressbar('TaRB-MH MCMC in Progress')
pk_last = -PostKer1(P.mod.para(P.mod.cvp)',P.mod.cvp,0,spec,P,V,data);    % always most current pk

for iter = 1:spec.M
    % Randomize blocks & recycle
    reopt = ~mod(iter,spec.jump);
    if reopt || iter==next %|| rej1/nb1>0.60
        tailor = 1;                         % tailor proposal density/recycle
        if reopt                            % mode-jumping
            block = RandBlock(P.mod.cvp,1);
            n = 1;
        else                                % random blocks of parameters
            block = RandBlock(P.mod.cvp,spec.prob);
            n = round(round(1/spec.freq-1,0)*2*rand,0)+1;
        end
        nb = length(block);                 % number of blocks within iteration
        xrec = cell(nb,1);                  % recycle theta
        Hrec = cell(nb,1);                  % recycle (-)inverse Hessian
        next = iter+n;                      % recycle for next n iterations
    else
        nb2 = nb2+nb;
    end
    nb1 = nb1+nb;
    
    % Sample blocks of non-sv parameters
    for j = 1:nb
        % Tailored proposal density
        if tailor
            mu = P.mod.para(block{j})';
            Sigma = Hopt(block{j},block{j});
            [xrec{j},Hrec{j},~] = TailorProp(@PostKer1,mu,Sigma,...
                [],spec.cs,block{j},0,spec,P,V,data);
        end
        
        % Metropolis-Hastings step
        theta = mvt_rnd_mex(xrec{j},Hrec{j},spec.pdof,1);  % candidate draw from Student-t
        [~,pk_next,alpha] = MoveProb(@PostKer1,P.mod.para(block{j})',theta,xrec{j},Hrec{j},...
            spec.pdof,pk_last,[],block{j},0,spec,P,V,data);
        if rand > alpha                     % reject
            rej1 = rej1+1;
            if ~tailor
                rej2 = rej2+1;
            end
        else                                % accept
            P.mod.para(block{j}) = theta';
            pk_last = pk_next;
            if reopt                        % update (-)inverse Hessian
                Hopt = Hrec{j};
            end
        end
    end
    
    % Update if DSGE with t shock or sv
    if strcmp(spec.mod,'dsge') && (~isinf(sdof) || sv)
        % Sample shocks & Gamma precisions
        [shock,spec.lamb] = PostShock(sdof,spec,P,V,data.Y);
        
        % Student-t shock
        if ~isinf(sdof)
            chain_lamb(iter,:) = spec.lamb;
        end
        
        % Stochastic volatility by KSC's method
        if sv
            data2 = log(shock.^2.*repmat(spec.lamb',1,V.mod.nshock)+1e-5);
            for k = 1:V.mod.nshock
                pid = [P.mod.svp_ss(k) P.mod.svp_ar(k) P.mod.svp_vv(k)];
                
                % Sample mixture indices
                for t = 1:T
                    for i = 1:10
                        w(i) = exp(mvt_pdf_mex(data2(t,k),spec.h(k,t)+m10(i),v10(i),Inf))*p10(i);
                    end
                    mixid(k,t) = randsample(10,1,true,w);
                end
                
                % Sample sv parameters (integration sampler)
                [xopt2,Hopt2,~] = TailorProp(@PostKer2,P.mod.para(pid)',eye(3),...
                    [],spec.cs,pid,P.mod,data2(:,k),mixid(k,:));
                theta = mvt_rnd_mex(xopt2,Hopt2,spec.pdof,1);
                [~,~,alpha] = MoveProb(@PostKer2,P.mod.para(pid)',theta,xopt2,Hopt2,...
                    spec.pdof,[],[],pid,P.mod,data2(:,k),mixid(k,:));
                if rand <= alpha
                    P.mod.para(pid) = theta';
                end

                % Sample volatilities
                spec.h(k,:) = PostVol(P.mod.para(pid),data2(:,k),mixid(k,:));
            end
            chain_h(:,:,iter) = spec.h;
        end
        
        % Update conditional posterior kernel
        pk_last = -PostKer1(P.mod.para(P.mod.cvp)',P.mod.cvp,0,spec,P,V,data);
    end

    % Update
    chain_para(iter,:) = P.mod.para';
    tailor = 0;
    progressbar(iter/spec.M)
    clc
end

% Display time elapsed
time = datestr(toc/(24*60*60),'DD:HH:MM:SS');
diary([spec.savepath filesep 'mylog.out']); % save screen
fprintf('Excuted on %s\n\n',char(datetime));
switch spec.mod
    case 'dsge'
        fprintf('***** Model = DSGE, sdof = %.1f, sv = %s *****\n\n',sdof,string(sv));
    case 'var'
        fprintf('***** Model = VAR, prior = %.2f, nlag = %d *****\n\n',spec.prior,spec.nlag);
    case 'bmm'
        fprintf('***** Model = BMM, mdof = %.1f, moments = %d *****\n\n',spec.mdof,length(data.M));
end
fprintf('Number of draws = %d after %d burn-in\n',spec.M-spec.N,spec.N);
fprintf('Elapsed time = %s [dd:hh:mm:ss]\n\n',time);

%% -------------------------------------------
%             Posterior Analysis
%---------------------------------------------

% Discard burn-in's
M = spec.M-spec.N;
chain_para = chain_para((spec.N+1):end,:);
save([spec.savepath filesep 'tarb_full.mat'],'chain_para');
if strcmp(spec.mod,'dsge')
    if ~isinf(sdof)
        chain_lamb = chain_lamb((spec.N+1):end,:);
        save([spec.savepath filesep 'tarb_full.mat'],'chain_lamb','-append');
    end
    if sv
        chain_h = chain_h(:,:,(spec.N+1):end);
        save([spec.savepath filesep 'tarb_full.mat'],'chain_h','-append');
    end
end

% Compute posterior statistics
stat = PostStat(chain_para,P.mod.name);
fprintf('\nOverall rejection rate (non-sv)  =  %.1f %%\n',rej1/nb1*100);
fprintf('Recycle rejection rate (non-sv)  =  %.1f %%\n',rej2/nb2*100);
fprintf('Average number of blocks (non-sv)  =  %.1f\n',nb1/spec.M);
fprintf('Average inefficiency factor (all)  =  %.1f\n',mean(stat(:,4)));
save([spec.savepath filesep 'tarb_full.mat'],'stat','-append');

% Convergence diagnostic
DiagPlot(chain_para,P.mod.name,spec.savepath)

% Sample posterior VAR parameters
if strcmp(spec.mod,'var')
    PostVAR(spec.prior,spec.nlag,data,spec.savepath)
end

% Sample predictive data
if strcmp(spec.mod,'dsge') && spec.pred>0
    PostPred(spec.pred,sdof,data.Y,spec.savepath)
end

% Marginal likelihood by harmonic mean
if spec.harm
    pk = zeros(M,1);
    spec.sdof = sdof; spec.sv = sv;         % activate particle filter
    spec.lamb = ones(1,T);
    progressbar('Computing posterior kernel at MCMC sample...')
    for k = 1:M
        pk(k) = -PostKer1(chain_para(k,:),1:P.mod.npara,0,spec,P,V,data);
        progressbar(k/M)
    end
    save([spec.savepath filesep 'tarb_full.mat'],'pk','-append');
    
    [logmlik,nse] = MarLik('harm',spec.savepath);
    fprintf('Truncation   Log marginal likelihood\n');
    p = 0.1:0.2:0.9;
    for k = 1:5
        fprintf('p = %.2f:    %.2f (n.s.e. %.4f)\n',p(k),logmlik(end,k),nse(k));
    end
end
fprintf('\n');
diary off

%-------------------- END --------------------
