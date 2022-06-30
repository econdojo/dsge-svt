function [num,den1,den2,rej1,rej2] = PostOrd(blk_fix,b,spec)
% Function POSTORD
%
% Purpose:    Implement reduced TaRB-MH MCMC algorithm; see
%             Chib & Jeliazkov (2001), Marginal likelihood from Metropolis
%             Hastings output, Journal of American Statistical Association
%
% Format:     [num,den1,den2,rej1,rej2] = PostOrd(blk_fix,b,spec)
%
% Input:      blk_fix   fixed blocks of non-sv parameter indices
%             b         current block
%             spec      TaRB-MH specification (structure)
%
% Output:     num       summand in current block numerator
%             den1      summand in previous block denominator
%             den2      summand in last block denominator
%             rej1      overall rejection rate
%             rej2      recycle rejection rate
%
% Written by Fei Tan, Saint Louis University
% Updated: October 20, 2017

%% -------------------------------------------
%          Reduced TaRB-MH Algorithm
%---------------------------------------------

% Initialization
load([spec.savepath filesep 'chain_init.mat'],'xopt','Hopt','pid')   % load post mode & (-)inverse Hessian
load([spec.savepath filesep 'tarb_full.mat'],'stat')  % load post mean
[P,V] = ParVar;                             % load parameters & variables
data.Y = importdata(spec.datpath);          % import data
data.Y = data.Y((spec.datrange(1)+1):(end-spec.datrange(2)),:);
sdof = spec.sdof; spec.sdof = Inf;          % shock degrees of freedom
sv = ~isempty(P.mod.svp); spec.sv = false;  % stochastic volatility
T = size(data.Y,1);                         % time span
P.mod.para(pid) = xopt';                    % always most current draw
P.mod.para(P.mod.svp) = stat(P.mod.svp,1);
B = length(blk_fix);                        % number of fixed blocks
ind = cat(2,blk_fix{(b+1):end});            % remaining indices
num = zeros(spec.M,1);                      % summand in current block numerator
den1 = zeros(spec.M,1);                     % summand in previous block denominator
den2 = zeros(spec.M,1);                     % summand in last block denominator
rej1 = 0;                                   % number of overall rejections
rej2 = 0;                                   % number of recycle rejections
nb1 = 0;                                    % number of overall blocks
nb2 = 1e-3;                                 % number of recycle blocks
next = 1;                                   % next round of randomization

% Model-specific setup
if strcmp(spec.mod,'dsge')
    % DSGE stochastic volatility
    if sv
        spec.homo = false;
        [p10,m10,v10] = mix10;
        w = zeros(10,1);
        mixid = repmat(randsample(10,T,true,p10)',V.mod.nshock,1);
    end
else
    % Construct VAR/BMM data matrices
    data = DatMat(data.Y,spec.mod,spec.nlag);
end

% Implement reduced TaRB-MH MCMC algorithm
update = floor((0.2:0.2:1)*spec.M); ii = 1;     % update every 20% progress
pk_last = -PostKer1(P.mod.para(P.mod.cvp)',P.mod.cvp,0,spec,P,V,data);    % always most current pk

for iter = 1:spec.M
    if b<=B
        % Randomize blocks & recycle
        reopt = ~mod(iter,spec.jump);
        if reopt || iter==next %|| rej1(L)/nb1(L)>0.60
            tailor = 1;                         % tailor proposal density/recycle
            if reopt                            % mode-jumping
                blk_rand = RandBlock(ind,1);
                n = 1;
            else                                % randomize remaining indices
                blk_rand = RandBlock(ind,spec.prob);
                n = round((1/spec.freq-1)*2*rand,0)+1;
            end
            blk_rand = cat(2,blk_rand,blk_fix(b));
            nb = length(blk_rand);              % number of blocks within iteration
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
                mu = P.mod.para(blk_rand{j})';
                Sigma = Hopt(blk_rand{j},blk_rand{j});
                [xrec{j},Hrec{j},~] = TailorProp(@PostKer1,mu,Sigma,...
                    [],spec.cs,blk_rand{j},0,spec,P,V,data);
            end

            % Metropolis-Hastings step
            theta = mvt_rnd_mex(xrec{j},Hrec{j},spec.pdof,1);  % candidate draw from Student-t
            [~,pk_next,alpha] = MoveProb(@PostKer1,P.mod.para(blk_rand{j})',theta,xrec{j},Hrec{j},...
                spec.pdof,pk_last,[],blk_rand{j},0,spec,P,V,data);
            if rand > alpha                     % reject
                rej1 = rej1+1;
                if ~tailor
                    rej2 = rej2+1;
                end
            else                                % accept
                P.mod.para(blk_rand{j}) = theta';
                pk_last = pk_next;
                if reopt                        % update (-)inverse Hessian
                    Hopt(blk_rand{j},blk_rand{j}) = Hrec{j};
                end
            end
        end

        % Compute block-b numerator
        [~,pk_opt,alpha] = MoveProb(@PostKer1,P.mod.para(blk_fix{b})',xopt(blk_fix{b}),xrec{end},Hrec{end},...
            spec.pdof,pk_last,[],blk_fix{b},0,spec,P,V,data);
        num(iter) = alpha*exp(mvt_pdf_mex(xopt(blk_fix{b}),xrec{end},Hrec{end},spec.pdof));
    end
    
    % Compute final block denominator
    if b==B && ~sv
        [~,~,den2(iter)] = MoveProb(@PostKer1,xopt(blk_fix{b}),theta,xrec{end},Hrec{end},...
            spec.pdof,pk_opt,pk_next,blk_fix{b},0,spec,P,V,data);
    end
    
    % Compute block-(b-1) denominator
    if b>1 && b<=B+1
        if b==B+1 || tailor
            mu = P.mod.para(blk_fix{b-1})';
            Sigma = Hopt(blk_fix{b-1},blk_fix{b-1});
            [xsup,Hsup,~] = TailorProp(@PostKer1,mu,Sigma,...
                [],spec.cs,blk_fix{b-1},0,spec,P,V,data);
        end
        theta = mvt_rnd_mex(xsup,Hsup,spec.pdof,1);   % supplementary draw
        [~,~,den1(iter)] = MoveProb(@PostKer1,xopt(blk_fix{b-1}),theta,xsup,Hsup,...
            spec.pdof,[],[],blk_fix{b-1},0,spec,P,V,data);
    end
    
    % Update if DSGE with t shock or sv
    if strcmp(spec.mod,'dsge') && (~isinf(sdof) || sv)
        % Sample shocks & Gamma precisions
        [shock,spec.lamb] = PostShock(sdof,spec,P,V,data.Y);
        
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
                if k>=b-B
                    [xopt2,Hopt2,~] = TailorProp(@PostKer2,P.mod.para(pid)',eye(3),...
                        [],spec.cs,pid,P.mod,data2(:,k),mixid(k,:));
                    theta = mvt_rnd_mex(xopt2,Hopt2,spec.pdof,1);
                    [~,~,alpha] = MoveProb(@PostKer2,P.mod.para(pid)',theta,xopt2,Hopt2,...
                        spec.pdof,[],[],pid,P.mod,data2(:,k),mixid(k,:));
                    if rand <= alpha
                        P.mod.para(pid) = theta';
                    end
                end
                
                % Compute block-b numerator
                if k==b-B
                    [~,~,alpha] = MoveProb(@PostKer2,P.mod.para(pid)',stat(pid,1)',xopt2,Hopt2,...
                        spec.pdof,[],[],pid,P.mod,data2(:,k),mixid(k,:));
                    num(iter) = alpha*exp(mvt_pdf_mex(stat(pid,1)',xopt2,Hopt2,spec.pdof));
                end
                
                % Compute final block denominator
                if k==V.mod.nshock && b==B+V.mod.nshock
                    [~,~,den2(iter)] = MoveProb(@PostKer2,stat(pid,1)',theta,xopt2,Hopt2,...
                        spec.pdof,[],[],pid,P.mod,data2(:,k),mixid(k,:));
                end
                
                % Compute block-(b-1) denominator
                if k>1 && k==b-B
                    pid2 = [P.mod.svp_ss(k-1) P.mod.svp_ar(k-1) P.mod.svp_vv(k-1)];
                    [xsup2,Hsup2,~] = TailorProp(@PostKer2,P.mod.para(pid2)',eye(3),...
                        [],spec.cs,pid2,P.mod,data2(:,k-1),mixid(k-1,:));
                    theta = mvt_rnd_mex(xsup2,Hsup2,spec.pdof,1); % supplementary draw
                    [~,~,den1(iter)] = MoveProb(@PostKer2,stat(pid2,1)',theta,xsup2,Hsup2,...
                        spec.pdof,[],[],pid2,P.mod,data2(:,k-1),mixid(k-1,:));
                end
                
                % Sample volatilities
                spec.h(k,:) = PostVol(P.mod.para(pid),data2(:,k),mixid(k,:));
            end
        end
        
        % Update conditional posterior kernel
        if b<=B
            pk_last = -PostKer1(P.mod.para(P.mod.cvp)',P.mod.cvp,0,spec,P,V,data);
        end
    end
    
    % Update
    tailor = 0;
    if iter==update(ii)
        dispstat(sprintf('%d%% of block %d/%d completed.',20*ii,b,B+sv*V.mod.nshock),'keepprev','timestamp');
        ii=ii+1;
    end
end

rej1 = rej1/nb1; rej2 = rej2/nb2;

%-------------------- END --------------------
