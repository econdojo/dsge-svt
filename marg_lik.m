function marg_lik(spec)
% Function MARG_LIK
%
% Purpose:    Estimate marginal likelihood with reduced TaRB-MH; see
%             Chib & Jeliazkov (2001), Marginal likelihood from Metropolis
%             Hastings output, Journal of American Statistical Association
%
% Format:     marg_lik(spec)
%
% Input:      spec      TaRB-MH specification (structure)
%
% Output:     none
%
% Written by Fei Tan, Saint Louis University
% Updated: October 20, 2017

%% -------------------------------------------
%          Reduced TaRB-MH Algorithm
%---------------------------------------------

% Find posterior mode under estimated volatility
load([spec.savepath filesep 'tarb_full.mat']) %#ok<LOAD>
if exist('chain_lamb','var')           % t shock: gamma precisions
    spec = tarb(@tarb_spec,spec,'lamb',mean(chain_lamb));
    clearvars chain_lamb
end
if exist('chain_h','var')              % sv: log volatilities
    spec = tarb(@tarb_spec,spec,'h',mean(chain_h,3));
    clearvars chain_h
end
find_mode(spec)
fprintf('Mode search done! ');

% Initialization
load([spec.savepath filesep 'chain_init.mat'],'xopt','pid')% load post mode
[P,V] = ParVar;                        % load parameters & variables
data.Y = importdata(spec.datpath);     % import data
data.Y = data.Y((spec.datrange(1)+1):(end-spec.datrange(2)),:);
T = size(data.Y,1);                    % time span
spec.sv = ~isempty(P.mod.svp);         % stochastic volatility
P.mod.para(pid) = xopt';               % point of evaluation (theta*)
P.mod.para(P.mod.svp) = stat(P.mod.svp,1);

% Check validity of inputs
if spec.M-spec.N<500
    error('Number of draws after burn-in < 500.')
elseif spec.pdof<=2
    error('Student-t degrees of freedom <= 2.')
elseif spec.prob<0 || spec.prob>1
    error('Blocking probability <0 or >1.')
elseif spec.freq<0 || spec.freq>1
    error('Tailoring frequency <0 or >1.')
elseif iscell(spec.B)                  % user-specified blocking
    block = spec.B;
    nb = length(block);
    pid = cat(2,block{:});
    if length(pid)~=P.mod.ncvp || ~isempty(setdiff(P.mod.cvp,pid))
        error('Incorrect blocking of parameter indices.')
    end
else                                   % almost equal-sized blocking
    nb = spec.B;
    block = RandBlock(P.mod.cvp',nb);
end

if P.mod.ncvp/nb>10
    error('Number of blocks is too small.')
end

% Model-specific setup
if strcmp(spec.mod,'dsge')
    if ~isfield(spec,'lamb') || isinf(spec.sdof)
        spec.lamb = ones(1,T);
    end
    if ~isfield(spec,'h') && spec.sv
        spec.h = repmat(P.mod.para(P.mod.svp_ss),1,T);
    end
else
    % Construct VAR/BMM data matrices
    data = DatMat(data.Y,spec.mod,spec.nlag);
end

% Implement reduced TaRB-MH MCMC algorithm
num = zeros(spec.M,nb+spec.sv*V.mod.nshock); % summand in numerator
den1 = zeros(spec.M,nb+spec.sv*V.mod.nshock);% summand in previous block denominator
den2 = zeros(spec.M,nb+spec.sv*V.mod.nshock);% summand in last block denominator
rej1 = zeros(nb+spec.sv*V.mod.nshock,1);     % overall rejection rate
rej2 = zeros(nb+spec.sv*V.mod.nshock,1);     % recycle rejection rate

tic;
dispstat('','init');
dispstat('Begining estimation...','keepthis','timespamp');
if spec.par
    mypool = parpool;                  % open parallel pool
    pctRunOnAll warning off            % suppress warning on all workers
    parfor b = 1:nb+spec.sv*V.mod.nshock
        [num(:,b),den1(:,b),den2(:,b),rej1(b),rej2(b)] = PostOrd(block,b,spec);
    end
    warning on all                     % turn warning on
    delete(mypool);                    % close parallel pool
else
    for b = 1:nb+spec.sv*V.mod.nshock
        [num(:,b),den1(:,b),den2(:,b),rej1(b),rej2(b)] = PostOrd(block,b,spec);
    end
end
dispstat('Ending estimation...','keepprev');
fprintf('\n');

% Display time elapsed
total = datestr(toc/(24*60*60),'DD:HH:MM:SS');
diary([spec.savepath filesep 'mylog.out']); % save screen
fprintf('Excuted on %s\n\n',char(datetime));
switch spec.mod
    case 'dsge'
        fprintf('***** Model = DSGE, sdof = %.1f, sv = %s *****\n\n',spec.sdof,string(spec.sv));
    case 'var'
        fprintf('***** Model = VAR, prior = %.2f, nlag = %d *****\n\n',spec.prior,spec.nlag);
    case 'bmm'
        fprintf('***** Model = BMM, mdof = %.1f, moments = %d *****\n\n',spec.mdof,length(data.M));
end
fprintf('Number of draws = %d after %d burn-in\n',spec.M-spec.N,spec.N);
fprintf('Elapsed time = %s [dd:hh:mm:ss]\n\n',total);

%% -------------------------------------------
%            Marginal Likelihood
%---------------------------------------------

% Discard burn-in's
num = num((spec.N+1):end,:);
den = [den1(:,2:end) den2(:,end)];
den = den((spec.N+1):end,:);

% Display rejection rate & run time
disp('Block No.     Overall     Recycle     Indices')
for k = 1:nb
    fprintf('Block: %d      %.1f %%      %.1f %%      (%s)\n',...
        k,rej1(k)*100,rej2(k)*100,strjoin(strtrim(cellstr(num2str(block{k}'))),','));
end

% Marginal likelihood
ns = 10;
pk = zeros(ns,1);
spec.lamb = ones(1,T);
fprintf('\nComputing posterior kernel at point of evaluation...\n');
for k = 1:ns                      % average out sampling error
    pk(k) = -PostKer1(P.mod.para',1:P.mod.npara,0,spec,P,V,data);
end
pk = mean(pk);
save([spec.savepath filesep 'tarb_reduce.mat'],'num','den','block','pk');

[logmlik,nse] = MarLik('chib',spec.savepath);
fprintf('Log marginal likelihood  =  %.2f (n.s.e. %.4f)\n\n',logmlik(end),nse);
diary off

%-------------------- END --------------------
