function find_mode(spec)
% Function FIND_MODE
%
% Purpose:    Initialize chain at posterior mode
%
% Format:     find_mode(spec)
%
% Input:      spec      TaRB-MH specification (structure)
%
% Output:     none
%
% Written by Fei Tan, Saint Louis University
% Updated: July 18, 2017

%% -------------------------------------------
%   Find Posterior Mode & (-)Inverse Hessian
%---------------------------------------------

% Initialization
[P,V] = ParVar;                        % load parameters & variables
data.Y = importdata(spec.datpath);     % import data
data.Y = data.Y((spec.datrange(1)+1):(end-spec.datrange(2)),:);
[T,n] = size(data.Y);                  % time span & number of data variables
sdof = spec.sdof; spec.sdof = Inf;     % shock t degrees of freedom
sv = ~isempty(P.mod.svp); spec.sv = false;  % stochastic volatility
pid = P.mod.cvp;                       % parameter indices

% Check validity of inputs
if spec.chi<0
    error('Optimization tuning parameter < 0.')
elseif spec.nopt<1
    error('Number of optimizations < 1.')
elseif sdof<=2
    error('Student-t degrees of freedom <= 2.')
elseif spec.prior<(1+n*spec.nlag+n)/(T-spec.nlag) && strcmp(spec.mod,'var')
    error('DSGE prior weight < lower bound.')
end

% Model-specific setup
if strcmp(spec.mod,'dsge')
    % Optimize under homoscedastic Gaussian shock
    if ~isfield(spec,'lamb') || isinf(sdof)
        spec.lamb = ones(1,T);
    end
    if ~isfield(spec,'h') && sv
        spec.homo = true;
        pid = [P.mod.cvp;P.mod.svp_ss];
    else
        spec.homo = false;
    end
else
    % Construct VAR/BMM data matrices
    data = DatMat(data.Y,spec.mod,spec.nlag);
end

% Re-optimization
xopt = P.mod.para(pid)';
fopt = -PostKer1(xopt,pid,0,spec,P,V,data);
np = length(pid);

for i = 1:spec.nopt
    % Initializer (change via user_parvar.m)
    xini = P.mod.para(pid)';
    fini = -PostKer1(xini,pid,0,spec,P,V,data);
    progressbar(['Prior Sampling in Progress (' num2str(i) '/' num2str(spec.nopt) ')'])
    for j = 1:spec.npd
        xini2 = zeros(1,np);
        for k = 1:np
            xini2(k) = prior_rnd(P.mod.para1(pid(k)),P.mod.para2(pid(k)),P.mod.type{pid(k)});
        end
        fini2 = -PostKer1(xini2,pid,0,spec,P,V,data);
        if fini2>fini
            xini = xini2; fini = fini2;
            disp('Initializer updated.')
        end
        progressbar(j/spec.npd)
    end
    
    % Find mode & (-)inverse Hessian
    xopt2 = ParTran(xini,P.mod.lb(pid)',P.mod.ub(pid)',spec.chi,'unc');
    Hopt2 = eye(np)*1e-4;
    [xopt2,Hopt2,fopt2] = TailorProp(@PostKer1,xopt2,Hopt2,spec.sa,{[],1e-5,1000,false},pid,spec.chi,spec,P,V,data);
    if fopt2>fopt
        xopt = xopt2; Hopt = Hopt2; fopt = fopt2;
    end
end

Jac = diag(JacTran(xopt,P.mod.lb(pid)',P.mod.ub(pid)',spec.chi));
xopt = ParTran(xopt,P.mod.lb(pid)',P.mod.ub(pid)',spec.chi,'con');
Hopt = Jac*Hopt*Jac;
R = cholmod(Hopt); Hopt = R'*R;

% Display optimization output
clc
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
disp('Para             Mode')
for k = 1:np
    fprintf('%s          %.3f\n',P.mod.name{pid(k)},xopt(k));
end
fprintf('\nBest posterior kernel found  =  %.3f\n',fopt);

if spec.hess
    fprintf('Computing Hessian at posterior mode... ');
    [Hopt,~] = hessian(@(theta) PostKer1(theta,pid,0,spec,P,V,data),xopt);
    R = eye(np)/cholmod(Hopt); Hopt = R*R';
    fprintf('Done!\n');
end
fprintf('\n');
diary off
save([spec.savepath filesep 'chain_init.mat'],'fopt','xopt','Hopt','pid');

%-------------------- END --------------------
