function [P,V] = ParVar
% Function PARVAR
%
% Purpose:    Construct parameter & variable structures
%
% Format:     [P,V] = ParVar
%
% Input:      none
%
% Output:     P         parameter structure
%             V         variable structure
%
% Written by Fei Tan, Saint Louis University
% Updated: December 15, 2016

user_parvar

%% -------------------------------------------
%             Parameter Structure           
%---------------------------------------------

% Initialization
para = [para;svp];                     %#ok<NODEF> consolidate parameters
P.mod.name = para(:,1);                % parameter name
P.mod.para = cell2mat(para(:,2));      % current parameter value
P.mod.lb = cell2mat(para(:,3));        % parameter lower bound
P.mod.ub = cell2mat(para(:,4));        % parameter upper bound
P.mod.type = para(:,5);                % prior type
P.mod.para1 = cell2mat(para(:,6));     % prior parameter 1
P.mod.para2 = cell2mat(para(:,7));     % prior parameter 2
P.mod.npara = length(P.mod.para);      % number of model parameters
P.mod.nssp = length(ssp);              %#ok<USENS> number of steady state/implied parameters

% Model parameter indices
for k = 1:P.mod.npara
    str = strrep(para{k,1},' ','');
    if strcmp(str,'mod')
        error('mod is reserved name.')
    else
        eval(['P.' str ' = k;']);
    end
end

% Steady state/implied parameter indices
for k = 1:P.mod.nssp
    str = strrep(ssp{k},' ','');
    if strcmp(str,'mod')
        error('mod is reserved name.')
    else
        eval(['P.' str ' = k;']);
    end   
end

% Keep parameters away from bounds
realsmall = 1e-5;
P.mod.lb = P.mod.lb+realsmall;
P.mod.ub = P.mod.ub-realsmall;

%% -------------------------------------------
%             Variable Structure           
%---------------------------------------------

% Initialization
V.mod.var = mvar;                      %#ok<USENS> model variable
V.mod.shock = shock;                   %#ok<USENS> structural shock
V.mod.nvar = length(mvar);             % number of model variables
V.mod.nshock = length(shock);          % number of structural shocks
V.mod.nfore = length(fore);            %#ok<USENS> number of forecast errors
V.mod.ndata = length(data);            %#ok<USENS> number of observables

% Model variable indices
for k = 1:V.mod.nvar
    str = strrep(mvar{k},' ','');
    if strcmp(str,'mod')
        error('mod is reserved name.')
    else
        eval(['V.' str ' = k;']);
    end
end

% Structural shock indices
for k = 1:V.mod.nshock
    str = strrep(shock{k},' ','');
    if strcmp(str,'mod')
        error('mod is reserved name.')
    else
        eval(['V.' str ' = k;']);
    end
end

% Forecast error indices
for k = 1:V.mod.nfore
    str = strrep(fore{k},' ','');
    if strcmp(str,'mod')
        error('mod is reserved name.')
    else
        eval(['V.' str ' = k;']);
    end
end

% Observable indices
for k = 1:V.mod.ndata
    str = strrep(data{k},' ','');
    if strcmp(str,'mod')
        error('mod is reserved name.')
    else
        eval(['V.' str ' = k;']);
    end
end

% Stochastic volatility parameter indices
svp_ss = zeros(V.mod.nshock,1);
svp_ar = zeros(V.mod.nshock,1);
svp_vv = zeros(V.mod.nshock,1);
user_svp
P.mod.svp_ss = svp_ss;
if any(svp_ar) && any(svp_vv)                    % sv parameters
    P.mod.svp_ar = svp_ar;
    P.mod.svp_vv = svp_vv;
    P.mod.svp = [svp_ss;svp_ar;svp_vv];
else
    P.mod.svp_ar = [];
    P.mod.svp_vv = [];
    P.mod.svp = [];
end
P.mod.cvp = setdiff(1:P.mod.npara,P.mod.svp)';   % non-sv parameters
P.mod.ncvp = length(P.mod.cvp);

%-------------------- END --------------------
