%% Housekeeping
clear
close all
clc

%% User search path & mex files
modpath = ['user' filesep 'as07'];
datpath = ['user' filesep 'as07' filesep 'data.txt'];
savepath = ['user' filesep 'as07'];
as07 = tarb(@tarb_spec,[],'modpath',modpath,'datpath',datpath,'savepath',savepath);  % TaRB specification
%OneFileToMexThemAll

%% Sample prior
sdof = Inf;        % shock degrees of freedom
npd = 1000;        % number of prior draws
T = 200;           % number of periods
as07 = tarb(@tarb_spec,as07,'sdof',sdof,'npd',npd,'periods',T);   % update specification
tarb(@sample_prior,as07);

%% MCMC (full run)
sa_spec = optimoptions(@simulannealbnd,...  % simulated annealing
    'TemperatureFcn',@temperaturefast,...
    'InitialTemperature',2,...
    'TolFun',1e-3,...
    'MaxTime',10,...
    'Display','iter',...
    'DisplayInterval',10);
cs_spec = {[],1e-3,100,true}; % csminwel
npd = 100;         % number of prior draws to initialize mode search
nopt = 2;          % number of mode searches
p = 0.7;           % blocking probability
w = 0.5;           % tailoring frequency
M = 1100;         % number of draws including burn-in
N = 100;          % number of burn-in
mod = 'var';
mdof = Inf;
nlag = 1;
as07 = tarb(@tarb_spec,as07,'sa',sa_spec,'cs',cs_spec,'npd',npd,'nopt',nopt,'prob',p,'freq',w,'draws',M,'burn',N,'mod',mod,'mdof',mdof,'nlag',nlag);
tarb(@sample_post,as07);

%% Marginal likelihood (reduced run)
B = 3;             % number of blocks for structural parameters
as07 = tarb(@tarb_spec,as07,'blocks',B);
tarb(@marg_lik,as07);
