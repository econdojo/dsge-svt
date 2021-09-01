% DSGE-SVt demo

%% Housekeeping
clear
close all
clc
readme

%% User search path & mex files
model = 'ltw17';
regime = 'reg_M';
modpath = ['user' filesep model filesep regime];
datpath = ['user' filesep model filesep 'data' filesep 'data_551to142.txt'];
savepath = ['user' filesep model filesep regime];
ltw17 = tarb(@tarb_spec,[],'modpath',modpath,'datpath',datpath,'savepath',savepath);  % TaRB specification
%OneFileToMexThemAll

%% Sample prior
fprintf('Step 1: sample prior. Press any key to continue...\n');
pause;

sdof = 5;          % shock degrees of freedom
npd = 1000;        % number of prior draws
T = 200;           % number of periods
ltw17 = tarb(@tarb_spec,ltw17,'sdof',sdof,'npd',npd,'periods',T);   % update specification
tarb(@sample_prior,ltw17);

%% MCMC (full run)
fprintf('Step 2: sample posterior. Press any key to continue...\n');
pause;

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
head = 50;         % training sample size
tail = 22;         % forecasting sample size
h = 1;             % forecasting horizon
ltw17 = tarb(@tarb_spec,ltw17,'sa',sa_spec,'cs',cs_spec,'npd',npd,'nopt',nopt,'prob',p,'freq',w,'draws',M,'burn',N,'datrange',[head tail],'pred',h);
tarb(@sample_post,ltw17);

%% Marginal likelihood (reduced run)
fprintf('Step 3: marginal likelihood. Press any key to continue...\n');
pause;

B = 7;             % number of blocks for structural parameters
ltw17 = tarb(@tarb_spec,ltw17,'blocks',B);
tarb(@marg_lik,ltw17);
