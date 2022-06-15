function test_bmm(index) % cluster use

if index==1
    mod = 'dsge';
elseif index==2
    mod = 'bmm';
end

%% User search path & mex files
modpath = ['user' filesep 'sw07'];
datpath = ['user' filesep 'sw07' filesep 'data_661to044.txt'];
savepath = ['user' filesep 'sw07' filesep num2str(index)];
sw07 = tarb(@tarb_spec,[],'modpath',modpath,'datpath',datpath,'savepath',savepath);  % TaRB specification
% OneFileToMexThemAll

%% Sample prior
sdof = Inf;        % shock degrees of freedom
npd = 1000;        % number of prior draws
T = 200;           % number of periods
sw07 = tarb(@tarb_spec,sw07,'sdof',sdof,'npd',npd,'periods',T);   % update specification
tarb(@sample_prior,sw07);

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
p = 0.8;           % blocking probability
w = 0.5;           % tailoring frequency
M = 11000;         % number of draws including burn-in
N = 1000;          % number of burn-in
mdof = 2.1;
nlag = 1;
sw07 = tarb(@tarb_spec,sw07,'sa',sa_spec,'cs',cs_spec,'npd',npd,'nopt',nopt,'prob',p,'freq',w,'draws',M,'burn',N,'mod',mod,'mdof',mdof,'nlag',nlag);
tarb(@sample_post,sw07);

%% Marginal likelihood (reduced run)
B = 6;             % number of blocks for structural parameters
sw07 = tarb(@tarb_spec,sw07,'blocks',B);
tarb(@marg_lik,sw07);