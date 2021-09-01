function spec = tarb(fun,spec,varargin)
% Function TARB
%
% Purpose:    DSGE-SVt main function
%
% Format:     spec = tarb(fun,spec,varargin)
%
% Input:      fun       function handle
%                       @tarb_spec
%                       @sample_prior
%                       @sample_post
%                       @marg_lik
%             spec      current TaRB-MH specification (structure)
%             varargin  (optional) string/value pair
%                       'sdof'     - DSGE: shock t degrees of freedom (Inf)
%                       'lamb'     - DSGE: t shock gamma precisions (1's)
%                       'h'        - DSGE: log-dev volatilities (0's)
%                       'pred'     - DSGE: forecasting horizon (0)
%                       'init'     - DSGE: number of initial observations to condition (0)
%                       'var'      - DSGE-VAR (false)
%                       'prior'    - VAR: DSGE prior weight (0.50)
%                       'nlag'     - VAR: number of VAR lags (4)
%                       'tune'     - optimization tuning parameter (0)
%                       'npd'      - number of prior draws (1000)
%                       'nopt'     - number of optimizations (1)
%                       'sa'       - optimization by simulated annealing (sa options)
%                       'cs'       - optimization by csminwel ({grad,crit,nit,nodisp})
%                       'hess'     - numerical Hessian at mode (false)
%                       'modpath'  - model path ('user')
%                       'datpath'  - data path ('user/data.txt')
%                       'savepath' - save path ('user')
%                       'datrange' - remove head (training) & tail (forecasting) sample ([0 0])
%                       'periods'  - number of periods (200)
%                       'draws'    - number of draws including burn-in (11000)
%                       'burn'     - number of burn-in (1000)
%                       'blocks'   - scalar: number of fixed blocks
%                                    cell w/ row vecs: fixed blocks of non-sv parameter indices
%                       'pdof'     - proposal t degrees of freedom (15); pdof>2
%                       'prob'     - blocking probability (0.80)
%                       'freq'     - tailoring frequency (0.20)
%                       'jump'     - number of iterations before mode-jumping (0)
%                       'harm'     - ML by harmonic mean (false)
%                       'par'      - ML by parallel computation (true)
%
% Output:     spec      updated TaRB-MH specification (structure)
%
% Written by Fei Tan, Saint Louis University
% Updated: February 1, 2020

%% -------------------------------------------
%             All-In-One Function
%---------------------------------------------

% One function to run them all
if isequal(fun,@sample_prior) || isequal(fun,@sample_post) || isequal(fun,@marg_lik)
    fun(spec)
    return
elseif ~isequal(fun,@tarb_spec)
    error('Function does not exist.')
end

% Default setting
if isempty(spec)
    spec.sdof = Inf;
    spec.pred = 0;
    spec.init = 0;
    spec.var = false;
    spec.prior = 0.5;
    spec.nlag = 4;
    spec.chi = 0;
    spec.npd = 1000;
    spec.nopt = 1;
    spec.sa = [];
    spec.cs = {[],1e-3,100,true};
    spec.hess = false;
    spec.modpath = 'user';
    spec.datpath = ['user' filesep 'data.txt'];
    spec.savepath = 'user';
    spec.datrange = [0 0];
    spec.T = 200;
    spec.M = 11000;
    spec.N = 1000;
    spec.B = 1;
    spec.pdof = 15;
    spec.prob = 0.8;
    spec.freq = 0.2;
    spec.jump = 0;
    spec.harm = false;
    spec.par = true;
end

% User setting
narg = length(varargin);
for k = 1:2:narg
    switch varargin{k}
        case 'sdof', spec.sdof = varargin{k+1};
        case 'lamb', spec.lamb = varargin{k+1};
        case 'h', spec.h = varargin{k+1};
        case 'pred', spec.pred = varargin{k+1};
        case 'init', spec.init = varargin{k+1};
        case 'var', spec.var = varargin{k+1};
        case 'prior', spec.prior = varargin{k+1};
        case 'nlag', spec.nlag = varargin{k+1};
        case 'tune', spec.chi = varargin{k+1};
        case 'npd', spec.npd = varargin{k+1};
        case 'nopt', spec.nopt = varargin{k+1};
        case 'sa', spec.sa = varargin{k+1};
        case 'cs', spec.cs = varargin{k+1};
        case 'hess', spec.hess = varargin{k+1};
        case 'modpath', spec.modpath = varargin{k+1};
        case 'datpath', spec.datpath = varargin{k+1};
        case 'savepath', spec.savepath = varargin{k+1};
        case 'datrange', spec.datrange = varargin{k+1};
        case 'periods', spec.T = varargin{k+1};
        case 'draws', spec.M = varargin{k+1};
        case 'burn', spec.N = varargin{k+1};
        case 'blocks', spec.B = varargin{k+1};
        case 'pdof', spec.pdof = varargin{k+1};
        case 'prob', spec.prob = varargin{k+1};
        case 'freq', spec.freq = varargin{k+1};
        case 'jump', spec.jump = varargin{k+1};
        case 'harm', spec.harm = varargin{k+1};
        case 'par', spec.par = varargin{k+1};
        otherwise
            error('Unrecognized optional argument.')
    end
end

% Set TaRB files on top of search path
addpath(genpath('utils'),genpath('mex'),spec.modpath);
chk_dir(spec.savepath);
save([spec.savepath filesep 'tarb_spec.mat'],'spec');

%-------------------- END --------------------
