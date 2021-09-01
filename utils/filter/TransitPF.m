function [swarm,w_imp] = TransitPF(y,SSR,swarm,sdof,prop)
% Function TRANSIT
%
% Purpose:    Generate new swarm from transition equation implied by linearized
%             DSGE model solution using approximately conditionally optimal
%             proposal distribution adapted to current observation via one-step
%             Kalman filter; see
%             Herbst & Schorfheide (2015), Bayesian Estimation of DSGE Models
%
% Format:     [swarm,w_imp] = Transit(y,SSR,swarm,sdof,prop)
%
% Input:      y         current observation (row vector)
%             SSR       state space representation (structure)
%
%                       s(t) = Transit(s(t-1),e(t)), E[ee'] = Sigma_e (diagonal)
%
%             swarm     collection of particles (dims x N)
%             sdof      shock t degrees of freedom
%             prop      basic/generic propagation
%
% Output:     swarm     collection of particles
%             w_imp     log importance weights
%
% Written by Fei Tan, Saint Louis University
% Updated: December 15, 2016

%% -------------------------------------------
%          Particle Swarm Propagation
%---------------------------------------------

% Initialization
[dims,N] = size(swarm);                % numbers of states, particles
dime = size(SSR.M,2);                  % number of shocks
dimy = length(y);                      % number of observables
pdof = 15;                             % proposal t degrees of freedom

% Choose propagation method
switch prop
    case 'basic'                       % naive forward propagation
        shock = mvt_rnd(zeros(1,dime),diag(SSR.Sigma_e),sdof,N)';
        w_imp = zeros(N,1);
    case 'generic'                     % generic importance sampler
        % Construct expanded linear SSR
        C = [SSR.C;zeros(dime,1)];
        G = [SSR.G zeros(dims,dime);zeros(dime,dims+dime)];
        M = [SSR.M;eye(dime)];
        D = SSR.D;
        Z = [SSR.Z zeros(dimy,dime)];

        % Period-(t-1) predictive density
        ps = repmat(C,1,N)+G*[swarm;zeros(dime,N)];
        Omega_ps = M*diag(SSR.Sigma_e)*M';

        % Period-t log likelihood
        py = repmat(D,1,N)+Z*ps;
        Omega_py = Z*Omega_ps*Z'+diag(SSR.Sigma_u);

        % Period-t filtering density
        fs = ps+(Omega_ps*Z')/Omega_py*(repmat(y',1,N)-py);
        Omega_fs = Omega_ps-(Omega_ps*Z')/Omega_py*Z*Omega_ps;

        % Generate new swarm & importance weights
        shock = fs((dims+1):end,:)+mvt_rnd(zeros(1,dime),Omega_fs((dims+1):end,(dims+1):end),pdof,N)';
        w_imp = mvt_pdf(shock',zeros(1,dime),diag(SSR.Sigma_e),sdof)...
            -mvt_pdf(shock'-fs((dims+1):end,:)',zeros(1,dime),Omega_fs((dims+1):end,(dims+1):end),pdof);
    otherwise
        error('Propagation does not exist!')
end

% Generate new swarm
swarm = repmat(SSR.C,1,N)+SSR.G*swarm+SSR.M*shock;

%-------------------- END --------------------
