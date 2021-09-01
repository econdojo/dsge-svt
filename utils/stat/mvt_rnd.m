function x = mvt_rnd(mu,Sigma,v,N)
% Function MVT_RND
%
% Purpose:    Generate multivariate Student-t distributed random vectors; note
%             common mean vector and scaling matrix as opposed to MATLAB
%             built-in function MVNRND
%
% Format:     x = mvt_rnd(mu,Sigma,v,N)
%
% Input:      mu        mean vector (1 x dim)
%             Sigma     scaling matrix (dim x dim); covariance = v/(v-2)*Sigma
%             v         degrees of freedom
%             N         number of realizations
%
% Output:     x         parameter row vectors (N x dim)
%
% Written by Fei Tan, Saint Louis University
% Updated: December 15, 2016

%% -------------------------------------------
%     Multivariate Student-t Random Vector
%---------------------------------------------

% Initialization
dim = length(mu);
Mu = repmat(mu,N,1);
R = cholmod(Sigma);

% Choose distribution type
if isinf(v)
    % Multivariate Normal
    x = Mu+randn(N,dim)*R;
else
    % Multivariate Student-t
    x = Mu+repmat(sqrt(v./chi2rnd(v,[N 1])),1,dim).*(randn(N,dim)*R);
end

%-------------------- END --------------------