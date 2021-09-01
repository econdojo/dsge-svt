function logpdf = mvt_pdf(x,mu,Sigma,v)
% Function MVT_PDF
%
% Purpose:    Evaluate log multivariate Student-t density; note
%             common mean vector and scaling matrix as opposed to MATLAB
%             built-in function MVNPDF
%
% Format:     logpdf = mvt_pdf(x,mu,Sigma,v)
%
% Input:      x         parameter row vectors (N x dim)
%             mu        mean vector (1 x dim)
%             Sigma     scaling matrix (dim x dim); covariance = v/(v-2)*Sigma
%             v         degrees of freedom
%
% Output:     logpdf    log multivariate Student-t density (N x 1)
%
% Written by Fei Tan, Saint Louis University
% Updated: December 15, 2016

%% -------------------------------------------
%        Multivariate Student-t Density
%---------------------------------------------

% Initialization
[N,dim] = size(x);
Mu = repmat(mu,N,1);
R = cholmod(Sigma);

% Choose density type
if isinf(v)
    % Multivariate Normal
    const = -dim/2*log(2*pi)-sum(log(diag(R)));
    logpdf = const-sum(((x-Mu)/R).^2,2)/2;
else
    % Multivariate Student-t
    const = v/2*log(v)+log(gamma((dim+v)/2))-dim/2*log(pi)-log(gamma(v/2))-sum(log(diag(R)));
    logpdf = const-(dim+v)/2*log(abs(v+sum(((x-Mu)/R).^2,2)));
end

%-------------------- END --------------------
