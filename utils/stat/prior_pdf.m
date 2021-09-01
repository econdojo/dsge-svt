 function logprior = prior_pdf(x,para1,para2,type)
% Function PRIOR_PDF
%
% Purpose:    Evaluate log prior density
%
% Format:     logprior = prior_pdf(x,para1,para2,type)
%
% Input:      x         parameter value
%             para1     prior parameter 1
%             para2     prior parameter 2
%             type      prior type
%
% Output:     logprior  log prior density
%
% Written by Fei Tan, Saint Louis University
% Updated: December 15, 2016

%% -------------------------------------------
%          Prior Density Evaluation
%---------------------------------------------

% Choose prior type
switch type
    case 't'       % Student-t distribution
        % para1 = mean, para2 = s.d.
        nu = 2.1;  % degrees of freedom
        logprior = mvt_pdf(x,para1,para2^2*(nu-2)/nu,nu);
    case 'G'       % Gamma distribution
        % para1 = mean, para2 = s.d.
        a = para1^2/para2^2;      % shape
        b = para2^2/para1;        % scale
        logprior = log(gampdf(x,a,b));
    case 'N'       % Normal distribution
        % para1 = mean, para2 = s.d.
        logprior = log(normpdf(x,para1,para2));
    case 'LN'      % Lognormal distribution
        % para1 = mean, para2 = s.d.
        a = log(para1^2/sqrt(para2^2+para1^2));  % log-mean
        b = sqrt(log(para2^2/para1^2+1));        % log-s.d.
        logprior = log(lognpdf(x,a,b));
    case 'B'       % Standard Beta distribution (0,1)
        % para1 = mean, para2 = s.d.
        a = -para1*(para2^2+para1^2-para1)/para2^2;   % shape
        b = (para1-1)*(para2^2+para1^2-para1)/para2^2;% shape
        logprior = log(betapdf(x,a,b));
    case 'NB'      % Non-standard Beta distribution (-1,1)
        % para1 = mean, para2 = s.d.
        para1 = (para1+1)/2; % normalize to standard B
        para2 = para2/2;
        a = -para1*(para2^2+para1^2-para1)/para2^2;   % shape
        b = (para1-1)*(para2^2+para1^2-para1)/para2^2;% shape
        logprior = log(betapdf((x+1)/2,a,b)/2);
    case 'I1'      % Inverse-Gamma type-1 distribution
        % para1 = nu, para2 = s; see An & Schorfheide (2007)
        a = para1/2;
        b = a*para2^2;
        logprior = log(2)+a*log(b)-gammaln(a)-(2*a+1)*log(x)-b/x^2;
    case 'I2'      % Inverse-Gamma type-2 distribution
        % para1 = shape, para2 = scale
        logprior = para1*log(para2)-gammaln(para1)-(para1+1)*log(x)-para2/x;
    case 'U'       % Uniform distribution
        % para1 = lb, para2 = ub
        logprior = log(unifpdf(x,para1,para2));
    otherwise
        error('Prior type does not exist!')
end

%-------------------- END --------------------
