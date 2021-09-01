 function x = prior_rnd(para1,para2,type)
% Function PRIOR_RND
%
% Purpose:    Generate random number from prior distribution
%
% Format:     x = prior_rnd(para1,para2,type)
%
% Input:      para1     prior parameter 1
%             para2     prior parameter 2
%             type      prior type
%
% Output:     x         random number
%
% Written by Fei Tan, Saint Louis University
% Updated:    April 30, 2017

%% -------------------------------------------
%           Prior Random Generator
%---------------------------------------------

% Choose prior type
switch type
    case 't'       % Student-t distribution
        % para1 = mean, para2 = s.d.
        nu = 2.1;  % degrees of freedom
        x = mvt_rnd(para1,para2^2*(nu-2)/nu,nu,1);
    case 'G'       % Gamma distribution
        % para1 = mean, para2 = s.d.
        a = para1^2/para2^2;      % shape
        b = para2^2/para1;        % scale
        x = gamrnd(a,b);
    case 'N'       % Normal distribution
        % para1 = mean, para2 = s.d.
        x = normrnd(para1,para2);
    case 'B'       % Beta distribution
        % para1 = mean, para2 = s.d.
        a = -para1*(para2^2+para1^2-para1)/para2^2;   % shape
        b = (para1-1)*(para2^2+para1^2-para1)/para2^2;% shape
        x = betarnd(a,b);
    case 'I1'      % Inverse-Gamma type-1 distribution
        % para1 = nu, para2 = s; see An & Schorfheide (2007)
        % Comment by Bing Li: s.d. being wildly large relative to mean would
        % affect numerical accuracy of simulated data. Simulated mean may be
        % ok, but simulated s.d. can be sharply smaller than the theoretical value.
        a = para1/2;
        b = a*para2^2;
        x = sqrt(1/gamrnd(a,1/b));
    case 'I2'      % Inverse-Gamma type-2 distribution
        % para1 = mean, para2 = s.d.
        a = 2+para1^2/para2^2;    % shape
        b = para1*(a-1);          % rate
        x = 1/gamrnd(a,1/b);
    case 'U'       % Uniform distribution
        % para1 = lb, para2 = ub
        x = unifrnd(para1,para2);
    case 'C'       % Categorical distribution
        % para1 = index vector, para2 = prob vector
        x = randsample(para1,1,true,para2);
    otherwise
        error('Prior type does not exist!')
end

%-------------------- END --------------------
