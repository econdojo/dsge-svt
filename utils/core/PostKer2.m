function postker = PostKer2(theta,block,P,data,mixid)
% Function POSTKER2
%
% Purpose:    Evaluate log posterior kernel density
%
% Format:     postker = PostKer2(theta,block,P,data,mixid)
%
% Input:      theta     sv parameters
%             block     sv parameter indices
%             P         structure of sv parameters
%             data      transformed data
%             mixid     mixture normal indices
%
% Output:     postker   (-)log posterior kernel density
%
% Written by Fei Tan, Saint Louis University
% Updated: January 6, 2020

%% -------------------------------------------
%          Evaluate Posterior Kernel
%---------------------------------------------

if any(theta'<P.lb(block)) || any(theta'>P.ub(block))
    postker = 1e10;
    return
else
    % Prior
    logprior = 0;
    for k = 1:3
        logprior = logprior+prior_pdf(theta(k),P.para1(block(k)),P.para2(block(k)),P.type{block(k)});
    end
    
    % Posterior (integrating out log vol)
    mu = theta(1);
    phi = theta(2);
    sig2 = theta(3);
    [p10,m10,v10] = mix10;
    for k = 1:10
        SSR.C{k} = (1-phi)*mu;
        SSR.G{k} = phi;
        SSR.Sigma_e{k} = sig2;
        SSR.D{k} = m10(k);
        SSR.Z{k} = 1;
        SSR.Sigma_u{k} = v10(k);
    end
    [~,~,loglik,~,~] = GibbsSmoother_mex(data,SSR,repmat(p10,10,1),mixid,false,false);
    postker = -(sum(loglik)+logprior);
end

%-------------------- END --------------------
