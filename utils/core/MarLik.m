function [logmlik,nse] = MarLik(method,savepath)
% Function MARLIK
%
% Purpose:    Compute log marginal likelihood
%
% Format:     [logmlik,nse] = MarLik(method,savepath)
%
% Input:      method    Chib identity/harmonic mean
%             savepath  save path
%
% Output:     logmlik   (log) marginal likelihood
%             nse       numerical standard error
%
% Written by Fei Tan, Saint Louis University
% Updated: July 18, 2017

%% -------------------------------------------
%             Marginal Likelihood
%---------------------------------------------

% Initialization
step = [0.02 0.1:0.1:1];               % percentage step size
load([savepath filesep 'tarb_prior.mat'],'c')    % normalization constant for truncated prior

% Choose method
switch method
    case 'chib'                        % see Chib (1995), Chib & Jeliazkov (2001)
        % Initialization
        load([savepath filesep 'tarb_reduce.mat'],'num','den','pk')  % load summands & pk at point of evaluation
        [M,B] = size(num);             % sample size & number of blocks
        s = zeros(11,B);               % recursive post ordinate
        logmlik = zeros(11,1);         % recursive log marginal likelihood
        
        % Recursive log marginal likelihood
        for k = 1:11
            for j = 1:B
                s(k,j) = log(mean(num(1:round(M*step(k)),j))/mean(den(1:round(M*step(k)),j)));
            end
            logmlik(k) = pk-sum(s(k,:))-log(c);
        end
        
        % Numerical standard error
        H = [num den];
        H = H-repmat(mean(H),M,1);
        dev = [1./mean(num) -1./mean(den)];
        nse = sqrt(dev*NeweyWest(H)*dev'/M);
        
        % Determine plot dimension
        if B<=2
            dim1 = 1;
        elseif B>=3 && B<=4
            dim1 = 2;
        elseif B>=5 && B<=6
            dim1 = 3;
        elseif B>=7
            dim1 = 4;
        end
        dim2 = ceil(B/dim1);
        
        % Plot recursive post ordinate
        figure('Name','Reduced Posterior Ordinate')
        for j = 1:B
            subplot(dim1,dim2,j);
            plot(100*step,s(:,j),'r','linewidth',1.5);
            hold on; grid on; box on;
            title(['Stage ',num2str(j)],'fontsize',12,'fontname','times');
            xlabel('% of sample')
        end
        saveas(gcf,[savepath filesep 'Fig_PostOrd.fig']);
        save([savepath filesep 'tarb_reduce.mat'],'logmlik','nse','-append');
    case 'harm'                        % see Gelfand & Dey (1994), Geweke (1999)
        % Initialization
        load([savepath filesep 'tarb_full.mat'],'chain_para','pk')   % load MCMC draws & post kernels
        [M,np] = size(chain_para);     % sample size & number of parameters
        smean = mean(chain_para);      % post mean vector
        svar = cov(chain_para);        % post covariance matrix
        pk_mean = mean(pk);            % post kernel mean
        R = cholmod(svar); svar = R'*R;
        tau = 0.1:0.2:0.9;             % truncation parameters
        bound = chi2inv(tau,np);       % bound for indicator function
        s = zeros(M,5);                % initialize summand
        logmlik = zeros(11,5);         % recursive log marginal likelihood
        nse = zeros(5,1);              % numerical standard error
        
        % Harmonic mean estimator
        for k = 1:M
            ind = find(((chain_para(k,:)-smean)/svar)*((chain_para(k,:)-smean)')<=bound);
            if ~isempty(ind)
                s(k,ind) = 1./tau(ind)*exp(mvt_pdf(chain_para(k,:),smean,svar,inf))/exp(pk(k)-pk_mean);
            end
        end
        
        % Recursive log marginal likelihood
        for k = 1:11
            logmlik(k,:) = -log(mean(s(1:round(M*step(k)),:)))+pk_mean-log(c);
        end
        
        % Numerical standard error
        for k = 1:5
            nse(k) = sqrt(NeweyWest(s(:,k)-mean(s(:,k)))/M)/mean(s(:,k));
        end
        save([savepath filesep 'tarb_full.mat'],'logmlik','nse','-append');
    otherwise
        error('Method does not exist.')
end

% Plot recursive marginal likelihood
figure('Name','Marginal Likelihood')
plot(100*step,logmlik,'linewidth',1.5);
hold on;   grid on;   box on;
title('Marginal Likelihood','fontsize',12,'fontname','times');
if strcmp(method,'harmean')
    legend('p = 0.1','p = 0.3','p = 0.5','p = 0.7','p = 0.9');
end
xlabel('Percent of sample')
saveas(gcf,[savepath filesep 'Fig_MarLik.fig']);

%-------------------- END --------------------
