function data = DatMat(Y,mod,nlag)
% Function DATMAT
%
% Purpose:    Construct data matrices for VAR/BMM
%
% Format:     data = DatMat(Y,mod,nlag)
%
% Input:      Y         data set with ROW observations
%             mod       model type
%             nlag      number of lags
%
% Output:     data      data matrices (structure)
%
% Written by Fei Tan, Saint Louis University
% Updated: June 10, 2022

%% -------------------------------------------
%           Construct Data Matrices
%---------------------------------------------

T = size(Y,1)-nlag;            % time span

switch mod
    case 'var'
        data.X = ones(T,1);            % regressors
        for k = 1:nlag
            data.X = cat(2,data.X,Y((nlag-k+1):(nlag-k+T),:));
        end
        data.Y = Y((nlag+1):end,:);    % exclude initial p lags
    case 'bmm'
        [P,V] = ParVar;                %#ok<ASGLU> parameter & variable indices
        n = size(Y,2);                 % number of observables
        E_y = zeros(n,1)-1;            % 1st-order moments
        E_yy = zeros(n,n,nlag+1)-1;    % 2nd-order raw moments
        user_bmm
        M = [];
        for k = 1:n
            if E_y(k)>=0
                M = cat(2,M,Y(nlag+1:end,k));
            end
        end
        for k = 1:nlag+1
            for i = 1:n
                for j = 1:n
                    if E_yy(i,j,k)>=0
                        M = cat(2,M,Y(nlag+1:end,i).*Y(nlag+2-k:end+1-k,j));
                    end
                end
            end
        end
        data.M = mean(M);
        data.V = NeweyWest(M-repmat(data.M,T,1));
        data.E_y = E_y;
        data.E_yy = E_yy;
end

%-------------------- END --------------------