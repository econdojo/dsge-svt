% Function USER_BMM
%
% Purpose:    User input - moment conditions
%
% Format:     user_bmm
%
% Input:      P         parameter structure
%             V         variable structure
%
% Output:     E_y       1st-order moments: E[y(t)] = E_y
%             E_yy      2nd-order raw moments: E[y(t)y(t-k)'] = E_yy(:,:,k+1)
%
% Written by Fei Tan, Saint Louis University
% Updated: June 10, 2022

%% User Input Start Here

% 1st-order moments
E_y(V.YGR) = 0;
E_y(V.CGR) = 0;
E_y(V.IGR) = 0;
E_y(V.WGR) = 0;
E_y(V.HOUR) = 0;
E_y(V.INF) = 0;
E_y(V.FFR) = 0;

% 2nd-order raw moments
E_yy(V.YGR,V.YGR,1) = 0;
E_yy(V.YGR,V.CGR,1) = 0;
E_yy(V.YGR,V.IGR,1) = 0;
E_yy(V.YGR,V.WGR,1) = 0;
E_yy(V.YGR,V.HOUR,1) = 0;
E_yy(V.YGR,V.INF,1) = 0;
E_yy(V.YGR,V.FFR,1) = 0;
E_yy(V.CGR,V.CGR,1) = 0;
E_yy(V.CGR,V.IGR,1) = 0;
E_yy(V.CGR,V.WGR,1) = 0;
E_yy(V.CGR,V.HOUR,1) = 0;
E_yy(V.CGR,V.INF,1) = 0;
E_yy(V.CGR,V.FFR,1) = 0;
E_yy(V.IGR,V.IGR,1) = 0;
E_yy(V.IGR,V.WGR,1) = 0;
E_yy(V.IGR,V.HOUR,1) = 0;
E_yy(V.IGR,V.INF,1) = 0;
E_yy(V.IGR,V.FFR,1) = 0;
E_yy(V.WGR,V.WGR,1) = 0;
E_yy(V.WGR,V.HOUR,1) = 0;
E_yy(V.WGR,V.INF,1) = 0;
E_yy(V.WGR,V.FFR,1) = 0;
E_yy(V.HOUR,V.HOUR,1) = 0;
E_yy(V.HOUR,V.INF,1) = 0;
E_yy(V.HOUR,V.FFR,1) = 0;
E_yy(V.INF,V.INF,1) = 0;
E_yy(V.INF,V.FFR,1) = 0;
E_yy(V.FFR,V.FFR,1) = 0;

% E_yy(V.YGR,V.YGR,2) = 0;
% E_yy(V.YGR,V.INF,2) = 0;
% E_yy(V.YGR,V.INT,2) = 0;
% E_yy(V.INF,V.YGR,2) = 0;
% E_yy(V.INF,V.INF,2) = 0;
% E_yy(V.INF,V.INT,2) = 0;
% E_yy(V.INT,V.YGR,2) = 0;
% E_yy(V.INT,V.INF,2) = 0;
% E_yy(V.INT,V.INT,2) = 0;

%% User Input End Here