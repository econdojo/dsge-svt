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
E_y(V.INF) = 0;
E_y(V.INT) = 0;

% 2nd-order raw moments
E_yy(V.YGR,V.YGR,1) = 0;
E_yy(V.YGR,V.INF,1) = 0;
E_yy(V.YGR,V.INT,1) = 0;
E_yy(V.INF,V.INF,1) = 0;
E_yy(V.INF,V.INT,1) = 0;
E_yy(V.INT,V.INT,1) = 0;

E_yy(V.YGR,V.YGR,2) = 0;
E_yy(V.YGR,V.INF,2) = 0;
E_yy(V.YGR,V.INT,2) = 0;
E_yy(V.INF,V.YGR,2) = 0;
E_yy(V.INF,V.INF,2) = 0;
E_yy(V.INF,V.INT,2) = 0;
E_yy(V.INT,V.YGR,2) = 0;
E_yy(V.INT,V.INF,2) = 0;
E_yy(V.INT,V.INT,2) = 0;

% E_yy(V.YGR,V.YGR,3) = 0;
% E_yy(V.YGR,V.INF,3) = 0;
% E_yy(V.YGR,V.INT,3) = 0;
% E_yy(V.INF,V.YGR,3) = 0;
% E_yy(V.INF,V.INF,3) = 0;
% E_yy(V.INF,V.INT,3) = 0;
% E_yy(V.INT,V.YGR,3) = 0;
% E_yy(V.INT,V.INF,3) = 0;
% E_yy(V.INT,V.INT,3) = 0;

% E_yy(V.YGR,V.YGR,4) = 0;
% E_yy(V.YGR,V.INF,4) = 0;
% E_yy(V.YGR,V.INT,4) = 0;
% E_yy(V.INF,V.YGR,4) = 0;
% E_yy(V.INF,V.INF,4) = 0;
% E_yy(V.INF,V.INT,4) = 0;
% E_yy(V.INT,V.YGR,4) = 0;
% E_yy(V.INT,V.INF,4) = 0;
% E_yy(V.INT,V.INT,4) = 0;

% E_yy(V.YGR,V.YGR,5) = 0;
% E_yy(V.YGR,V.INF,5) = 0;
% E_yy(V.YGR,V.INT,5) = 0;
% E_yy(V.INF,V.YGR,5) = 0;
% E_yy(V.INF,V.INF,5) = 0;
% E_yy(V.INF,V.INT,5) = 0;
% E_yy(V.INT,V.YGR,5) = 0;
% E_yy(V.INT,V.INF,5) = 0;
% E_yy(V.INT,V.INT,5) = 0;

%% User Input End Here