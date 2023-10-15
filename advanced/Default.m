% blendent file
%blendent = [L1_color,L2_color,L3_color,L4_color,L5_color,L6_color,L7_color,L8_color];

%% well-mixed population
% population size
%n = 50;
%type_num = 3;
%Lratio = 1; Cratio = 1; Dratio = 1; 
%Lnum = n*Lratio; Cnum = n*Cratio; Dnum = n*Dratio; 
%N = Lnum+Cnum+Dnum; 
%pop = 1:N; 

%% parameter setup
repItnum = 4e5; 
%strItnum = 1e5; 
err = 0; % error probability
q = 1; % observing probability
kmax = 10; 
kmin = 2; 
nk = kmax-kmin+1; 

%% payoff parameter
%alpha = 1;
beta = 1;
%delta = 1;
c = 1; % cost
%mu = 0.1; % mutation rate
s = 1; 

%% strategy setup
% leading eight and ALLC/ALLD
L_1 = [0,0,1,1,1,0,1,1,1,1,0,1];
L_2 = [0,0,1,1,1,0,0,1,1,1,0,1];
L_3 = [1,0,1,1,1,0,1,1,0,1,0,1];
L_4 = [1,0,0,1,1,0,1,1,0,1,0,1];
L_5 = [1,0,1,1,1,0,0,1,0,1,0,1];
L_6 = [1,0,0,1,1,0,0,1,0,1,0,1];
L_7 = [0,0,0,1,1,0,1,1,0,1,0,1];
L_8 = [0,0,0,1,1,0,0,1,0,1,0,1];
ALLC = [1,1,1,1,1,1,1,1,1,1,1,1];
ALLD = [0,0,0,0,0,0,0,0,0,0,0,0];
strategy_list = [L_1; L_2; L_3; L_4; L_5; L_6; L_7; L_8; ALLC; ALLD];
%situation_judge_alpha = [0,0,0;0,0,1;0,1,0;0,1,1;1,0,0;1,0,1;1,1,0;1,1,1];
%situation_judge_beta = [0,0;0,1;1,0;1,1];
