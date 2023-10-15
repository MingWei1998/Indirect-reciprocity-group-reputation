%% Reputation Evolution in Higher-Order Games
% advanced code
% 2022.5.18
function[Lcr,Lcrspt] = rep_evol_sym_cr(strategy_list,strNum,Lnum,Cnum,Dnum,N,ngs,repItnum,kmin,kmax,Gamestc,lambda_i,lambda_e,q,err)
%% strategy distribution
%strNum = input('Type of Leading Strategy:'); % as written
%strNum = 1;
str = strategy_list(strNum,:);
%tic % start timing
Lpop = 1:Lnum; 
Cpop = Lnum+1:Lnum+Cnum;
Dpop = Lnum+Cnum+1:N;  % distribute strategies
%strRt = [Lnum;Cnum;Dnum]/N; % strategy ratio init
%mutnum = 0; % number of mutations 

%% reputation init
% parameter setup
%k = 6; % size of group
%M = [randi([0,1],Lnum,N);ones(Cnum,N);zeros(Dnum,N)]; % initial reputation matrix(random)
M = [ones(Lnum,N);ones(Cnum,N);zeros(Dnum,N)]; % initial reputation matrix(good)
%M = [zeros(Lnum,N);ones(Cnum,N);zeros(Dnum,N)]; % initial reputation matrix(bad)

%lambda_e = 0.1; % external parameter
%lambda_i = 0.1; % internal parameter

%% cooperation details
% GameNum = 0;
% CoopNum = 0;
% CoopRat = zeros(1,Itnum); %

%% payoff & cooperation calculation
% ngs = size(Gamestc30,1); 
Gs = 1:ngs; 
% pof = zeros(1,N); 
% pofspt = zeros(ngs,N); 
GameTime = zeros(1,N); 
GameTimespt = zeros(ngs,N); 
CoopTime = zeros(1,N); 
CoopTimespt = zeros(ngs,N); 

%% game process: reputation&strategy evolution
for it  = 1:repItnum    
    k = randi([kmin,kmax]); %random size of group
    %k = 2;    
    % agent selection
    temp = randperm(N);
    agtlist = sort(temp(1:k)); 
    Lagt = intersect(agtlist,Lpop); 
    Cagt = intersect(agtlist,Cpop); 
    Dagt = intersect(agtlist,Dpop); 
    
    gamestc = Gs(Gamestc(:,1) == length(Lagt) & Gamestc(:,2) == length(Cagt) & Gamestc(:,3) == length(Dagt));
    restagt = setdiff(Lpop,agtlist); 
    
    % acting stage & inner reputation update
    
    [acts,inre] = inner_game2(agtlist,M,lambda_i,Lpop,Cpop,Dpop,str,k);
    M(Lagt,agtlist) = inre; 
        
    % outer reputation update    
    
    exre = outer_eval(agtlist,M,lambda_e,acts,Lpop,str,k,q,err);
    M(restagt,agtlist) = exre; 

    % payoff & cooperation update
    
    if it >= repItnum/2 
%         R = alpha*(k^beta); %synergy factor
%         Defpof = R*c*sum(acts)/k;
%         Cooppof = Defpof - c;    
%         agtpof = acts*Cooppof + (~acts)*Defpof; 
%         pof(agtlist) = pof(agtlist)+agtpof; 
        GameTime(agtlist) = GameTime(agtlist)+1; 
        CoopTime(agtlist) = CoopTime(agtlist)+acts; 
        
%         pofspt(gamestc,agtlist) = pofspt(gamestc,agtlist)+agtpof; 
        GameTimespt(gamestc,agtlist) = GameTimespt(gamestc,agtlist)+1; 
        CoopTimespt(gamestc,agtlist) = CoopTimespt(gamestc,agtlist)+acts; 
    end   
end

%% cooperation rate
% meanpof = pof./GameTime; 
% strmeanpof = [mean(meanpof(Lpop)),mean(meanpof(Cpop)),mean(meanpof(Dpop))]; 
% meanpofspt = pofspt./GameTimespt; 
% meanpofspt(isnan(meanpofspt)) = 0; 
% strmeanpofspt = [mean(meanpofspt(:,Lpop),2),mean(meanpofspt(:,Cpop),2),mean(meanpofspt(:,Dpop),2)]; 

popcooprat = CoopTime(Lpop)./GameTime(Lpop); 
Lcr = mean(popcooprat); 
popcoopratspt = CoopTimespt(:,Lpop)./GameTimespt(:,Lpop); 
popcoopratspt(isnan(popcoopratspt)) = 0; 
Lcrspt = mean(popcoopratspt,2)'; 
end
