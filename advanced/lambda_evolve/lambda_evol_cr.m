% co-evolution
% 
function[CR,crspt1,crspt2,crspt3] = lambda_evol_cr(str,popstrc,N,ngs,repItnum,kmin,kmax,Gamestc,q,err,lambda_list)
%% init
%la1 = lambda_list(1); la2 = lambda_list(2); la3 = lambda_list(3);% 
num1 = popstrc(1); num2 = popstrc(2); num3 = popstrc(3); % 
pop = 1:N;pop1 = 1:num1; pop2 = num1+1:num1+num2; pop3 = num1+num2+1:N; % 
%popla = [la1*ones(1,num1),la2*ones(1,num2),la3*ones(1,num3)]; % 

%% reputation init
%M = randi([0,1],N,N); % 
M = ones(N,N); % 
%M = zeros(N,N); % 

%% cooperation calculate init
Gs = 1:ngs; %
GameTime = zeros(1,N); %
GameTimespt = zeros(ngs,N); %
CoopTime = zeros(1,N); %
CoopTimespt = zeros(ngs,N); %

%% game process: reputation & lambda evolution
for it  = 1:repItnum  
    k = randi([kmin,kmax]); %random size of group
    temp = randperm(N);
    agtlist = sort(temp(1:k)); %
    %
    agt1 = intersect(agtlist,pop1); 
    agt2 = intersect(agtlist,pop2); 
    agt3 = intersect(agtlist,pop3); 
    K = [length(agt1),length(agt2),length(agt3)];rK = [num1,num2,num3]-K; 
    
    %
    gamestc = Gs(Gamestc(:,1) == length(agt1) & Gamestc(:,2) == length(agt2) & Gamestc(:,3) == length(agt3));
    restagt = setdiff(pop,agtlist); %

    % acting stage & inner reputation update
    % 
    [acts,inre] = lainner_game(agtlist,M,lambda_list,K,str,k);
    M(agtlist,agtlist) = inre; %

    % outer reputation update    
    % 
    exre = laouter_eval(agtlist,M,lambda_list,acts,pop,rK,str,k,q,err);
    M(restagt,agtlist) = exre; %

    if it >= repItnum/2 %
        GameTime(agtlist) = GameTime(agtlist)+1; %
        CoopTime(agtlist) = CoopTime(agtlist)+acts; %
        GameTimespt(gamestc,agtlist) = GameTimespt(gamestc,agtlist)+1; %
        CoopTimespt(gamestc,agtlist) = CoopTimespt(gamestc,agtlist)+acts; %
    end
end

%% 
popcooprat = CoopTime(pop)./GameTime(pop); %
cr1 = mean(popcooprat(pop1)); cr2 = mean(popcooprat(pop2));cr3 = mean(popcooprat(pop3));%
CR = [cr1,cr2,cr3]; 
popcoopratspt = CoopTimespt./GameTimespt; %
popcoopratspt(isnan(popcoopratspt)) = 0; %
crspt1 = mean(popcoopratspt(:,pop1),2)';crspt2 = mean(popcoopratspt(:,pop2),2)'; crspt3 = mean(popcoopratspt(:,pop3),2)'; %
%CRspt = [crspt1,crspt2,crspt3];


end