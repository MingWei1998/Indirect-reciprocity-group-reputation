%% selection-mutation equilibrium for lambda evolve
function[PopCoop,AvFr] = laSMEq(mu,s,Pop,SMP,LCR)
N = max(Pop(:,1)); %
nPopstc = size(Pop,1); %
W = zeros(nPopstc,nPopstc); %
for nPs = 1:nPopstc %
    num1 = Pop(nPs,1); num2 = Pop(nPs,2); num3 = Pop(nPs,3); %
    %
    wL2C = num1/N * (mu/2+(1-mu)*num2/(N-1)/(1+exp(-s*(SMP(nPs,2)-SMP(nPs,1)))));
    wL2D = num1/N * (mu/2+(1-mu)*num3/(N-1)/(1+exp(-s*(SMP(nPs,3)-SMP(nPs,1)))));
    wC2L = num2/N * (mu/2+(1-mu)*num1/(N-1)/(1+exp(-s*(SMP(nPs,1)-SMP(nPs,2)))));
    wC2D = num2/N * (mu/2+(1-mu)*num3/(N-1)/(1+exp(-s*(SMP(nPs,3)-SMP(nPs,2)))));
    wD2L = num3/N * (mu/2+(1-mu)*num1/(N-1)/(1+exp(-s*(SMP(nPs,1)-SMP(nPs,3)))));
    wD2C = num3/N * (mu/2+(1-mu)*num2/(N-1)/(1+exp(-s*(SMP(nPs,2)-SMP(nPs,3)))));
     %
    if num1>0 %
        W(nPs,Pop(:,1) == num1-1 & Pop(:,2) == num2+1) = wL2C; %
        W(nPs,Pop(:,1) == num1-1 & Pop(:,2) == num2) = wL2D; %     
    end
    if num2>0 %
        W(nPs,Pop(:,1) == num1+1 & Pop(:,2) == num2-1) = wC2L; %
        W(nPs,Pop(:,1) == num1 & Pop(:,2) == num2-1) = wC2D; %       
    end
    if num3>0 %
        W(nPs,Pop(:,1) == num1+1 & Pop(:,2) == num2) = wD2L; %
        W(nPs,Pop(:,1) == num1 & Pop(:,2) == num2+1) = wD2C; %      
    end
    W(nPs,nPs) = 1-(wL2C+wL2D+wC2L+wC2D+wD2L+wD2C); %
end
v = null(W'-eye(nPopstc)); %
if size(v,2) == 1 %
    SMEq = v/sum(v); %
    PopCoop = sum(Pop'.*LCR')/N*SMEq;
    %PopCoop = (Pop(:,1)'.*LCR(:,1)'+Pop(:,2)')/N*SMEq; %
    AvFr = Pop'*SMEq/N; %
else
    disp('error');
end
end