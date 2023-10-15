%% selection-mutation equilibrium for 3 strategies
function[PopCoop,AvFr] = slc_mut_equ(mu,s,Pop,SMP,LCR)
N = max(Pop(:,1)); 
nPopstc = size(Pop,1); 
W = zeros(nPopstc,nPopstc); 
for nPs = 1:nPopstc 
    Lnum = Pop(nPs,1); Cnum = Pop(nPs,2); Dnum = Pop(nPs,3); 
    
    wL2C = Lnum/N * (mu/2+(1-mu)*Cnum/(N-1)/(1+exp(-s*(SMP(nPs,2)-SMP(nPs,1)))));
    wL2D = Lnum/N * (mu/2+(1-mu)*Dnum/(N-1)/(1+exp(-s*(SMP(nPs,3)-SMP(nPs,1)))));
    wC2L = Cnum/N * (mu/2+(1-mu)*Lnum/(N-1)/(1+exp(-s*(SMP(nPs,1)-SMP(nPs,2)))));
    wC2D = Cnum/N * (mu/2+(1-mu)*Dnum/(N-1)/(1+exp(-s*(SMP(nPs,3)-SMP(nPs,2)))));
    wD2L = Dnum/N * (mu/2+(1-mu)*Lnum/(N-1)/(1+exp(-s*(SMP(nPs,1)-SMP(nPs,3)))));
    wD2C = Dnum/N * (mu/2+(1-mu)*Cnum/(N-1)/(1+exp(-s*(SMP(nPs,2)-SMP(nPs,3)))));
    
%     wL2C = Lnum/N * (mu/2+(1-mu)*Cnum/(N-1));
%     wL2D = Lnum/N * (mu/2+(1-mu)*Dnum/(N-1));
%     wC2L = Cnum/N * (mu/2+(1-mu)*Lnum/(N-1));
%     wC2D = Cnum/N * (mu/2+(1-mu)*Dnum/(N-1));
%     wD2L = Dnum/N * (mu/2+(1-mu)*Lnum/(N-1));
%     wD2C = Dnum/N * (mu/2+(1-mu)*Cnum/(N-1));

    
    if Lnum>0 
        W(nPs,Pop(:,1) == Lnum-1 & Pop(:,2) == Cnum+1) = wL2C; 
        W(nPs,Pop(:,1) == Lnum-1 & Pop(:,2) == Cnum) = wL2D;        
    end
    if Cnum>0 
        W(nPs,Pop(:,1) == Lnum+1 & Pop(:,2) == Cnum-1) = wC2L; 
        W(nPs,Pop(:,1) == Lnum & Pop(:,2) == Cnum-1) = wC2D;        
    end
    if Dnum>0 
        W(nPs,Pop(:,1) == Lnum+1 & Pop(:,2) == Cnum) = wD2L; 
        W(nPs,Pop(:,1) == Lnum & Pop(:,2) == Cnum+1) = wD2C;       
    end
    W(nPs,nPs) = 1-(wL2C+wL2D+wC2L+wC2D+wD2L+wD2C); 
end
v = null(W'-eye(nPopstc)); 
if size(v,2) == 1 
    %nv = v./sum(v); 
    SMEq = v/sum(v); 
    %LiCoop = LCR(:,1)'*SMEq; 
    PopCoop = (Pop(:,1)'.*LCR(:,1)'+Pop(:,2)')/N*SMEq; 
    AvFr = Pop'*SMEq/N; 
else
    disp('error');
end
end