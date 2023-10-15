%% coop_rate2payoff 
function[SMP] = cr2pofv(Pop,Game,LCRspt,alpha,beta,c,pk,Lstcprob,Cstcprob,Dstcprob)
nPopstc = size(Pop,1);
nGamestc = size(Game,1);
SMP = zeros(nPopstc,3);
LCRspt(isnan(LCRspt)) = 0;
k = sum(Game,2)'; 
k1 = Game(:,1); k2 = Game(:,2);
for nPs = 1:nPopstc    
    cr = LCRspt(nPs,:)'; 
    pi = zeros(nGamestc,4);
    for nGs = 1:nGamestc
        lnc = 0:k1(nGs);
        bino1 = binopdf(lnc(1:k1(nGs)),k1(nGs)-1,cr(nGs));
        bino2 = binopdf(lnc,k1(nGs),cr(nGs));
        payoff1 = payoff(k(nGs),lnc+k2(nGs)+1,alpha,beta,c,1);
        payoff2 = payoff(k(nGs),lnc+k2(nGs),alpha,beta,c,1);
        pi(nGs,1) = bino1*payoff1(1:k1(nGs))';
        pi(nGs,2) = bino1*(payoff2(1:k1(nGs))+c)';
        pi(nGs,3) = bino2*payoff2';
        pi(nGs,4) = bino2*(payoff2+c)';
    end
    SMP(nPs,1) = (pk(k-1).*Lstcprob(nPs,:))*(cr.*pi(:,1)+(1-cr).*pi(:,2));
    SMP(nPs,2) = (pk(k-1).*Cstcprob(nPs,:))*pi(:,3);
    SMP(nPs,3) = (pk(k-1).*Dstcprob(nPs,:))*pi(:,4);
end
end