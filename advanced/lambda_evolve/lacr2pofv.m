% 
function[SMP] = lacr2pofv(Pop,Game,CRspt1,CRspt2,CRspt3,alpha,beta,c,pk,Lstcprob1,Lstcprob2,Lstcprob3)
nPopstc = size(Pop,1);
nGamestc = size(Game,1);
SMP = zeros(nPopstc,3);
%
CRspt1(isnan(CRspt1)) = 0;CRspt2(isnan(CRspt2)) = 0;CRspt3(isnan(CRspt3)) = 0;
k = sum(Game,2)'; %
k1 = Game(:,1); k2 = Game(:,2); k3 = Game(:,3);%
for nPs = 1:nPopstc
    cr = [CRspt1(nPs,:)',CRspt2(nPs,:)',CRspt3(nPs,:)']; %
    pi = zeros(nGamestc,6);
    for nGs = 1:nGamestc
        lnc1 = 0:k1(nGs);lnc2 = 0:k2(nGs);lnc3 = 0:k3(nGs); %
        binoex1 = binopdf(lnc1(1:k1(nGs)),k1(nGs)-1,cr(nGs,1));
        binoex2 = binopdf(lnc2(1:k2(nGs)),k2(nGs)-1,cr(nGs,2));
        binoex3 = binopdf(lnc3(1:k3(nGs)),k3(nGs)-1,cr(nGs,3));%
        binoin1 = binopdf(lnc1,k1(nGs),cr(nGs,1));
        binoin2 = binopdf(lnc2,k2(nGs),cr(nGs,2));
        binoin3 = binopdf(lnc3,k3(nGs),cr(nGs,3)); %
        payoffc = payoff(k(nGs),1:k(nGs),alpha,beta,c,1); %
        payoffd = [0,payoffc(1:end-1)+c]; %
        %
        %
        for j1 = 0:k1(nGs)-1
            for j2 = 0:k2(nGs)
                for j3 = 0:k3(nGs)
                    pi(nGs,1) = pi(nGs,1)+binoex1(j1+1)*binoin2(j2+1)*binoin3(j3+1)*payoffc(j1+j2+j3+1);
                    pi(nGs,2) = pi(nGs,2)+binoex1(j1+1)*binoin2(j2+1)*binoin3(j3+1)*payoffd(j1+j2+j3+1);
                end
            end
        end
        %
        for j1 = 0:k2(nGs)-1
            for j2 = 0:k1(nGs)
                for j3 = 0:k3(nGs)
                    pi(nGs,3) = pi(nGs,3)+binoex2(j1+1)*binoin1(j2+1)*binoin3(j3+1)*payoffc(j1+j2+j3+1);
                    pi(nGs,4) = pi(nGs,4)+binoex2(j1+1)*binoin1(j2+1)*binoin3(j3+1)*payoffd(j1+j2+j3+1);
                end
            end
        end
        %
        for j1 = 0:k3(nGs)-1
            for j2 = 0:k1(nGs)
                for j3 = 0:k2(nGs)
                    pi(nGs,5) = pi(nGs,5)+binoex3(j1+1)*binoin1(j2+1)*binoin2(j3+1)*payoffc(j1+j2+j3+1);
                    pi(nGs,6) = pi(nGs,6)+binoex3(j1+1)*binoin1(j2+1)*binoin2(j3+1)*payoffd(j1+j2+j3+1);
                end
            end
        end       
    end
    SMP(nPs,1) = (pk(k-1).*Lstcprob1(nPs,:))*(cr(:,1).*pi(:,1)+(1-cr(:,1)).*pi(:,2));
    SMP(nPs,2) = (pk(k-1).*Lstcprob2(nPs,:))*(cr(:,2).*pi(:,3)+(1-cr(:,2)).*pi(:,4));
    SMP(nPs,3) = (pk(k-1).*Lstcprob3(nPs,:))*(cr(:,3).*pi(:,5)+(1-cr(:,3)).*pi(:,6));
end
end