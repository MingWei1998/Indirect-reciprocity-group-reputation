%% μ>>0
pk = [kmin:kmax]/sum(kmin:kmax);
for strNum = 1:8
    allcoop = zeros(9,1);
    allAvFr = zeros(9,3);
    for l = 1:9
        lambda = 0.1*l;
        Dataname1 = strcat('Indirect-reciprocity-group-reputation\results\formal\strategy evolve\cr2pof\μ→0&SMP\β=',num2str(beta),',s=',num2str(s),'\α=',num2str(alpha),',β=',num2str(beta),',s=',num2str(s),'\SMP\L',num2str(strNum),'_λ=',num2str(lambda),'_SMP.mat');
        load(Dataname1);
        Dataname2 = strcat('Indirect-reciprocity-group-reputation\results\formal\strategy evolve\cr2pof\coop_rate\L',num2str(strNum),'\λ=',num2str(lambda),'_L',num2str(strNum),'CR.mat');
        load(Dataname2);
        LCR(isnan(LCR)) = 0;
        [PopCoop,AvFr] = slc_mut_equ(mu,s,Popstc30,SMP,LCR);
        allcoop(l) = PopCoop;
        allAvFr(l,:) = AvFr';
    end
    filename = strcat('Indirect-reciprocity-group-reputation\results\formal\strategy evolve\cr2pof\μgg0\μ=',num2str(mu),'\β=',num2str(beta),',s=',num2str(s),'\α=',num2str(alpha),'\L',num2str(strNum),'coop&AvFr.mat');
    save(filename,'allAvFr','allcoop');
end