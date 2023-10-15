%% calculate fixation probability
tic
%Default;
pk = [kmin:kmax]/sum(kmin:kmax);
Gametype = [1,2;1,3;2,3];
for strNum = 1:8
    FXP = zeros(9,6);
    for nGs = 1:9
        lambda = 0.1*nGs;
        %Dataname = strcat('Indirect-reciprocity-group-reputation\results\formal\strategy evolve\cr2pof\coop_rate\L',num2str(strNum),'\λ=',num2str(lambda),'_L',num2str(strNum),'CR.mat');
        %load(Dataname);        
        %SMP = cr2pofv(Popstc30,Gamestc30,LCRspt,alpha,beta,c,pk,Lstcprob,Cstcprob,Dstcprob);
        SMPname = strcat('Indirect-reciprocity-group-reputation\results\formal\strategy evolve\cr2pof\μ→0&SMP\β=',num2str(beta),',s=1\α=',num2str(alpha),',β=',num2str(beta),',s=1\SMP\L',num2str(strNum),'_λ=',num2str(lambda),'_SMP.mat');
        %save(SMPname,'SMP');
        load(SMPname);
        for gametype = 1:3
            str1 = Gametype(gametype,1);
            str2 = Gametype(gametype,2);
            fxp = fixprob(s,str1,str2,Popstc30,SMP);
            FXP(nGs,2*gametype-1:2*gametype) = fxp;
        end
    end
    %filename = strcat('Indirect-reciprocity-group-reputation\results\formal\strategy evolve\cr2pof\μ→0&SMP\β=',num2str(beta),',s=',num2str(s),'\α=',num2str(alpha),',β=',num2str(beta),',s=',num2str(s),'\fix_prob\L',num2str(strNum),'fixprob.mat');
    filename = strcat('Indirect-reciprocity-group-reputation\results\formal\strategy evolve\cr2pof\s_changing\s=',num2str(s),'\α=',num2str(alpha),'\fix_prob\L',num2str(strNum),'fixprob.mat');
    save(filename,'FXP');
    %disp(strNum);
end
%
%Default;
for strNum = 1:8
    allSMEq = zeros(9,3);
    allcoop = zeros(9,1);
    %Dataname1 = strcat('Indirect-reciprocity-group-reputation\results\formal\strategy evolve\formal\fixation probability\L',num2str(strNum),'fxp.mat');
    %Dataname1 = strcat('Indirect-reciprocity-group-reputation\results\formal\strategy evolve\cr2pof\μ→0\α=',num2str(alpha),',β=1,s=1\fix_prob\L',num2str(strNum),'fixprob.mat');
    %Dataname1 = strcat('Indirect-reciprocity-group-reputation\results\formal\strategy evolve\cr2pof\μ→0&SMP\β=',num2str(beta),',s=',num2str(s),'\α=',num2str(alpha),',β=',num2str(beta),',s=',num2str(s),'\fix_prob\L',num2str(strNum),'fixprob.mat');
    load(filename);    
    for nGs = 1:9
        tFXP = zeros(3,3);
        lambda = 0.1*nGs;
        %Dataname2 = strcat('Indirect-reciprocity-group-reputation\results\formal\strategy evolve\formal\L',num2str(strNum),'\λ=',num2str(lambda),'.mat');
        Dataname2 = strcat('Indirect-reciprocity-group-reputation\results\formal\strategy evolve\cr2pof\coop_rate\L',num2str(strNum),'\λ=',num2str(lambda),'_L',num2str(strNum),'CR.mat');
        load(Dataname2);
        lcr = LCR(1);%
        tFXP(2,1) = FXP(nGs,1);%C trans2 Ld
        tFXP(1,2) = FXP(nGs,2);%Ld trans2 C
        tFXP(3,1) = FXP(nGs,3);%D trans2 Ld
        tFXP(1,3) = FXP(nGs,4);%Ld trans2 D
        tFXP(3,2) = FXP(nGs,5);%D trans2 C
        tFXP(2,3) = FXP(nGs,6);%C trans2 D
        T = tFXP/2;
        for gametype = 1:3
            T(gametype,gametype) = 1-sum(T(gametype,:));
        end
        v=null(T'-eye(3)); 
        SMEq=v'/sum(v);
        allSMEq(nGs,:) = SMEq;
        coop = SMEq*[lcr;1;0];
        allcoop(nGs) = coop;
    end
    filename = strcat('Indirect-reciprocity-group-reputation\results\formal\strategy evolve\cr2pof\μ→0&SMP\β=',num2str(beta),',s=',num2str(s),'\α=',num2str(alpha),',β=',num2str(beta),',s=',num2str(s),'\SMEq&popcr\L',num2str(strNum),'_SMEq&coop.mat');
    save(filename,'allSMEq','allcoop');
end
toc