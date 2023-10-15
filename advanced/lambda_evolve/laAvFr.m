% 
erthd = 1e-14;
Default;
pk = [kmin:kmax]/sum(kmin:kmax);
%% 
Gametype = [1,2;1,3;2,3];
FXP = zeros(8,6);
for strNum = 1:8    
    Dataname = strcat('Indirect-reciprocity-group-reputation\results\formal\lambda evolve\coop_rate\147\initrand\L',num2str(strNum),'cooprate.mat');
    load(Dataname); 
    SMP = lacr2pofv(Popstc30,Gamestc30,CRspt1,CRspt2,CRspt3,alpha,beta,c,pk,Lstcprob,Cstcprob,Dstcprob);
    SMPname = strcat('Indirect-reciprocity-group-reputation\results\formal\lambda evolve\μ→0&SMP\β=',num2str(beta),'\α=',num2str(alpha),'\SMP\L',num2str(strNum),'_SMP.mat');
    save(SMPname,'SMP');
    for gametype = 1:3
        str1 = Gametype(gametype,1);
        str2 = Gametype(gametype,2);
        fxp = fixprob(s,str1,str2,Popstc30,SMP);
        FXP(strNum,2*gametype-1:2*gametype) = fxp;
    end
end
filename = strcat('Indirect-reciprocity-group-reputation\results\formal\lambda evolve\μ→0&SMP\β=',num2str(beta),'\α=',num2str(alpha),'\fixprob.mat');
save(filename,'FXP');

%% 
allSMEq = zeros(8,3);
allcoop = zeros(8,1);
Dataname1 = strcat(filename);
load(Dataname1);
for strNum = 1:8
    tFXP = zeros(3,3);
    Dataname = strcat('Indirect-reciprocity-group-reputation\results\formal\lambda evolve\coop_rate\147\initrand\L',num2str(strNum),'cooprate.mat');
    load(Dataname);
    lcr = [allCR(1,1),allCR(31,2),allCR(496,3)]; %
    tFXP(2,1) = FXP(strNum,1);%la2 trans2 la1
    tFXP(1,2) = FXP(strNum,2);%la1 trans2 la2
    tFXP(3,1) = FXP(strNum,3);%la3 trans2 la1
    tFXP(1,3) = FXP(strNum,4);%la1 trans2 la3
    tFXP(3,2) = FXP(strNum,5);%la3 trans2 la2
    tFXP(2,3) = FXP(strNum,6);%la2 trans2 la3
    T = tFXP/2;
    for gametype = 1:3
        T(gametype,gametype) = 1-sum(T(gametype,:));
    end
    [x,y] = eig(T');
    errdet = abs(diag(y)-1);
    if min(errdet)<erthd
        v = real(x(:,errdet==min(errdet)));
        SMEq=v'/sum(v);
    else 
        SMEq = [0,0,0];
        disp([strNum,lambda]);
    end
%         v=null(T'-eye(3));
%         SMEq = v'/sum(v);
    allSMEq(strNum,:) = SMEq;
    coop = SMEq*[lcr(1);lcr(2);lcr(3)];
    allcoop(strNum) = coop;
end
filename = strcat('Indirect-reciprocity-group-reputation\results\formal\lambda evolve\μ→0&SMP\β=',num2str(beta),'\α=',num2str(alpha),'\SMEq&coop.mat');
save(filename,'allSMEq','allcoop');




