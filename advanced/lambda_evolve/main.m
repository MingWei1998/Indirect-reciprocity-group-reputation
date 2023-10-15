% 
Default;
%strNum = 1;
N = max(Popstc30(:,1));
str = strategy_list(strNum,:); %
nPopstc = size(Popstc30,1); %
ngs = size(Gamestc30,1); %
lambda_list = [0.1,0.4,0.7];
%init
allCR = zeros(nPopstc,3);
CRspt1 = zeros(nPopstc,ngs);CRspt2 = zeros(nPopstc,ngs);CRspt3 = zeros(nPopstc,ngs);
for nPs = 1:nPopstc
    popstrc = Popstc30(nPs,:);
    [CR,crspt1,crspt2,crspt3] = lambda_evol_cr(str,popstrc,N,ngs,repItnum,kmin,kmax,Gamestc30,q,err,lambda_list);
    allCR(nPs,:) = CR;
    CRspt1(nPs,:) = crspt1; CRspt2(nPs,:) = crspt2; CRspt3(nPs,:) = crspt3;
end
filename = strcat('Indirect-reciprocity-group-reputation\results\formal\lambda evolve\coop_rate\147\initgood\L',num2str(strNum),'cooprate');
save(filename,'allCR','CRspt1','CRspt2','CRspt3');


