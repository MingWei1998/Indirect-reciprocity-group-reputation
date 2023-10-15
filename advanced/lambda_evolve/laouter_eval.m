% 
function[exre] = laouter_eval(agtlist,RM,lambda_list,acts,pop,rK,str,k,q,err)
restagt = setdiff(pop,agtlist); %
rk = sum(rK); %
exgre1 = zeros(rK(1),k); exgre2 = zeros(rK(2),k);exgre3 = zeros(rK(3),k);exre = zeros(rk,k); %init
%ragtla = [lambda_list(1)*ones(1,rK(1)),lambda_list(2)*ones(1,rK(2)),lambda_list(3)*ones(1,rK(3))];

for j1 = 1:k %
    foc = agtlist(j1); %
    eval = agtlist(agtlist~=foc); %
    resum = sum(RM(restagt,eval),2); %
    exgre1(resum(1:rK(1))>=lambda_list(1)*(k-1),j1) = 1;
    exgre2(resum(rK(1)+1:rK(1)+rK(2))>=lambda_list(2)*(k-1),j1) = 1;
    exgre3(resum(rK(1)+rK(2)+1:rk)>=lambda_list(3)*(k-1),j1) = 1;%    
end
exgre = [exgre1;exgre2;exgre3];

for j2 = 1:rk
    ob = restagt(j2); %
    if rand(1)<q %
        err1 = rand(1,k); %
        obacts = acts; %
        obacts(err1<err) = ~acts(err1<err); %
        idx = RM(ob,agtlist)*4+obacts*2+exgre(j2,:); %
        exre(j2,:) = str(idx+1); %
    else
        exre(j2,:) = RM(ob,agtlist); %
    end
end




end