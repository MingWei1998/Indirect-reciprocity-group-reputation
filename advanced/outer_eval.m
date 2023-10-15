% 
function[exre] = outer_eval(agtlist,RM,para,acts,Lpop,str,k,q,err)
restagt = setdiff(Lpop,agtlist); %
koL = length(restagt);
exgre = zeros(koL,k); exre = zeros(koL,k);%

for j = 1:k %
    foc = agtlist(j); %
    eval = agtlist(agtlist~=foc); %
    resum = sum(RM(restagt,eval),2); %
    exgre(resum>=para*(k-1),j) = 1;%
end

for j = 1:koL %
    ob = restagt(j); %
    if rand(1)<q %
        err1 = rand(1,k); %
        obacts = acts; %
        obacts(err1<err) = ~acts(err1<err); %
        idx = RM(ob,agtlist)*4+obacts*2+exgre(j,:); %
        exre(j,:) = str(idx+1); %
    else
        exre(j,:) = RM(ob,agtlist); %
    end
end
end