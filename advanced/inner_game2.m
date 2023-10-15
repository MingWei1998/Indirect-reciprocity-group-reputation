% 群内博弈过程
function[acts,inre] = inner_game2(agtlist,RM,para,Lpop,Cpop,Dpop,str,k)
%vtagt = agtlist(agtlist<=Lnum); 
Lagt = intersect(agtlist,Lpop); kL = length(Lagt); 
Cagt = intersect(agtlist,Cpop); kC = length(Cagt); 
Dagt = intersect(agtlist,Dpop); kD = length(Dagt); 
acts = [ones(1,kL+kC),zeros(1,kD)];ingre = zeros(kL,kL); inre = zeros(kL,k);
%acts(agtlist<=Lnum+Cnum & agtlist>Lnum) = 1; 
%acts(agtlist<=Lnum+Cnum+Dnum & agtlist>Lnum+Cnum) = 0; 

for j = 1:kL 
    foc = Lagt(j);
    eval = agtlist(agtlist ~= foc); 
    resum = sum(RM(Lagt,eval),2); 
    ingre(resum>=para*(k-1),j) = 1;
    idx = RM(foc,foc)*2 + ingre(j,j);
    acts(j) = str(idx+9); 
end

for j = 1:kL 
    foc = Lagt(j);
    idx1 = RM(foc,agtlist)*4+acts*2+ingre(j,:); 
    inre(j,:) = str(idx1+1);
end
end