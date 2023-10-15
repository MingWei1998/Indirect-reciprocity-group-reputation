% 
function[acts,inre] = lainner_game(agtlist,RM,lambda_list,K,str,k)
% 
agtla = [lambda_list(1)*ones(1,K(1)),lambda_list(2)*ones(1,K(2)),lambda_list(3)*ones(1,K(3))];
% init
acts = ones(1,k); %
ingre = zeros(1,k); % 
inre = zeros(k,k); %

for j = 1:k %
    foc = agtlist(j); %
    eval = agtlist(agtlist ~= foc); %  
    %
    if sum(RM(foc,eval)) >= agtla(j)*(k-1) 
        ingre(j) = 1;
    else
        ingre(j) = 0;
    end
    idx = RM(foc,foc)*2 + ingre(j); %
    acts(j) = str(idx+9); %
end

for j = 1:k % 
    foc = agtlist(j); %
    idx1 = RM(foc,agtlist)*4+acts*2+ingre(j); %
    inre(j,:) = str(idx1+1);
end

end