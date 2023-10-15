%% fixation probability
function[fxp] = fixprob(s,str1,str2,Pop,SMP)
fxp = zeros(1,2);
if str1==2 && str2==3
    ind = Pop(:,1)==0 & Pop(:,2)~=0 & Pop(:,3)~=0;
    nSMP = SMP(ind,[str1,str2]);
elseif str1==1 && str2==2
    ind = Pop(:,1)~=0 & Pop(:,2)~=0 & Pop(:,3)==0;
    nSMP = SMP(ind,[str1,str2]);
elseif str1==1 && str2==3
    ind = Pop(:,1)~=0 & Pop(:,2)==0 & Pop(:,3)~=0;
    nSMP = SMP(ind,[str1,str2]);
end
payoff1 = nSMP(:,1);
payoff2 = nSMP(:,2);
z = exp(-s*(payoff1(end:-1:1)-payoff2(end:-1:1)));
fxp(1) = 1/(sum(cumprod(z))+1);
z = exp(-s*(payoff2-payoff1));
fxp(2) = 1/(sum(cumprod(z))+1);
end