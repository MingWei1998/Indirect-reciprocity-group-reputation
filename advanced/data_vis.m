% test and data-visualization
%% 1 reputation
%markers = ['v','^','d','<','x','s','o','>'];
A = table2array(conRe);
h = figure(1);
x = 0.1:0.1:0.9;
for L = 1:8    
    temp = A(:,3*L+1);
    y = temp';
    color = L_color(L,:);
    plot(x,y,'color',color,'LineWidth',2.5,'LineStyle',':','Marker','o','MarkerSize',12.5,'MarkerFaceColor',color);
    hold on
end
box off;
%xlabel('λ');
%ylabel('{r_{L_i}}');
set(gca, 'xTick', [0.1:0.2:0.9]);
set(gca,'XTickLabel',{'','','','',''});
set(gca, 'yTick', [0:0.2:1]);
set(gca,'YTickLabel',{'','','','','',''});
set(gca,'xminortick','on'); %
set(gca,'yminortick','on'); %
ax = gca;
ax.XAxis.MinorTickValues = 0.1:0.1:0.9;
ax.YAxis.MinorTickValues = 0:0.1:1; %
set(gca,'TickLength',[0.018,0.009]); %
set(gca,'linewidth',1); %
set(gca,'XAxisLocation','bottom');  
xlim([0.05,0.95]);
ylim([-0.05,1]);
set(gca, 'LooseInset', [0.005,0.005,0.003,0.02]); %
%legend('L1','L2','L3','L4','L5','L6','L7','L8','Location','northoutside','Orientation','horizontal','boxoff');

%% 1-1 legend
h = figure(1);
%x = [1:8];
x = 1*ones(1,8);
y = 10*ones(1,8);
for L = 1:8
    color = L_color(L,:);
    plot(x(L),y(L),'color',color,'LineWidth',1.5,'LineStyle',':','Marker','o','MarkerSize',3,'MarkerFaceColor',color);
    hold on
end
box off;
set(gca, 'xTick', [0.1:0.2:0.9]);
set(gca,'XTickLabel',{'','','','',''});
set(gca, 'yTick', [0:0.2:1]);
set(gca,'YTickLabel',{'','','','','',''});
xlim([0.05,0.95]);
ylim([-0.05,1]);
legend('                     ','   ','   ','   ','   ','   ','   ','   ','Location','best','Orientation','horizontal');
legend('boxoff');
set(gca, 'LooseInset', [0.005,0.005,0.003,0.02]); %
ax = get(gca);
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off'; %
set(gcf,'unit','centimeters','position',[7 20 20 1]);%


%% 1-1-1 legend double rows
h = figure(1);
%x = [1:8];
x = 1*ones(1,4);
y = 10*ones(1,4);
for L = 1:4
    color = L_color(L+4,:);
    plot(x(L),y(L),'color',color,'LineWidth',2.5,'LineStyle',':','Marker','o','MarkerSize',12.5,'MarkerFaceColor',color);
    hold on
end
box off;
set(gca, 'xTick', [0.1:0.2:0.9]);
set(gca,'XTickLabel',{'','','','',''});
set(gca, 'yTick', [0:0.2:1]);
set(gca,'YTickLabel',{'','','','','',''});
xlim([0.05,0.95]);
ylim([-0.05,1]);
legend('   ','   ','   ','   ','Location','best','Orientation','horizontal');
legend('boxoff');
set(gca, 'LooseInset', [0.005,0.005,0.003,0.02]); %
ax = get(gca);
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off'; %
set(gcf,'unit','centimeters','position',[7 20 8 1]);%

%% 2 cooperation rate of reputation dynamics
%markers = ['v','^','d','<','x','s','o','>'];
x = 0.1:0.1:0.9;
A = cprt;
for t = 1:8
    color = L_color(t,:);
    marker = markers(t);
    plot(x,A(t,:),'color',color,'LineWidth',2.5,'LineStyle',':','Marker','o','MarkerSize',12.5,'MarkerFaceColor',color);
    hold on
end
box off;
%xlabel('λ');
%ylabel('Cooperation Rate');
set(gca, 'xTick', [0.1:0.2:0.9]);
set(gca,'XTickLabel',{'','','','',''});
set(gca, 'yTick', [0.3:0.1:0.7]);
set(gca,'YTickLabel',{'','','','','',''});
set(gca,'xminortick','on'); %
set(gca,'yminortick','on'); %
ax = gca;
ax.XAxis.MinorTickValues = 0.1:0.1:0.9;
ax.YAxis.MinorTickValues = 0.3:0.05:07; %
set(gca,'TickLength',[0.018,0.009]); %
set(gca,'linewidth',1); %
set(gca,'XAxisLocation','bottom');  
xlim([0.05,0.95]);
ylim([0.28,0.7]);
set(gca, 'LooseInset', [0.005,0.005,0.003,0.02]); %
%legend('L1','L2','L3','L4','L5','L6','L7','L8','Location','best');

%% 3
tic
Default;
lambda = 0.7;
lambda_e = lambda; lambda_i = lambda;
PS = 200;
Lnum = Popstc30(PS,1); Cnum = Popstc30(PS,2); Dnum = Popstc30(PS,3);
N = Lnum+Cnum+Dnum;
rep_evol_sym_cr;
% LCR = [Lcr,Lcrspt];
% disp(strmeanpof);
% disp(strmeanpofspt);
% disp(LCR);
toc

%% 4
% %nPopstc = size(Popstc,1);
% N = 30;
% %expgamemb = zeros(nPopstc,27);
% %for i = 1:nPopstc
%     Lnum = 10; Cnum = 12; Dnum = 8;
%     n = [Lnum,Cnum,Dnum];    
%     for k = 2:10
%         exp = zeros(1,3);
%         for l = 0:min([k,Lnum])
%             for c = 0:min([k-l,Cnum])
%                 d = k-c-l;
%                 if d<=Dnum
%                     exp = exp+[l,c,d]*nchoosek(Lnum,l)*nchoosek(Cnum,c)*nchoosek(Dnum,d)/nchoosek(N,k);
%                 end
%             end
%         end
%         disp(exp);
%         %disp(sum(exp));
%         %l = cross(exp,n)
%     end
% %end
% %l = hygepdf()

%% 5
% pi = zeros(1,nk);
% %nc = zeros(1,nk);
% %pk = [2:10]/sum(2:10);
% for k = 2:10
%     %nc(k-1) = ([Lcr(k-1),1,0]*[Lnum;Cnum;Dnum])*k/N;
%     cpof = Lcr(k-1)*(alpha*nc(k-1)-1);
%     dpof = (1-Lcr(k-1))*alpha*nc(k-1);
%     pi(k-1) = (cpof+dpof);    
% end
% %disp(pi);

%% 6
%kGameTime = sum(GameTimespt,2)./[2:10]';
%(repItnum/2)*[2:10]/(N*nk)
%mean(GameTimespt,2)
%pi1 = pk*meanpofspt(:,21);
% k=2;
% cpof = Lcr(k-1)*(alpha*nc(k-1)-1);
% dpof = (1-Lcr(k-1))*alpha*nc(k-1);
% pi = (cpof+dpof);

%% 7 
%markers = ['v','^','d','<','x','s','o','>'];
x = 0.1:0.1:0.9;
%A = table2array(cprt);
A = allAvFr;
color = L_color(8,:);
plot(x,A(:,1),'color',color,'LineWidth',1,'LineStyle','--','Marker','o','MarkerSize',6,'MarkerFaceColor',color);
hold on
plot(x,A(:,2),'color',L_color(9,:),'LineWidth',1,'LineStyle','--','Marker','o','MarkerSize',6,'MarkerFaceColor',L_color(9,:));
hold on
plot(x,A(:,3),'color','k','LineWidth',1,'LineStyle','--','Marker','o','MarkerSize',6,'MarkerFaceColor','k');
xlabel('λ');
ylabel('Strategy Abundance');
xlim([0.05,0.95]);
ylim([-0.1,1.1]);
title('μ=0.0001');
legend('L8','ALLC','ALLD','Location','best');

%% 8
markers = ['v','^','d','<','x','s','o','>'];
x = 0.1:0.1:0.9;
A = huizong;
for t = 1:8
    color = L_color(t,:);
    marker = markers(t);
    plot(x,A(:,t),'Marker',marker,'color',color,'LineWidth',3);
    hold on
end
xlabel('λ');
ylabel('Cooperation Rate');
xlim([0.05,0.95]);
ylim([-0.05,0.5]);
title('μ=0.01');
legend('L1','L2','L3','L4','L5','L6','L7','L8','Location','best');

%% 9 
hold on
xlim = get(gca,'Xlim');
ylim = get(gca,'Ylim');
num = [50.5,100.5];
for nGs = 1:2
    plot(xlim,[num(nGs),num(nGs)],'Color','k','LineStyle',':','LineWidth',1);
    hold on
    plot([num(nGs),num(nGs)],ylim,'Color','k','LineStyle',':','LineWidth',1);
    hold on
end

%% 10 Gamestc30
Gamestc30 = [];
% nPopstc = size(Popstc30,1); %
% N = 30;
for k = 2:10
    for k1 = 0:k
        for k2 = 0:k-k1
            k3 = k-k1-k2;
            vk = [k3,k2,k1];
            Gamestc30 = [Gamestc30;vk];
        end
    end
end

%% 11 PosPopGamestc
nPopstc = size(Popstc30,1); %
nGamestc = size(Gamestc30,1);
PosPopGamestc = zeros(nPopstc,nGamestc);
for nGs = 1:nPopstc
    for nGs = 1:nGamestc
        if Popstc30(nGs,1)>=Gamestc30(nGs,1) && Popstc30(nGs,2)>=Gamestc30(nGs,2) && Popstc30(nGs,3)>=Gamestc30(nGs,3)
            PosPopGamestc(nGs,nGs) = 1;
        end
    end
end

%% 12 GamestcProb
Lstcprob = zeros(nPopstc,nGamestc);
Cstcprob = zeros(nPopstc,nGamestc);
Dstcprob = zeros(nPopstc,nGamestc);
N=30;
for nGs = 1:nPopstc
    for nGs = 1:nGamestc
        if PosPopGamestc(nGs,nGs) == 1
            n1 = Popstc30(nGs,1); n2 = Popstc30(nGs,2); n3 = Popstc30(nGs,3);
            k1 = Gamestc30(nGs,1); k2 = Gamestc30(nGs,2); k3 = Gamestc30(nGs,3);
            k = sum(Gamestc30(nGs,:));
            if k1>=1
                Lstcprob(nGs,nGs) = nchoosek(n1-1,k1-1)*nchoosek(n2,k2)*nchoosek(n3,k3)/nchoosek(N-1,k-1); 
            end
            if k2>=1
                Cstcprob(nGs,nGs) = nchoosek(n1,k1)*nchoosek(n2-1,k2-1)*nchoosek(n3,k3)/nchoosek(N-1,k-1);
            end
            if k3>=1
                Dstcprob(nGs,nGs) = nchoosek(n1,k1)*nchoosek(n2,k2)*nchoosek(n3-1,k3-1)/nchoosek(N-1,k-1);
            end
        end
    end
end
filename = 'GamestcProb.mat';
save(filename,'Lstcprob','Cstcprob','Dstcprob');

%% 13 cooperation rate 2 payoff
tic
Default;
pk = [kmin:kmax]/sum(kmin:kmax);
Dataname = strcat('Indirect-reciprocity-group-reputation\results\formal\strategy evolve\cr2pof\coop_rate\L2\λ=0.1_L2CR.mat');
load(Dataname);
LCR(isnan(LCR)) = 0;
LCRspt(isnan(LCRspt)) = 0;
SMP = cr2pof(Popstc30,Gamestc30,LCRspt,alpha,beta,c,pk,Lstcprob,Cstcprob,Dstcprob);
[PopCoop,AvFr] = slc_mut_equ(mu,s,Popstc30,SMP,LCR);
toc

%% 14 
tic
beta = 2;
s = 1;
Mu = [1,0.5,0.1,0.05,0.01,0.005,0.001,0.0005,0.0001];
for m = 1:9
    mu = Mu(m);
    for strNum = 1:8
        rscoop = [];
        rsAvFr = [];
        for alpha = 0:0.025:0.6
            dataname = strcat('Indirect-reciprocity-group-reputation\results\formal\strategy evolve\cr2pof\μgg0\μ=',num2str(mu),'\β=',num2str(beta),',s=',num2str(s),'\α=',num2str(alpha),'\L',num2str(strNum),'coop&AvFr.mat');
            %dataname = strcat('Indirect-reciprocity-group-reputation\results\formal\strategy evolve\cr2pof\μ→0&SMP\β=',num2str(beta),',s=',num2str(s),'\α=',num2str(alpha),',β=',num2str(beta),',s=',num2str(s),'\SMEq&popcr\L',num2str(strNum),'_SMEq&coop.mat');
            load(dataname);
            rscoop = [allcoop';rscoop];
            rsAvFr = [allAvFr';rsAvFr];
        end
        filename = strcat('Indirect-reciprocity-group-reputation\results\formal\strategy evolve\cr2pof\μgg0\μ=',num2str(mu),'\resort\β=',num2str(beta),'\L',num2str(strNum),'.mat');
        %filename = strcat('Indirect-reciprocity-group-reputation\results\formal\strategy evolve\cr2pof\μ→0&SMP\β=',num2str(beta),',s=',num2str(s),'\resort\L',num2str(strNum),'.mat');
        save(filename,'rscoop','rsAvFr');
    end
end
toc

%% 15 
beta = 2;s = 1;
lw=0.7;
%figure;
for strNum = 8
    figure;
    dataname = strcat('Indirect-reciprocity-group-reputation\results\formal\strategy evolve\cr2pof\μ→0&SMP\β=',num2str(beta),',s=',num2str(s),'\resort\L',num2str(strNum),'.mat');
    load(dataname);
    %sub = subplot(2,4,strNum);
    colorname = strcat('L',num2str(strNum),'_color');
    color_map = eval(colorname);
    imagesc(rscoop(9:25,:));
    colorbar;
    caxis([0,1]);
    cb = colorbar;
    set(cb,'linewidth',lw);
    set(cb,'tickdir','out')  % 
    set(cb,'YTick',[0,0.2,0.4,0.6,0.8,1]); %
    set(cb,'YTickLabel',{'','','','','',''}) %
    colormap(color_map);
    %xlabel('λ'); ylabel('α');
    xticks([1,3,5,7,9]);xticklabels({'','','','',''});
    %yticks([1,5,9,13,17]);yticklabels({'','','','',''});
    %yticks([2,7,12,17,22,27]);yticklabels({'','','','','',''});
    yticks([1,5,9,13,17]);yticklabels({'','','','',''});
    %subtitle(strcat('L',num2str(strNum)));
    %axis xy;
    %xlim(sub,[0.05,0.95]);
    %ylim(sub,[-0.05,1.15]);
%     set(gca, 'xTick', [0.1:0.2:0.9]);
%     set(gca,'XTickLabel',{'','','','',''});
%     set(gca, 'yTick', [0:0.5:2.5]);
%     set(gca,'YTickLabel',{'','','','','',''});
    hold on;
    box off;
    set(gca,'xminortick','on'); %
    set(gca,'yminortick','on'); %
    ax = gca;
    ax.XAxis.MinorTickValues = 1:1:9;
    ax.YAxis.MinorTickValues = 1:1:17; %
    set(gca,'TickLength',[0.018,0.008]); %
    set(gca,'TickDir','out'); %
    set(gca,'linewidth',lw); %
    %
    XL = get(gca,'xlim'); XR = XL(2);
    YL = get(gca,'ylim'); YB = YL(1);
    plot(XL,YB*ones(size(XL)),'color','k','linewidth',lw);
    plot(XR*ones(size(YL)),YL,'color','k','linewidth',lw);
    set(gca,'LooseInset', [0.04,0.025,0.02,0.02]); %
    set(gcf,'unit','centimeters','position',[7 20 4.2 5.5]);%
end

%% 16 
beta = 2;s = 1;
figure;

A = [183,21,64;255,194,64;130,204,221]/255;


for strNum = 8
    dataname = strcat('Indirect-reciprocity-group-reputation\results\formal\strategy evolve\cr2pof\μ→0&SMP\β=',num2str(beta),',s=',num2str(s),'\resort\L',num2str(strNum),'.mat');
    load(dataname);
    rsAvFr(rsAvFr<0) = 0;
    %sub = subplot(3,3,strNum);
    for l = 1:9
        lambda = l*0.1;
        for a = 0:16
            alpha = a*0.025;
            %color_base = rsAvFr(end-a*3-2:end-a*3,l)';
            color_base = rsAvFr(end-a*3-2:end-a*3,l)'/sum([rsAvFr(end-a*3-2:end-a*3,l)]);
            color = color_base*A;
            %color = color/sum(color);
            %fill([lambda-0.05,lambda+0.05,lambda+0.05,lambda-0.05],[alpha-0.5,alpha-0.5,alpha+0.5,alpha+0.5],color,'EdgeAlpha',0);
            %fill([lambda-0.05,lambda+0.05,lambda+0.05,lambda-0.05],[alpha-0.05,alpha-0.05,alpha+0.05,alpha+0.05],color,'EdgeAlpha',0); 
            fill([lambda-0.05,lambda+0.05,lambda+0.05,lambda-0.05],[alpha-0.0125,alpha-0.0125,alpha+0.0125,alpha+0.0125],color,'EdgeAlpha',0);
            hold on
        end
    end
    axis([0.05,0.95,-0.0125,0.4125]);
    %axis square
    %xlabel('λ'); ylabel('α');
    %subtitle(strcat('L',num2str(strNum)));
    set(gca, 'xTick', [0.1:0.2:0.9]);
    set(gca,'XTickLabel',{'','','','',''});
    set(gca, 'yTick', [0:0.1:0.4]);
    set(gca,'YTickLabel',{'','','','',''});
    hold on;
    box off;
    set(gca,'xminortick','on'); %
    set(gca,'yminortick','on'); %
    ax = gca;
    ax.XAxis.MinorTickValues = 0.1:0.1:0.9;
    ax.YAxis.MinorTickValues = 0:0.025:0.4; %
    set(gca,'TickLength',[0.018,0.008]); %
    set(gca,'TickDir','out'); %
    set(gca,'linewidth',lw); %
    %
    XL = get(gca,'xlim'); XR = XL(2);
    YL = get(gca,'ylim'); YT = YL(2);
    plot(XL,YT*ones(size(XL)),'color','k','linewidth',lw);
    plot(XR*ones(size(YL)),YL,'color','k','linewidth',lw);
    set(gca,'LooseInset', [0.04,0.025,0.02,0.02]); %
    set(gcf,'unit','centimeters','position',[7 20 3.2 3.5]);%
end
%sub = subplot(3,3,9);
%% 17 
%figure;
lw=0.85;
l = 1; w = l*sqrt(3)/2; N = 100;
square_list = [];
color_list = [];
for i = 1:N
    for j = 1:N-i+1
        position = [0.5*l*(i-1)+j-1,w*(i-1)];
        square_list = [square_list;position];
        color = [(N-i)*(N-i-j+2)/(N*(N-i+1)),(i-1)/N,(N-i)*(j-1)/(N*(N-i+1))]*A;
        color_list = [color_list;color];
    end
end
K = size(square_list,1);
for i = 1:K
    fill([square_list(i,1),square_list(i,1),square_list(i,1)+l,square_list(i,1)+l],[square_list(i,2),square_list(i,2)+w,square_list(i,2)+w,square_list(i,2)],color_list(i,:),'EdgeAlpha',0); 
    hold on
end
axis equal;
axis([-1,101,-1,88]);
%txt1 = 'Li';txt2 = 'ALLC';txt3 = 'ALLD';
%text(0,0-4,txt1,'HorizontalAlignment','center','FontSize',12);
%text(N/2,N*w+4,txt2,'HorizontalAlignment','center','FontSize',12);
%text(N,0-4,txt3,'HorizontalAlignment','center','FontSize',12);
line([0,N],[0,0],'Linewidth',lw,'Color','black');
line([0,N/2],[0,N*w],'Linewidth',lw,'Color','black');
line([N,N/2],[0,N*w],'Linewidth',lw,'Color','black');
ax = get(gca);
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
set(gca,'LooseInset', [0.002,0.002,0.002,0.002]); %
set(gcf,'unit','centimeters','position',[7 20 3.35 3.5]);%

%% 18 μ>0
beta = 2;
Mu = [1,0.5,0.1,0.05,0.01,0.005,0.001,0.0005,0.0001];
for m = 1:9
    tic
    mu = Mu(m);
    for alpha = 0.275:0.025:0.6
        mugg0;
    end
    toc
end

%% 19 
lw=0.7;
beta = 1;
s=0.3;
Slist = [0.3,0.1,0.03,0.01,0.003,0.001];
Mu = [1,0.5,0.1,0.05,0.01,0.005,0.001,0.0005,0.0001];
%for i = 1:6
    %s = Slist(i);
    mu = 0.05;
    figure;
    for strNum = 4
        %dataname = strcat('Indirect-reciprocity-group-reputation\results\formal\strategy evolve\cr2pof\μgg0\μ=',num2str(mu),'\resort\β=',num2str(beta),'\L',num2str(strNum),'.mat');
        dataname = strcat('Indirect-reciprocity-group-reputation\results\formal\strategy evolve\cr2pof\s_changing\s=',num2str(s),'\resort\L',num2str(strNum),'.mat');
        load(dataname);
        %sub = subplot(2,4,strNum);
        colorname = strcat('L',num2str(strNum),'_color');
        color_map = eval(colorname);
        %imagesc(rscoop(11:27,:));
        imagesc(rscoop);
        cb = colorbar;
        caxis([0,1]);
        set(cb,'linewidth',lw);
        set(cb,'tickdir','out')  %
        set(cb,'YTick',[0,0.2,0.4,0.6,0.8,1]); %
        set(cb,'YTickLabel',{'','','','','',''}) %
        colormap(color_map);
        %xlabel('λ'); ylabel('α');
        xticks([1,3,5,7,9]);xticklabels({'','','','',''});
        %yticks([1,3,5,7,9,11,13,15,17]);yticklabels({16,14,12,10,8,6,4,2,0});
        yticks([2,7,12,17,22,27]);yticklabels({'','','','','',''});
        %yticks([1,3,5,7,9,11,13,15,17,19,21,23,25]);yticklabels({0.6,0.55,0.5,0.45,0.4,0.35,0.3,0.25,0.2,0.15,0.1,0.05,0});
        %subtitle(strcat('L',num2str(strNum)));
        %axis xy;
        %xlim(sub,[0.05,0.95]);
        %ylim(sub,[-0.05,1.15]);
    %     set(gca, 'xTick', [0.1:0.2:0.9]);
    %     set(gca,'XTickLabel',{'','','','',''});
    %     set(gca, 'yTick', [0:0.5:2.5]);
    %     set(gca,'YTickLabel',{'','','','','',''});
        hold on;
        box off;
        set(gca,'xminortick','on'); %
        set(gca,'yminortick','on'); %
        ax = gca;
        ax.XAxis.MinorTickValues = 1:1:9;
        ax.YAxis.MinorTickValues = 1:1:27; %
        set(gca,'TickLength',[0.018,0.008]); %
        set(gca,'TickDir','out'); %
        set(gca,'linewidth',lw); %
        %
        XL = get(gca,'xlim'); XR = XL(2);
        YL = get(gca,'ylim'); YB = YL(1);
        plot(XL,YB*ones(size(XL)),'color','k','linewidth',lw);
        plot(XR*ones(size(YL)),YL,'color','k','linewidth',lw);
        set(gca,'LooseInset', [0.04,0.025,0.02,0.02]); %
        set(gcf,'unit','centimeters','position',[7 20 4.2 4]);%
    end
    %filename = strcat('Indirect-reciprocity-group-reputation\results\formal\strategy evolve\cr2pof\μgg0\μ=',num2str(mu),'\resort\β=',num2str(beta),'\cooprate.fig');
    %filename = strcat('Indirect-reciprocity-group-reputation\results\formal\strategy evolve\cr2pof\s_changing\s=',num2str(s),'\resort\cooprate.fig');
    %savefig(filename);
    %hold off
%end

%% 20 
beta = 1;s = 0.3;
A = [183,21,64;255,194,64;130,204,221]/255;
%for m = 9
    mu = 0.05;
    figure;
    for strNum = 8
        %dataname = strcat('Indirect-reciprocity-group-reputation\results\formal\strategy evolve\cr2pof\μgg0\μ=',num2str(mu),'\resort\β=',num2str(beta),'\L',num2str(strNum),'.mat');
        dataname = strcat('Indirect-reciprocity-group-reputation\results\formal\strategy evolve\cr2pof\s_changing\s=',num2str(s),'\resort\L',num2str(strNum),'.mat');
        load(dataname);
        %sub = subplot(2,4,strNum);
        for l = 1:9
            lambda = l*0.1;
            for a = 0:26
                alpha = a*0.1;
                color_base = rsAvFr(end-a*3-2:end-a*3,l)';
                color = color_base*A;
                %fill([lambda-0.05,lambda+0.05,lambda+0.05,lambda-0.05],[alpha-0.5,alpha-0.5,alpha+0.5,alpha+0.5],color,'EdgeAlpha',0);
                fill([lambda-0.05,lambda+0.05,lambda+0.05,lambda-0.05],[alpha-0.05,alpha-0.05,alpha+0.05,alpha+0.05],color,'EdgeAlpha',0);
                %fill([lambda-0.05,lambda+0.05,lambda+0.05,lambda-0.05],[alpha-0.0125,alpha-0.0125,alpha+0.0125,alpha+0.0125],color,'EdgeAlpha',0);
                hold on
            end
        end
        %axis([0.05,0.95,-0.5,16.5]);
        axis([0.05,0.95,-0.05,2.65]);
        %axis([0.05,0.95,-0.0125,0.6125]);
        %axis square
        %xlabel('λ'); ylabel('α');
        %subtitle(strcat('L',num2str(strNum)));
        %axis square
        set(gca, 'xTick', [0.1:0.2:0.9]);
        set(gca,'XTickLabel',{'','','','',''});
        set(gca, 'yTick', [0:0.5:2.5]);
        set(gca,'YTickLabel',{'','','','','',''});
        hold on;
        box off;
        set(gca,'xminortick','on'); %
        set(gca,'yminortick','on'); %
        ax = gca;
        ax.XAxis.MinorTickValues = 0.1:0.1:0.9;
        ax.YAxis.MinorTickValues = 0:0.1:2.6; %
        set(gca,'TickLength',[0.018,0.008]); %
        set(gca,'TickDir','out'); %
        set(gca,'linewidth',lw); %
        %
        XL = get(gca,'xlim'); XR = XL(2);
        YL = get(gca,'ylim'); YT = YL(2);
        plot(XL,YT*ones(size(XL)),'color','k','linewidth',lw);
        plot(XR*ones(size(YL)),YL,'color','k','linewidth',lw);
        set(gca,'LooseInset', [0.04,0.025,0.02,0.02]); %
        set(gcf,'unit','centimeters','position',[7 20 3.2 3.5]);%
    end
%     filename = strcat('Indirect-reciprocity-group-reputation\results\formal\strategy evolve\cr2pof\μgg0\μ=',num2str(mu),'\resort\β=',num2str(beta),'\AvFr.fig');
%     savefig(filename);
%     hold off
%end

%% 21-1 s resort
beta = 1;
Slist = [0.3,0.1,0.03,0.01,0.003,0.001];
for i = 1:6
    s = Slist(i);
    for strNum = 1:8
        rscoop = [];
        rsAvFr = [];
        for alpha = 0:0.1:2.6
            dataname = strcat('Indirect-reciprocity-group-reputation\results\formal\strategy evolve\cr2pof\s_changing\s=',num2str(s),'\α=',num2str(alpha),'\SMEq&popcr\L',num2str(strNum),'_SMEq&coop.mat');
            load(dataname);
            rscoop = [allcoop';rscoop];
            rsAvFr = [allSMEq';rsAvFr];
        end
        filename = strcat('Indirect-reciprocity-group-reputation\results\formal\strategy evolve\cr2pof\s_changing\s=',num2str(s),'\resort\L',num2str(strNum),'.mat');
        %filename = strcat('Indirect-reciprocity-group-reputation\results\formal\strategy evolve\cr2pof\μ→0&SMP\β=',num2str(beta),',s=',num2str(s),'\resort\L',num2str(strNum),'.mat');
        save(filename,'rscoop','rsAvFr');
    end
end
%% 21-2 
Slist = [0.3,0.1,0.03,0.01,0.003,0.001];
for i = 1:6
    s = Slist(i);
    figure;
%     for strNum = 1:8
%         dataname = strcat('Indirect-reciprocity-group-reputation\results\formal\strategy evolve\cr2pof\s_changing\s=',num2str(s),'\resort\L',num2str(strNum),'.mat');
%         load(dataname);
%         sub = subplot(2,4,strNum);
%         colorname = strcat('L',num2str(strNum),'_color');
%         color_map = eval(colorname);
%         %imagesc(rscoop(11:27,:));
%         imagesc(rscoop);
%         colorbar;
%         caxis([0,1]);
%         colormap(sub,color_map);
%         xlabel('λ'); ylabel('α');
%         xticks([2,4,6,8]);xticklabels({0.2,0.4,0.6,0.8});
%         %yticks([1,3,5,7,9,11,13,15,17]);yticklabels({16,14,12,10,8,6,4,2,0});
%         yticks([1,3,5,7,9,11,13,15,17,19,21,23,25,27]);yticklabels({2.6,2.4,2.2,2.0,1.8,1.6,1.4,1.2,1.0,0.8,0.6,0.4,0.2,0});
%         %yticks([1,3,5,7,9,11,13,15,17,19,21,23,25]);yticklabels({0.6,0.55,0.5,0.45,0.4,0.35,0.3,0.25,0.2,0.15,0.1,0.05,0});
%         %axis([0.05,0.95,-0.5,16.5]);
%         axis([0.05,0.95,-0.05,2.65]);
%         %axis([0.05,0.95,-0.0125,0.6125]);
%         subtitle(strcat('L',num2str(strNum)));
%         %axis xy;
%     end
%     filename = strcat('E:\Lab\北航研究生\process\Results\formal\strategy evolve\cr2pof\s_changing\s=',num2str(s),'\resort\cooprate.fig');
%     savefig(filename);
%     hold off

    for strNum = 1:8
        dataname = strcat('Indirect-reciprocity-group-reputation\results\formal\strategy evolve\cr2pof\s_changing\s=',num2str(s),'\resort\L',num2str(strNum),'.mat');
        load(dataname);
        sub = subplot(2,4,strNum);
        for l = 1:9
            lambda = l*0.1;
            for a = 0:26
                alpha = a*0.1;
                color = rsAvFr(end-a*3-2:end-a*3,l)';
                %fill([lambda-0.05,lambda+0.05,lambda+0.05,lambda-0.05],[alpha-0.5,alpha-0.5,alpha+0.5,alpha+0.5],color,'EdgeAlpha',0);
                fill([lambda-0.05,lambda+0.05,lambda+0.05,lambda-0.05],[alpha-0.05,alpha-0.05,alpha+0.05,alpha+0.05],color,'EdgeAlpha',0);
                %fill([lambda-0.05,lambda+0.05,lambda+0.05,lambda-0.05],[alpha-0.0125,alpha-0.0125,alpha+0.0125,alpha+0.0125],color,'EdgeAlpha',0);
                hold on
            end
        end
        %axis([0.05,0.95,-0.5,16.5]);
        axis([0.05,0.95,-0.05,2.65]);
        %axis([0.05,0.95,-0.0125,0.6125]);
        %axis square
        xlabel('λ'); ylabel('α');
        subtitle(strcat('L',num2str(strNum)));
    end
    filename = strcat('Indirect-reciprocity-group-reputation\results\formal\strategy evolve\cr2pof\s_changing\s=',num2str(s),'\resort\AvFr.fig');
    savefig(filename);
    hold off

end

%% 22 fixation prob
alpha = 0.9;
beta = 1; s = 1;
Lambda = 0.1:0.1:0.9;
xlim = [0.05,0.95];
lw=1.2;ms=4;
figure;
for strNum = 1
    dataname = strcat('Indirect-reciprocity-group-reputation\results\formal\strategy evolve\cr2pof\μ→0&SMP\β=',num2str(beta),',s=',num2str(s),'\α=',num2str(alpha),',β=',num2str(beta),',s=',num2str(s),'\fix_prob\L',num2str(strNum),'fixprob.mat');
    load(dataname);
    %sub = subplot(2,4,strNum);
    axis([0.05,0.95,-0.02,0.3]);
    hold on
    plot(xlim,[1/30,1/30],'Color',[0.75,0.75,0.75],'LineStyle','-','LineWidth',1);
    plot(Lambda,FXP(:,4),'Marker','^','Color','k','MarkerSize',ms,'Linewidth',lw,'MarkerFaceColor','k');
    plot(Lambda,FXP(:,2),'Marker','<','Color','k','MarkerSize',ms,'Linewidth',lw,'MarkerFaceColor','k','LineStyle',':');
    plot(Lambda,FXP(:,1),'Marker','square','Color',L_color(strNum,:),'MarkerSize',ms,'Linewidth',lw,'MarkerFaceColor',L_color(strNum,:),'LineStyle',':');    
    plot(Lambda,FXP(:,3),'Marker','o','Color',L_color(strNum,:),'MarkerSize',ms,'Linewidth',lw,'MarkerFaceColor',L_color(strNum,:));        
    %legend('',strcat('L',num2str(strNum),'→ALLD'),strcat('L',num2str(strNum),'→ALLC'),strcat('ALLC→','L',num2str(strNum)),strcat('ALLD→','L',num2str(strNum)),'Location','northwest');
    %xlabel('λ'); ylabel('fixation probability');
    %subtitle(strcat('L',num2str(strNum)));
    set(gca, 'xTick', [0.1:0.2:0.9]);
    set(gca,'XTickLabel',{'','','','',''});
    set(gca, 'yTick', [0:0.1:0.3]);
    set(gca,'YTickLabel',{'','','',''});
    set(gca,'xminortick','on'); %
    set(gca,'yminortick','on'); %
    ax = gca;
    ax.XAxis.MinorTickValues = 0.1:0.1:0.9;
    ax.YAxis.MinorTickValues = 0:0.05:0.3; %
    set(gca,'TickLength',[0.018,0.009]); %
    set(gca,'linewidth',0.5); %
    set(gca,'XAxisLocation','bottom');
    set(gca,'TickDir','out'); %
    set(gca, 'LooseInset', [0.02,0.02,0.003,0.02]); %
    set(gcf,'unit','centimeters','position',[7 20 2.7 3.5]);%
end

%% 22-1 
h = figure(1);
lw=1.2;ms=4;
%x = [1:8];
x = 1*ones(1,4);
y = 10*ones(1,4);
hold on;
comcolor = [0.7,0.7,0.7];
%plot(1,10,'Marker','^','Color','k','MarkerSize',ms,'Linewidth',lw,'MarkerFaceColor','k');
%plot(1,10,'Marker','<','Color','k','MarkerSize',ms,'Linewidth',lw,'MarkerFaceColor','k','LineStyle',':');
%plot(1,10,'Marker','square','Color',comcolor,'MarkerSize',ms,'Linewidth',lw,'MarkerFaceColor',comcolor,'LineStyle',':');    
plot(1,10,'Marker','o','Color',comcolor,'MarkerSize',ms,'Linewidth',lw,'MarkerFaceColor',comcolor);
box off;
set(gca, 'xTick', [0.1:0.2:0.9]);
set(gca,'XTickLabel',{'','','','',''});
set(gca, 'yTick', [0:0.2:1]);
set(gca,'YTickLabel',{'','','','','',''});
%xlim([0.05,0.95]);
ylim([-0.05,1]);
Leg = legend('   ','Location','north','Orientation','vertical');
legend('boxoff');
Leg.ItemTokenSize = [18,1];% 
set(gca, 'LooseInset', [0.005,0.005,0.003,0.02]); %
ax = get(gca);
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off'; %
set(gcf,'unit','centimeters','position',[7 20 0.7 0.5]);%

%% 23 
thr = 1/30;
beta = 1;s = 1;
figure;
dimgrey = [0.41176,0.41176,0.41176];lightgrey = [0.82745,0.82745,0.82745];black = [0,0,0];
for strNum = 1:8
    sub = subplot(2,4,strNum);
    hold on
    for alpha = 0:0.1:2.6
        a = alpha*10;
        %dataname = strcat('Indirect-reciprocity-group-reputation\results\formal\strategy evolve\cr2pof\μ→0&SMP\β=',num2str(beta),',s=',num2str(s),'\α=',num2str(alpha),',β=',num2str(beta),',s=',num2str(s),'\fix_prob\L',num2str(strNum),'fixprob.mat');
        for l = 1:9
            lambda = l*0.1;
            dataname = strcat('Indirect-reciprocity-group-reputation\results\formal\strategy evolve\cr2pof\μ→0&SMP\β=',num2str(beta),',s=1\α=',num2str(alpha),',β=',num2str(beta),',s=1\SMP\L',num2str(strNum),'_λ=',num2str(lambda),'_SMP.mat');
            load(dataname);            
            if SMP(2,1)<=SMP(2,2) && SMP(32,1)>SMP(32,3) %
                color = lightgrey;
            elseif SMP(2,1)>SMP(2,2) && SMP(32,1)<=SMP(32,3) %
                color = dimgrey;
            elseif SMP(2,1)<=SMP(2,2) && SMP(32,1)<=SMP(32,3) %
                color = black;
            else
                color = L_color(strNum,:);
            end            
            fill([lambda-0.05,lambda+0.05,lambda+0.05,lambda-0.05],[alpha-0.05,alpha-0.05,alpha+0.05,alpha+0.05],color,'EdgeAlpha',0);
        end
    end
    axis([0.05,0.95,-0.05,2.65]);
    xlabel('λ'); ylabel('α');
    subtitle(strcat('L',num2str(strNum)));
end
%legend()

%% 24
figure;
A = [1,0,0;0.196,0.803,0.196;0.2549,0.41176,0.88235];
l = 1; w = l*sqrt(3)/2; N = 50;
square_list = [];
color_list = [];
for i = 1:N
    for j = 1:N-i+1
        position = [0.5*l*(i-1)+j-1,w*(i-1)];
        square_list = [square_list;position];
        color = [(N-i)*(N-i-j+2)/(N*(N-i+1)),(i-1)/N,(N-i)*(j-1)/(N*(N-i+1))]*A;
        color_list = [color_list;color];
    end
end
K = size(square_list,1);
for i = 1:K
    fill([square_list(i,1),square_list(i,1),square_list(i,1)+l,square_list(i,1)+l],[square_list(i,2),square_list(i,2)+w,square_list(i,2)+w,square_list(i,2)],color_list(i,:),'EdgeAlpha',0); 
    hold on
end

%% 25 
for lambda = 0.1:0.1:0.9
    crpt = [];
    figname = strcat('Indirect-reciprocity-group-reputation\results\formal\reputation evolve\k=rand, sym\cooperation rate data\',num2str(lambda),'.fig');
    uiopen(figname,1);
    hl = get(gca,'Children');
    for strNum = 1:8
        crpt = [hl(strNum).YData;crpt];
    end
    filename = strcat('Indirect-reciprocity-group-reputation\results\formal\reputation evolve\k=rand, sym\cooperation rate data\fig2data\',num2str(lambda),'.mat');
    save(filename,"crpt");
    close
end
cprt = [];
figname = strcat('Indirect-reciprocity-group-reputation\results\formal\reputation evolve\k=rand, sym\cooperation rate data\cprt.fig');
uiopen(figname,1);
hl = get(gca,'Children');
for strNum = 1:8
    cprt = [hl(strNum).YData;cprt];
end
filename = strcat('Indirect-reciprocity-group-reputation\results\formal\reputation evolve\k=rand, sym\cooperation rate data\fig2data\cprt.mat');
save(filename,"cprt");
close

%% 26 
beta = 1;alpha = 0.5;s = 1; %
strNum = 5; popNum = 333; %
Lcolor = L_color(strNum,:);%
dimgrey = [0.51176,0.51176,0.51176];lightgrey = [0.82745,0.82745,0.82745];black = [0,0,0];
Gdata = []; C2pdata = [];
Lambda = 0.1:0.1:0.9;
for l = 1:9
    lambda = l*0.1;
    gddata = strcat('Indirect-reciprocity-group-reputation\results\formal\strategy evolve\sim_payoff\L',num2str(strNum),'\λ=',num2str(lambda),'.mat');
    load(gddata); SMP(isnan(SMP)) = 0; %
    Gdata = [Gdata;SMP(popNum,:)];
    c2pdata = strcat('Indirect-reciprocity-group-reputation\results\formal\strategy evolve\cr2pof\μ→0&SMP\β=',num2str(beta),',s=',num2str(s),'\α=',num2str(alpha),',β=',num2str(beta),',s=',num2str(s),'\SMP\L',num2str(strNum),'_λ=',num2str(lambda),'_SMP.mat');
    load(c2pdata);
    C2pdata = [C2pdata;SMP(popNum,:)];
end
figure(1);
hold on

plot(Lambda,Gdata(:,1),'Color',Lcolor,'LineStyle','none','Marker','^','MarkerFaceColor',Lcolor,'MarkerSize',6);
plot(Lambda,Gdata(:,2),'Color',lightgrey,'LineStyle','none','Marker','<','MarkerFaceColor',lightgrey,'MarkerSize',6);
plot(Lambda,Gdata(:,3),'Color',black,'LineStyle','none','Marker','>','MarkerFaceColor',black,'MarkerSize',6);

plot(Lambda,C2pdata(:,1),'Color',Lcolor,'LineStyle','--','Marker','.','MarkerFaceColor',Lcolor,'LineWidth',2);
plot(Lambda,C2pdata(:,2),'Color',lightgrey,'LineStyle','--','Marker','.','MarkerFaceColor',lightgrey,'LineWidth',2);
plot(Lambda,C2pdata(:,3),'Color',black,'LineStyle','--','Marker','.','MarkerFaceColor',black,'LineWidth',2);

%legend(strcat('L',num2str(strNum)),'ALLC','ALLD',strcat('c2pL',num2str(strNum)),'c2pALLC','c2pALLD');
%xlabel('λ'); ylabel('payoff');
set(gca, 'xTick', [0.1:0.2:0.9]);
set(gca,'XTickLabel',{'','','','',''});
set(gca, 'yTick', [0.2:0.5:1.7]);
set(gca,'YTickLabel',{'','','',''});
set(gca,'xminortick','on'); %
set(gca,'yminortick','on'); %
ax = gca;
ax.XAxis.MinorTickValues = 0.1:0.1:0.9;
ax.YAxis.MinorTickValues = 0.2:0.125:1.7; %
set(gca,'TickLength',[0.018,0.009]); %
set(gca,'linewidth',1); %
set(gca,'XAxisLocation','bottom');  
xlim([0.05,0.95]);
ylim([0.2,1.82]);
set(gca, 'LooseInset', [0.005,0.005,0.003,0.02]); %
set(gcf,'unit','centimeters','position',[7 20 5 6]);%

%% 27 legend
orange = [1,0.5,0];
figure;
hold on

plot(Lambda,Gdata(:,1),'Color',dimgrey,'LineStyle','none','Marker','^','MarkerFaceColor',dimgrey,'MarkerSize',6);
plot(Lambda,Gdata(:,2),'Color',lightgrey,'LineStyle','none','Marker','<','MarkerFaceColor',lightgrey,'MarkerSize',6);
plot(Lambda,Gdata(:,3),'Color',black,'LineStyle','none','Marker','>','MarkerFaceColor',black,'MarkerSize',6);

plot(Lambda,C2pdata(:,1),'Color',dimgrey,'LineStyle','--','Marker','.','MarkerFaceColor',dimgrey,'LineWidth',2);
plot(Lambda,C2pdata(:,2),'Color',lightgrey,'LineStyle','--','Marker','.','MarkerFaceColor',lightgrey,'LineWidth',2);
plot(Lambda,C2pdata(:,3),'Color',black,'LineStyle','--','Marker','.','MarkerFaceColor',black,'LineWidth',2);

box off;
set(gca, 'xTick', [0.1:0.2:0.9]);
set(gca,'XTickLabel',{'','','','',''});
set(gca, 'yTick', [0:0.2:1]);
set(gca,'YTickLabel',{'','','','','',''});
xlim([5,9]);
ylim([10,11]);
legend('   ','   ','   ','   ','   ','   ','Location','north','Orientation','horizontal');
legend('boxoff');
set(gca, 'LooseInset', [0.005,0.005,0.003,0.02]); %
ax = get(gca);
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off'; %
set(gcf,'unit','centimeters','position',[7 20 20 1]);%


%% 28 
beta = 2;s = 1;mu=0.05;
figure;

A = [183,21,64;255,194,64;130,204,221]/255;

for strNum = 1:8
    dataname = strcat('Indirect-reciprocity-group-reputation\results\formal\strategy evolve\cr2pof\μ→0&SMP\β=2,s=1\resort\L',num2str(strNum),'.mat');
    load(dataname);
    rsAvFr(rsAvFr<0) = 0;
    sub = subplot(3,3,strNum);
    for l = 1:9
        lambda = l*0.1;
        for a = 0:24
            alpha = a*0.025;
            color_base = rsAvFr(end-a*3-2:end-a*3,l)'/sum(rsAvFr(end-a*3-2:end-a*3,l));
            color = color_base*A;
            %color = color/sum(color);
            %fill([lambda-0.05,lambda+0.05,lambda+0.05,lambda-0.05],[alpha-0.5,alpha-0.5,alpha+0.5,alpha+0.5],color,'EdgeAlpha',0);
            %fill([lambda-0.05,lambda+0.05,lambda+0.05,lambda-0.05],[alpha-0.05,alpha-0.05,alpha+0.05,alpha+0.05],color,'EdgeAlpha',0); 
            fill([lambda-0.05,lambda+0.05,lambda+0.05,lambda-0.05],[alpha-0.0125,alpha-0.0125,alpha+0.0125,alpha+0.0125],color,'EdgeAlpha',0);
            hold on
        end
    end
    %axis([0.05,0.95,-0.5,16.5]);
    %axis([0.05,0.95,-0.05,2.65]);
    axis([0.05,0.95,-0.0125,0.6125]);
    %axis square
    xlabel('λ'); ylabel('α');
    subtitle(strcat('L',num2str(strNum)));
end
sub = subplot(3,3,9);
%%% 25-3 triangle legend
%figure;
l = 1; w = l*sqrt(3)/2; N = 100;
square_list = [];
color_list = [];
for i = 1:N
    for j = 1:N-i+1
        position = [0.5*l*(i-1)+j-1,w*(i-1)];
        square_list = [square_list;position];
        color = [(N-i)*(N-i-j+2)/(N*(N-i+1)),(i-1)/N,(N-i)*(j-1)/(N*(N-i+1))]*A;
        color_list = [color_list;color];
    end
end
K = size(square_list,1);
for i = 1:K
    fill([square_list(i,1),square_list(i,1),square_list(i,1)+l,square_list(i,1)+l],[square_list(i,2),square_list(i,2)+w,square_list(i,2)+w,square_list(i,2)],color_list(i,:),'EdgeAlpha',0); 
    hold on
end
axis equal;
axis([-10,110,-5,95]);
txt1 = 'Ld';txt2 = 'ALLC';txt3 = 'ALLD';
text(0,0-4,txt1,'HorizontalAlignment','center','FontSize',12);
text(N/2,N*w+4,txt2,'HorizontalAlignment','center','FontSize',12);
text(N,0-4,txt3,'HorizontalAlignment','center','FontSize',12);
line([0,N],[0,0],'Linewidth',2,'Color','black');
line([0,N/2],[0,N*w],'Linewidth',2,'Color','black');
line([N,N/2],[0,N*w],'Linewidth',2,'Color','black');
ax = get(gca);
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';










