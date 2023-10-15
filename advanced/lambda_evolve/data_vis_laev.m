% data_vis_laev
%% 1 
beta = 1;
Default;
for alpha = 2:0.1:2.6
    tic
    laAvFr;
    toc
end

%% 2 
beta = 1;
for strNum = 1:8
    rsAvFr = [];
    rscoop = [];
    for alpha = 0:0.1:2.6
        dataname = strcat('Indirect-reciprocity-group-reputation\results\formal\lambda evolve\μ→0&SMP\β=',num2str(beta),'\α=',num2str(alpha),'\SMEq&coop.mat');
        load(dataname);
        rsAvFr = [rsAvFr;allSMEq(strNum,:)];
        rscoop = [rscoop,allcoop(strNum)];        
    end
    filename = strcat('Indirect-reciprocity-group-reputation\results\formal\lambda evolve\μ→0&SMP\β=',num2str(beta),'\resort\L',num2str(strNum),'.mat');
    %filename = strcat('Indirect-reciprocity-group-reputation\results\formal\strategy evolve\cr2pof\μ→0&SMP\β=',num2str(beta),',s=',num2str(s),'\resort\L',num2str(strNum),'.mat');
    save(filename,'rscoop','rsAvFr');
end

%% 3-1 
beta=1;
%markers = ['v','^','d','<','x','s','o','>'];
mu=0.05;
figure;
Alpha = 0:0.1:2.6;
for strNum = 1:8
    %dataname = strcat('Indirect-reciprocity-group-reputation\results\formal\lambda evolve\μ→0&SMP\β=',num2str(beta),'\resort\L',num2str(strNum),'.mat');
    dataname = strcat('Indirect-reciprocity-group-reputation\results\formal\lambda evolve\μgg0\μ=',num2str(mu),'\β=',num2str(beta),'\L',num2str(strNum),'coop&AvFr.mat');
    load(dataname);
    color = L_color(strNum,:);
    %marker = markers(strNum);
    %plot(Alpha,rscoop,'color',color,'LineWidth',1,'LineStyle',':','Marker','o','MarkerSize',5.8,'MarkerFaceColor',color);
    plot(Alpha,allcoop,'color',color,'LineWidth',1,'LineStyle',':','Marker','o','MarkerSize',5.8,'MarkerFaceColor',color);
    hold on
end
%xlabel('α');
%ylabel('Cooperation Rate');
box off;
set(gca, 'xTick', [0:0.5:2.5]);
set(gca,'XTickLabel',{'','','','','',''});
set(gca, 'yTick', [0:0.2:1]);
set(gca,'YTickLabel',{'','','','','',''});
set(gca,'xminortick','on'); %
set(gca,'yminortick','on'); %
ax = gca;
ax.XAxis.MinorTickValues = 0:0.1:2.6;
ax.YAxis.MinorTickValues = 0:0.1:1; %
set(gca,'TickLength',[0.018,0.009]); %
set(gca,'linewidth',0.5); %
set(gca,'XAxisLocation','bottom');
xlim([-0.12,2.65]);
ylim([-0.06,1]);
set(gca, 'LooseInset', [0.005,0.005,0.003,0.02]); %
%title('μ=0.01');
%legend('L1','L2','L3','L4','L5','L6','L7','L8','Location','best');
set(gcf,'unit','centimeters','position',[7 20 6.93 5.2]);%

%% 3-2 
beta=1;
lw=0.7;
figure;
Alpha = 0:0.1:2.6;
for strNum = 8
    %dataname = strcat('Indirect-reciprocity-group-reputation\results\formal\lambda evolve\μ→0&SMP\β=',num2str(beta),'\resort\L',num2str(strNum),'.mat');
    dataname = strcat('Indirect-reciprocity-group-reputation\results\formal\lambda evolve\μgg0\μ=',num2str(mu),'\β=',num2str(beta),'\L',num2str(strNum),'coop&AvFr.mat');
    load(dataname);
    %sub = subplot(2,4,strNum);
    %X = rsAvFr;
    %bar(Alpha,rsAvFr,'stacked');
    bar(Alpha,allAvFr,'stacked');
    %xlabel('α');
    %ylabel('Parameter Abundance');
    %xticks([0,5,10,15,20,25]);xticklabels({0,0.5,1,1.5,2,2.5});
    %xlim([0,2.7]);
    ylim([0,1]);
    %subtitle(strcat('L',num2str(strNum)));
    %legend('λ=0.1','λ=0.4','λ=0.7','Location','northwest');
    hold on 
    box off;
    set(gca, 'xTick', [0:1:2]);
    set(gca,'XTickLabel',{'','','','','',''});
    set(gca, 'yTick', [0:0.2:1]);
    set(gca,'YTickLabel',{'','','','','',''});
    set(gca,'xminortick','on'); %
    set(gca,'yminortick','on'); %
    ax = gca;
    ax.XAxis.MinorTickValues = 0:0.1:2.6;
    ax.YAxis.MinorTickValues = 0:0.1:1; %
    set(gca,'TickLength',[0.018,0.009]); %
    set(gca,'linewidth',lw); %
    set(gca,'XAxisLocation','bottom');
    set(gca,'TickDir','out'); %
    %
    XL = get(gca,'xlim'); XR = XL(2);
    YL = get(gca,'ylim'); YT = YL(2);
    plot(XL,YT*ones(size(XL)),'color','k','linewidth',lw);
    plot(XR*ones(size(YL)),YL,'color','k','linewidth',lw);

    set(gca, 'LooseInset', [0.02,0.02,0.003,0.02]); %
    set(gcf,'unit','centimeters','position',[7 20 3.19 3.49]);%
end

%% 3-2-1
color1 = [0.00,0.45,0.74];
color4 = [0.85,0.33,0.10];
color7 = [0.93,0.69,0.13];
figure;
bar(rsAvFr,'FaceColor',color7);
hold on;
ylim([10,11]);
Leg = legend('   ','Location','north','Orientation','vertical');
legend('boxoff');
Leg.ItemTokenSize = [20,2.5];%
ax = get(gca);
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off'; %
set(gca, 'LooseInset', [0.02,0.02,0.003,0.02]); %
set(gcf,'unit','centimeters','position',[7 20 2 2]);%


%% 4 
beta = 1;
s = 1;
Mu = [1,0.5,0.1,0.05,0.01,0.005,0.001,0.0005,0.0001];
for m = 2:9
    tic
    mu = Mu(m);
    for strNum = 1:8
        LCRname = strcat('Indirect-reciprocity-group-reputation\results\formal\lambda evolve\coop_rate\147\initrand\L',num2str(strNum),'cooprate.mat');
        load(LCRname);
        allCR(isnan(allCR)) = 0;
        allcoop = []; allAvFr = [];
        for alpha = 0:0.1:2.6
            SMPname = strcat('Indirect-reciprocity-group-reputation\results\formal\lambda evolve\μ→0&SMP\β=',num2str(beta),'\α=',num2str(alpha),'\SMP\L',num2str(strNum),'_SMP.mat');
            load(SMPname);
            [coop,AvFr] = laSMEq(mu,s,Popstc30,SMP,allCR);
            allcoop = [allcoop,coop];
            allAvFr = [allAvFr;AvFr'];
        end
        filename = strcat('Indirect-reciprocity-group-reputation\results\formal\lambda evolve\μgg0\μ=',num2str(mu),'\β=',num2str(beta),'\L',num2str(strNum),'coop&AvFr.mat');
        save(filename,'allAvFr','allcoop');
    end
    toc
end

%% 5-1 
for m = 1:9
    mu = Mu(m);
    markers = ['v','^','d','<','x','s','o','>'];
    figure;
    Alpha = 0:0.1:2.6;
    for strNum = 1:8
        dataname = strcat('Indirect-reciprocity-group-reputation\results\formal\lambda evolve\μgg0\μ=',num2str(mu),'\β=',num2str(beta),'\L',num2str(strNum),'coop&AvFr.mat');
        load(dataname);
        color = L_color(strNum,:);
        marker = markers(strNum);
        plot(Alpha,allcoop,'Marker',marker,'color',color,'LineWidth',3);
        hold on
    end
    xlabel('α');
    ylabel('Cooperation Rate');
    xlim([-0.05,2.65]);
    ylim([-0.05,1.05]);
    %title('μ=0.01');
    legend('L1','L2','L3','L4','L5','L6','L7','L8','Location','best');
    figpath = strcat('Indirect-reciprocity-group-reputation\results\formal\lambda evolve\μgg0\μ=',num2str(mu),'\β=',num2str(beta),'\coop.fig');
    savefig(figpath);
end

%% 5-2 μ>>0
beta=1;
Alpha = 0:0.1:2.6;
for m = 1:9
    mu = Mu(m);
    figure;
    for strNum = 1:8
        dataname = strcat('Indirect-reciprocity-group-reputation\results\formal\lambda evolve\μgg0\μ=',num2str(mu),'\β=',num2str(beta),'\L',num2str(strNum),'coop&AvFr.mat');
        load(dataname);
        sub = subplot(2,4,strNum);
        %X = rsAvFr;
        bar(Alpha,allAvFr,'stacked');
        xlabel('α');
        ylabel('Parameter Abundance');
        %xticks([0,5,10,15,20,25]);xticklabels({0,0.5,1,1.5,2,2.5});
        %xlim([0,2.7]);
        ylim([0,1]);
        subtitle(strcat('L',num2str(strNum)));
        legend('λ=0.1','λ=0.4','λ=0.7','Location','northwest');
        hold on   
    end
    hold off
    figpath = strcat('Indirect-reciprocity-group-reputation\results\formal\lambda evolve\μgg0\μ=',num2str(mu),'\β=',num2str(beta),'\AvFr.fig');
    savefig(figpath);
end

%% 6-1 resort fixation probability
beta = 1;
L1 = [];L2 = [];L3 = [];L4 = [];L5 = [];L6 = [];L7 = [];L8 = [];
for alpha = 0:0.1:2.6
    dataname = strcat('Indirect-reciprocity-group-reputation\results\formal\lambda evolve\μ→0&SMP\β=',num2str(beta),'\α=',num2str(alpha),'\fixprob.mat');
    load(dataname);
    L1 = [L1,FXP(1,:)'];L2 = [L2,FXP(2,:)'];L3 = [L3,FXP(3,:)'];L4 = [L4,FXP(4,:)'];L5 = [L5,FXP(5,:)'];L6 = [L6,FXP(6,:)'];L7 = [L7,FXP(7,:)'];L8 = [L8,FXP(8,:)'];
    filename = strcat('Indirect-reciprocity-group-reputation\results\formal\lambda evolve\μ→0&SMP\β=',num2str(beta),'\resort\fixprob.mat');
    save(filename,'L1','L2','L3','L4','L5','L6','L7','L8');
end
%% 6-2 resort fixation probability visualization
lw1=0.5;lw2=1.2;ms=4;
beta=1;
dataname = strcat('Indirect-reciprocity-group-reputation\results\formal\lambda evolve\μ→0&SMP\β=',num2str(beta),'\resort\fixprob.mat');
load(dataname);
Alpha = 0:0.2:2.6;
figure;
for strNum = 4
    %sub = subplot(2,4,strNum);
    data = eval(strcat('L',num2str(strNum)));
    hold on
    axis([-0.1,2.7,-0.05,0.7]);
    %plot([-0.1,2.7],[1/30,1/30],'Color',[0.75,0.75,0.75],'LineStyle','--','LineWidth',lw1);
    plot(Alpha,data(1,1:2:27),'Color',[0.00,0.45,0.74],'LineStyle',':','LineWidth',lw2,'Marker','^','MarkerFaceColor',[0.00,0.45,0.74],'MarkerSize',ms);
    plot(Alpha,data(2,1:2:27),'Color',[0.85,0.33,0.10],'LineStyle',':','LineWidth',lw2,'Marker','o','MarkerFaceColor',[0.85,0.33,0.10],'MarkerSize',ms);
    plot(Alpha,data(3,1:2:27),'Color',[0.00,0.45,0.74],'LineStyle',':','LineWidth',lw2,'Marker','square','MarkerFaceColor',[0.00,0.45,0.74],'MarkerSize',ms);
    plot(Alpha,data(4,1:2:27),'Color',[0.93,0.69,0.13],'LineStyle',':','LineWidth',lw2,'Marker','o','MarkerFaceColor',[0.93,0.69,0.13],'MarkerSize',ms);
    plot(Alpha,data(5,1:2:27),'Color',[0.85,0.33,0.10],'LineStyle',':','LineWidth',lw2,'Marker','square','MarkerFaceColor',[0.85,0.33,0.10],'MarkerSize',ms);
    plot(Alpha,data(6,1:2:27),'Color',[0.93,0.69,0.13],'LineStyle',':','LineWidth',lw2,'Marker','^','MarkerFaceColor',[0.93,0.69,0.13],'MarkerSize',ms);
    plot([-0.1,2.7],[1/30,1/30],'Color',[0.75,0.75,0.75],'LineStyle','--','LineWidth',lw1);
    %legend(strcat('λ=0.4→λ=0.1'),strcat('λ=0.1→λ=0.4'),strcat('λ=0.7→λ=0.1'),strcat('λ=0.1→λ=0.7'),strcat('λ=0.7→λ=0.4'),strcat('λ=0.4→λ=0.7'),'','Location','north');
    %xlabel('α'); ylabel('fixation probability');
    %subtitle(strcat('L',num2str(strNum)));
    box off;
    set(gca, 'xTick', [0:1:2]);
    set(gca,'XTickLabel',{'','',''});
    set(gca, 'yTick', [0:0.1:7]);
    set(gca,'YTickLabel',{'','','','','','','',''});
    set(gca,'xminortick','on'); %
    set(gca,'yminortick','on'); %
    ax = gca;
    ax.XAxis.MinorTickValues = 0:0.2:2.6;
    ax.YAxis.MinorTickValues = 0:0.1:0.7; %
    set(gca,'TickLength',[0.018,0.009]); %
    set(gca,'linewidth',lw1); %
    set(gca,'XAxisLocation','bottom');
    set(gca,'TickDir','out'); %
    %
%     XL = get(gca,'xlim'); XR = XL(2);
%     YL = get(gca,'ylim'); YT = YL(2);
%     plot(XL,YT*ones(size(XL)),'color','k','linewidth',lw);
%     plot(XR*ones(size(YL)),YL,'color','k','linewidth',lw);
    set(gca, 'LooseInset', [0.02,0.02,0.003,0.02]); %
    set(gcf,'unit','centimeters','position',[7 20 3.19 3.49]);%

end

%% 6-2-1 
figure;
% plot(Alpha,data(1,1:2:27),'Color',[0.00,0.45,0.74],'LineStyle',':','LineWidth',lw2,'Marker','^','MarkerFaceColor',[0.00,0.45,0.74],'MarkerSize',ms);
% plot(Alpha,data(2,1:2:27),'Color',[0.85,0.33,0.10],'LineStyle',':','LineWidth',lw2,'Marker','o','MarkerFaceColor',[0.85,0.33,0.10],'MarkerSize',ms);
% plot(Alpha,data(3,1:2:27),'Color',[0.00,0.45,0.74],'LineStyle',':','LineWidth',lw2,'Marker','square','MarkerFaceColor',[0.00,0.45,0.74],'MarkerSize',ms);
% plot(Alpha,data(4,1:2:27),'Color',[0.93,0.69,0.13],'LineStyle',':','LineWidth',lw2,'Marker','o','MarkerFaceColor',[0.93,0.69,0.13],'MarkerSize',ms);
% plot(Alpha,data(5,1:2:27),'Color',[0.85,0.33,0.10],'LineStyle',':','LineWidth',lw2,'Marker','square','MarkerFaceColor',[0.85,0.33,0.10],'MarkerSize',ms);
plot(Alpha,data(6,1:2:27),'Color',[0.93,0.69,0.13],'LineStyle',':','LineWidth',lw2,'Marker','^','MarkerFaceColor',[0.93,0.69,0.13],'MarkerSize',ms);
hold on;
ylim([10,11]);
Leg = legend('   ','Location','north','Orientation','vertical');
legend('boxoff');
Leg.ItemTokenSize = [18,2.5];%
ax = get(gca);
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off'; %
set(gca, 'LooseInset', [0.02,0.02,0.003,0.02]); %
set(gcf,'unit','centimeters','position',[7 20 0.7 1.3]);%






