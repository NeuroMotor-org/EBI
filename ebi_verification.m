function []=ebi_verification(varargin)

% ebi_verification(test_n), n: vector with the tests/figures to generate, e.g.: % 1,2,3, [1 2], [1 3 4]
% Syntax: ebi_verification() : does all tests and generated all figures.

% 1: Sample Report Figure (fig 1)
% 2: Comparison againts Truth for 1 (fig 2)
% 3: General Performance Test Againts Truth (data file, no figures)
% 4: Plot the selected data for 3 (fig. 3)
% 5. Uni-Variate Test example

if isempty(varargin)
    test_n=1:4;
else
    test_n=varargin{1};
end

% basics
FixedRandomSeed=1849872129;
rng(FixedRandomSeed); % initialise randome number egneration with fixed seed (for reproduciability)

%% I-II. Run a multivariate example and compare it againts Real Truth
if any(ismember([1 2],test_n))
    % define variables and groups
    Gold=[zeros(1,1600),ones(1,400)]; % define reality as 1000 varibales with real efects out of 2000.
    NGold=[sum(Gold==0) sum(Gold==1)*3/4 sum(Gold==1)*1/4]; % Number of variables
    Ndim=sum(NGold); % number of varibales
    mu=[0 -1 1.5]; % mean difference of varibales from 0: null group are really zero, non-null group are either centred at -1 (750 varibales) or  4.7 (250 variables)
    sigma=[1 1 1]; % standard deviation of the null group and and non-null subgroup data.
    % Observations/Subjects
    N=[20 60]; % number of observations in the first and second sample group. For example, can be considered as 20 control (healthy) subjects and 80 patients.
    Nsubject=sum(N); % total subjects (observations): 100
    % create data chunks and merge them
    X1=mu(1)+sigma(1).*randn(sum(N),NGold(1));
    X2=[mvnrnd(mu(1).*ones(N(1),sum(NGold(2:3))),sigma(1).*eye(sum(NGold(2:3))));[mvnrnd(mu(2)*ones(N(2),NGold(2)),sigma(2).*eye(NGold(2))),mvnrnd(mu(3)*ones(N(2),NGold(3)),sigma(3).*eye(NGold(3)))]]; % build non-null data
    %Noise=mvnrnd(zeros(Nsubject,Ndim),1*eye(Ndim)); % add some noise to all data
    % Define data and group labels.
    X=[X1,X2];%+Noise; % define the whole dataset
    g=[zeros(N(1),1);ones(N(2),1)];% define group labels
    rng('shuffle'); 
    % Analyse
    tic;O=ebi(X,g);toc % Run ebi on data
    [~,~,~,ttest2stat]=ttest2(X(g==0,:),X(g>0,:),'Vartype','unequal'); % run ttest on data to prepare teh results for the R implementation
    zz=icdf('Normal',cdf('T',ttest2stat.tstat,ttest2stat.df),0,1); % find z-scores from t's
    [ locfdrOutputs, locfdrcmdOut ] = locfdrWrap( zz, 120, 7, 0, [0.4 0.6], 3,  0, 0, [0.5 1 1.5],[0 1],'label',0); %#ok<ASGLU> % Run the R impelmentation of
end

if ismember(1,test_n)
    ebi_report(O,'MS_fig1_example') % generate report, this will be MS Fig. 1
end

if ismember(2,test_n)
    % calculate the real FDR and Beta for the EBI and R's implementation, with respect to the Real Truth
    RealFDR=O.Prlist.*0;
    RealAlpha=O.Prlist.*0;
    RealBeta=O.Prlist.*0;
    locfdrReAlpha=O.Prlist.*0;
    locfdrReFDR=O.Prlist.*0;
    for Pri=1:length(O.Prlist)
        DecisionPr=O.Prlist(Pri);
        DetectionSamples=O.Posterior(2,:)>=DecisionPr;
        ConfusionMatrixEBI=[sum(Gold.*DetectionSamples),sum(~Gold.*DetectionSamples);sum(Gold.*~DetectionSamples),sum(~Gold.*~DetectionSamples)]; % [TP, FP; FN, TN]
        RealFDR(Pri)=ConfusionMatrixEBI(1,2)./(ConfusionMatrixEBI(1,1)+ConfusionMatrixEBI(1,2)); % FP/(FP+TP)
        RealAlpha(Pri)=ConfusionMatrixEBI(1,2)./(ConfusionMatrixEBI(2,2)+ConfusionMatrixEBI(1,2)); % FP/(FP+TN)
        RealBeta(Pri)=ConfusionMatrixEBI(2,1)./(ConfusionMatrixEBI(1,1)+ConfusionMatrixEBI(2,1)); % FN/(FN+TP)
        DetectionLOCFDR=(1-locfdrOutputs.mat(:,2))>DecisionPr;
        locfdrReAlpha(Pri)=sum(locfdrOutputs.mat(DetectionLOCFDR,6))./sum(locfdrOutputs.mat(:,6)); % checked: OK. Classic Frequentist alpha at each decision criterion
        locfdrReFDR(Pri)=mean(locfdrOutputs.mat(DetectionLOCFDR,2).*locfdrOutputs.mat(DetectionLOCFDR,5))./mean(locfdrOutputs.mat(DetectionLOCFDR,5));
    end
    
    % Generate comparison figure: MS Fig. 2.
    MyLineWidth=1;
    MyFontSize=8;
    MyTitleFontSize=10;
    RealColor=[0.75 0 0];
    RColor=[0.5 0.25 0.125];
    close all
    fh1=figure('Visible','off');
    set(gcf, 'PaperPositionMode', 'manual','PaperUnits', 'centimeters','PaperPosition', [0 0 16 6]);
    
    % FDR
    subplot(1,2,1);
    plot(O.Prlist,RealFDR,'Color',RealColor,'LineStyle','-','LineWidth',MyLineWidth);
    hold on;
    plot(O.Prlist,O.EstimatedFDR,'Color',[0 0 0],'LineStyle','-','LineWidth',MyLineWidth);
    plot(O.Prlist,locfdrReFDR,'Color',RColor,'LineStyle',':','LineWidth',MyLineWidth);
    xlim([0 1])
    ylim([0 1])
    legend({'Real','Current EBI','EBI (R)'},'Location','best')
    set(gca,'FontSize', MyFontSize);
    xlabel('P_{1,crit}','FontSize', MyFontSize,'FontWeight','bold');
    title('FDR','FontSize', MyTitleFontSize,'FontWeight','bold');
    
    subplot(1,2,2);
    plot(O.Prlist,RealBeta,'Color',RealColor,'LineStyle','-','LineWidth',MyLineWidth);
    hold on;
    plot(O.Prlist,O.EstimatedBeta,'Color',[0 0 0],'LineStyle','-','LineWidth',MyLineWidth);
    plot(1-locfdrOutputs.cdf1(:,1),1-locfdrOutputs.cdf1(:,2),'Color',RColor,'LineStyle',':','LineWidth',MyLineWidth);
    xlim([0 1])
    ylim([0 1])
    legend({'Real','Current EBI','EBI (R)'},'Location','best')
    set(gca,'FontSize', MyFontSize);
    xlabel('P_{1,crit}','FontSize', MyFontSize,'FontWeight','bold');
    title('{\boldmath$\beta_m$}','FontSize', MyTitleFontSize,'FontWeight','bold','Interpreter','latex');
    
    tifffilename='MS_fig2_compare';
    print(fh1,tifffilename,'-dpng', '-r600');
    disp([tifffilename ' was generated.'])
    close(fh1);
end
%% III. Test Performance in different Conditions.
if ismember(3,test_n)
    % Nvar (200-2000) and p0 ratio (0.25, 0.75), Ncon(25, 100), Npat (25, 100).
    % distribution types (normal, Chi), Difference Effect Size (Cohen's d=0.2,
    % d=0.9).
    NNlist=[200, 2000];
    p0list=[0.25 0.75];
    m0list=[25 100];
    m1list=[25 100];
    distIDlist=[1 2]; % Normal, Chi
    CohenDlist=[0.5 0.9];
    Rep=5;
    
    Nconditions=length(NNlist).*length(p0list).*length(m0list).*length(m1list).*length(distIDlist).*length(CohenDlist);
    ResultsBagGold=nan(Nconditions,7,Rep);
    ResultsBagEst=nan(Nconditions,7,Rep);
    
    condi=0;
    for NNi=1:length(NNlist)
        NN=NNlist(NNi);
        for p0i=1:length(p0list)
            p0=p0list(p0i);
            for m0i=1:length(m0list)
                m0=m0list(m0i);
                for m1i=1:length(m1list)
                    m1=m1list(m1i);
                    for disti=1:length(distIDlist)
                        distID=distIDlist(disti);
                        for cohi=1:length(CohenDlist)
                            CohenD=CohenDlist(cohi);
                            condi=condi+1;
                            for Repi=1:Rep
                                %%%%%
                                
                                m=m0+m1;
                                NGold=[round(NN.*p0) NN-round(NN.*p0)];
                                Gold=[zeros(1,NGold(1)),ones(1,NGold(2))];
                                g=[zeros(m0,1);ones(m1,1)];
                                
                                rng(FixedRandomSeed);
                                X=nan(m, NN);
                                for ci=1:NN
                                    switch Gold(ci)
                                        case 0 % No real difference
                                            switch disti
                                                case 1
                                                    X(:,ci)=randn(m,1);
                                                case 2
                                                    X(:,ci)=random('beta',2,10,m,1);
                                            end
                                        case 1 % Real Difference
                                            switch disti
                                                case 1
                                                    X(1:m0,ci)=randn(m0,1);
                                                    X((m0+1):end,ci)=randn(m1,1)+CohenD;
                                                case 2
                                                    X(1:m0,ci)=random('beta',2,10,m0,1);
                                                    X((m0+1):end,ci)=random('beta',2,10,m1,1)+(0.1.*CohenD);
                                            end
                                    end
                                end
                                rng('shuffle');
                                O=ebi(X,g);
                                Prlist=O.Prlist;
                                O.IX=X;
                                O.Icond=[NN,p0,m0,m1,distID,CohenD];
                                O.IcondLab={'NN','p0','m0','m1','distID','CohenD'};
                                O.IGold=Gold;
                                O.Ig=g;
                                
                                RealFDR=Prlist.*0;
                                RealAlpha=Prlist.*0;
                                RealBeta=Prlist.*0;
                                for Pri=1:length(Prlist)
                                    DecisionPr=Prlist(Pri);
                                    DetectionSamples=O.Posterior(2,:)>=DecisionPr;
                                    ConfusionMatrixEBI=[sum(Gold.*DetectionSamples),sum(~Gold.*DetectionSamples);sum(Gold.*~DetectionSamples),sum(~Gold.*~DetectionSamples)]; % [TP, FP; FN, TN]
                                    RealFDR(Pri)=ConfusionMatrixEBI(1,2)./(ConfusionMatrixEBI(1,1)+ConfusionMatrixEBI(1,2)); % FP/(FP+TP)
                                    RealAlpha(Pri)=ConfusionMatrixEBI(1,2)./(ConfusionMatrixEBI(2,2)+ConfusionMatrixEBI(1,2)); % FP/(FP+TN) % Perfect!
                                    RealBeta(Pri)=ConfusionMatrixEBI(2,1)./(ConfusionMatrixEBI(1,1)+ConfusionMatrixEBI(2,1)); % FN/(FN+TP) % Perfect
                                end
                                
                                O.RealFDR=RealFDR;
                                O.RealAlpha=RealAlpha;
                                O.RealBeta=RealBeta;
                                
                                ResultsGold=nan(1,7);
                                ResultsEst=nan(1,7);
                                % @ p1=p1 (The estimated ratio of variables with effect)
                                ResultsGold(1,1)=NGold(2)./sum(NGold);
                                ResultsEst(1,1)=O.EstimatedPriors(2);
                                % @ p1, FDR (The FDR at the estimated p1)
                                ResultsGold(1,2)=RealFDR(round(O.ThresholdMatrix(4,1).*1000+1));
                                ResultsEst(1,2)=O.ThresholdMatrix(4,2);
                                % @ p1, Beta (The Beta at the estimated p1)
                                ResultsGold(1,3)=RealBeta(round(O.ThresholdMatrix(4,1).*1000+1));
                                ResultsEst(1,3)=O.ThresholdMatrix(4,3);
                                
                                % @ FDR=0.05, FDR (The estimatied FDR vs. 0.05)
                                ResultsGold(1,4)=RealFDR(round(O.ThresholdMatrix(2,1).*1000+1));
                                ResultsEst(1,4)=O.ThresholdMatrix(2,2);
                                
                                % @ FDR=0.05, Beta (The estimatied Beta at FDR=0.05)
                                ResultsGold(1,5)=RealBeta(round(O.ThresholdMatrix(2,1).*1000+1));
                                ResultsEst(1,5)=O.ThresholdMatrix(2,3);
                                
                                % @ Beta=0.2, Beta (The estimatied Beta at Beta=0.2)
                                ResultsGold(1,6)=RealBeta(round(O.ThresholdMatrix(3,1).*1000+1));
                                ResultsEst(1,6)=O.ThresholdMatrix(3,3);
                                
                                % @ Beta=0.2, FDR (The estimatied FDR at Beta=0.2)
                                ResultsGold(1,7)=RealBeta(round(O.ThresholdMatrix(3,1).*1000+1));
                                ResultsEst(1,7)=O.ThresholdMatrix(3,2);
                                
                                O.ResultsGold=ResultsGold;
                                O.ResultsEst=ResultsEst;
                                save(['TestResults_' num2str(condi) '_Data_NN' num2str(NN) 'p0' num2str(p0) 'm0' num2str(m0) 'm1' num2str(m1) 'distID' num2str(distID) 'CohenD' num2str(CohenD) 'Rep' num2str(Repi) '.mat'],'O');
                                ResultsBagGold(condi,:,Repi)=ResultsGold;
                                ResultsBagEst(condi,:,Repi)=ResultsEst;
                                %%%%%
                            end
                        end
                    end
                end
            end
        end
    end
    
    save('TestResults_All.mat','ResultsBagGold','ResultsBagEst');
    % TP FP
    % FN TN
end

if ismember(4,test_n)
    NNlist=[200, 2000];
    p0list=[0.25 0.75];
    m0list=[25 100];
    m1list=[25 100];
    CohenDlist=[0.5 0.9];
    
    MyColorList=[0 0 0;0.6 0.6 0.6;0 0 0.85;0 0.85 0];
    load('TestResults_All.mat','ResultsBagGold','ResultsBagEst');
    
    FigLabelList={'Estimated Prior p_1', 'FDR @ 0.05','Beta @ 0.2'};
    MyLegendLabels={['m_0=' num2str(m0list(1)) ', m_1=' num2str(m1list(1))],['m_0=' num2str(m0list(1)) ', m_1=' num2str(m1list(2))],['m_0=' num2str(m0list(2)) ', m_1=' num2str(m1list(1))],['m_0=' num2str(m0list(2)) ', m_1=' num2str(m1list(2))]};
    
    MyFontSize=6;
    MyTitleFontSize=10;
    close all
    
    plotX=16; plotY=22;
    DX=ones(1,4);
    dx=[0.1 0.1.*ones(1,3) 0.05];
    TempSum=sum([DX dx]);
    DX=DX./TempSum; % normalise
    dx=dx./TempSum; % normalsie
    DY=DX(1).*0.62*ones(1,6).*plotX./plotY;
    dy=[0.06 0.04 0.07 0.04 0.07 0.04 0.065];
    dy(end)=dy(end)+(1-sum([DY dy])); % trim to 1.
    ah=zeros(length(DX),length(DY));
    mysubplotpos=bsubplotpos(DX,DY,dx,dy,1);
    
    fh1=figure('Visible','off');
    set(gcf, 'PaperPositionMode', 'manual','PaperUnits', 'centimeters','PaperPosition', [0 0 plotX plotY],'PaperSize', [plotX plotY]);
    figlist=[1 4 6];
    for figii=1:3
        figi=figlist(figii);
        % Visual Inspection Figures
        for ploti=1:8
            ah(ploti+(figii-1)*8)=axes('position',mysubplotpos{ploti+(figii-1)*8},'XGrid','off','XMinorGrid','off');
            cID0=arrayfun(@(x) str2double(x),dec2bin(ploti-1));
            cID=zeros(1,3);
            cID((end-length(cID0)+1):end)=cID0;
            PlotDataG1=squeeze(ResultsBagGold(32*cID(1)+16*cID(2)+cID(3)+(1:4:16),figi,:));
            PlotDataE1=squeeze(ResultsBagEst(32*cID(1)+16*cID(2)+cID(3)+(1:4:16),figi,:));
            PlotDataG2=squeeze(ResultsBagGold(32*cID(1)+16*cID(2)+cID(3)+(3:4:16),figi,:));
            PlotDataE2=squeeze(ResultsBagEst(32*cID(1)+16*cID(2)+cID(3)+(3:4:16),figi,:));
            
            switch figii
                case 1
                    plot([0.65 1.35],[1 1].*PlotDataG1(1),'LineWidth',0.5,'Color',[0.75 0 0]);
                    hold on
                    plot([1.65 2.35],[1 1].*PlotDataG1(1),'LineWidth',0.5,'Color',[0.75 0 0]);
                    
                    for ci=1:4
                        plot(0.75+ci/10,PlotDataE1(ci,:),'LineStyle','none','Marker','.','Color',MyColorList(ci,:));
                        plot(1.75+ci/10,PlotDataE2(ci,:),'LineStyle','none','Marker','.','Color',MyColorList(ci,:)); %,'DisplayName',MyLegendLabels{ci}
                        
                    end
                case {2,3}
                    plot([0.65 1.35],[1 1].*PlotDataE1(1),'LineWidth',0.5,'Color',[0.75 0 0]);
                    hold on
                    plot([1.65 2.35],[1 1].*PlotDataE1(1),'LineWidth',0.5,'Color',[0.75 0 0]);
                    
                    for ci=1:4
                        plot(0.75+ci/10,PlotDataG1(ci,:),'LineStyle','none','Marker','.','Color',MyColorList(ci,:));
                        plot(1.75+ci/10,PlotDataG2(ci,:),'LineStyle','none','Marker','.','Color',MyColorList(ci,:));
                    end
            end

            xlim([0.5 2.5]);
            ylim([0 1])
            title(['N=' num2str(NNlist(1+cID(1))) ', p_1=' num2str(p0list(1+cID(2))) ', d_c=' num2str(CohenDlist(1+cID(3)))],'FontSize', MyFontSize,'FontWeight','bold');
            if ploti==5 && figii==3
                xlabel('Distribution','FontSize', MyFontSize,'FontWeight','bold');
            end
            set(ah(ploti+(figii-1)*8),'position',mysubplotpos{ploti+(figii-1)*8},'XGrid','off','XMinorGrid','off','FontSize', MyFontSize,'XTick',[1 2],'XTickLabel',{'Normal','Beta'});
        end
        
    end % figi
    
    titleah = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    for hi=1:4
        plot(-1,-1,'LineStyle','none','Marker','.','MarkerSize',20,'Color',MyColorList(hi,:));
        hold on
    end
    legend(MyLegendLabels,'Location','southeast','Orientation','horizontal');
    set(titleah, 'YTick', [], 'XTick', [], 'YTickLabel', [], 'XTickLabel', [],'Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    for figii=1:3
        text(0.5,1-0.32*(figii-1),0,FigLabelList{figii},'HorizontalAlignment','center','VerticalAlignment','top','FontWeight','bold','FontSize',MyTitleFontSize+2, 'Color', [0 0 0]);
    end
    
    tifffilename='MS_fig3_PerformanceTest';
    print(fh1,tifffilename,'-dpng', '-r600');
    disp([tifffilename ' was generated.'])
    close(fh1);
end
%% V. UniVariate Case
if ismember(5,test_n)
    rng(FixedRandomSeed);
    p0=0;
    N1=15;
    N2=25;
    N2a=round(p0*N2);
    N2b=N2-N2a;
    X1=0+1.*randn(N1,1);
    X2=[0+1.*randn(N2a,1);-0.5+1.*randn(N2b,1)];
    X=[X1;X2];
    g=[zeros(N1,1);ones(N2,1)];
    rng('shuffle');
    O=ebi_uv(X,g);
    
    %%
    MyFontSize=6;
    MyTitleFontSize=10;
    plotX=8;
    plotY=6;
    fh1=figure('Visible','off');
    set(gcf, 'PaperPositionMode', 'manual','PaperUnits', 'centimeters','PaperPosition', [0 0 plotX plotY],'PaperSize', [plotX plotY]);
    plot(O.EstimatedBeta,O.EstimatedAlpha,'Color','k','LineWidth',1);
    set(gca,'FontSize',MyFontSize,'FontWeight','normal')
    xlabel('\beta','FontSize',MyTitleFontSize,'FontWeight','bold')
    ylabel('\alpha','FontSize',MyTitleFontSize,'FontWeight','bold')
    tifffilename='MS_fig4_univariate';
    print(fh1,tifffilename,'-dpng', '-r600');
    disp([tifffilename ' was generated.'])
    close(fh1);
end
end


%% Functions
function [ Outputs, cmdOut ] = locfdrWrap(zz, bre, df, pct, pct0, nulltype, type, plot, mult, mlests, main,sw) %#ok<*INUSD>
% Calls R and runs locfdr on the data and brings back the results to MATLAB
%   Detailed explanation goes here

% debug
% [ Outputs, cmdOut ] = locfdrWrap( [randn(1,180) 1+randn(1,10) -0.7+randn(1,10)], 10, 10, 0, [0.47 53], 1,  0, 0, [0.5 1 1.5],[0 1], 'label',2)

save('locfdrdatain.mat','zz', 'bre', 'df', 'pct', 'pct0', 'nulltype', 'type', 'plot', 'mult', 'mlests','main','sw','-v4');
[cmdStatus,cmdOut]=system('"C:\Program Files\R\R-3.4.4\bin\Rscript.exe" locfdrscript.r');
if cmdStatus==0
    Outputs=load('locfdrdataout.mat');
else
    disp('Error !!!')
end
if ~exist('Outputs','var')
    Outputs=[];
    % mat:    x, fdr, Fdrleft, Fdrright, f, f0, f0theo, fdrtheo, counts, lfdrse, p1f1
end
end

%%
% Copyright (c) 2018 Bahman Nasseroleslami, All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 
%     Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
%     Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
%     Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission. 
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.