function []=ebi_verification(varargin)
% Empirical Bayesian Inference for Multiple Comparisons
%
% Usage:
% ebi_verification(test_n)
% ebi_verification()
%
% n: vector with the tests/figures to generate, e.g.: % 1,2,3, [1 2], [1 3 4], % no argument ebi_verification() does all tests and generated all figures.
% Output: Data/figure files
%
% n:
% 1. Sample Report Figure (fig 1)
% 2. Comparison againts Truth for 1 (fig 2)
% 3. General Performance Test Againts Truth (data file, no figures)
% 4. Plot the selected data for 3 (fig. 3)
% 5. Uni-Variate Test example (fig. 4)
% 6. Experimnetal data example (fig. 5)

% Example:
% ebi_verification(1);
% ebi_verification(test_n)([1 3 4]);
%
% Written by:
% Bahman Nasseroleslami, Trinity College Dublin, the University of Dublin, nasserob@tcd.ie, bahman@neuromotor.org
% Part of the Emprical Bayesian Inference (EBI) toolbox for MATLAB
% Revision: 26/7/2019. 

if isempty(varargin)
    test_n=1:6;
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
    mu=[0 -1 1.5]; % mean difference of varibales from 0: null group are really zero, non-null group are either centred at -1 (750 varibales) or  4.7 (250 variables)
    sigma=[1 1 1]; % standard deviation of the null group and and non-null subgroup data.
    % Observations/Subjects
    N=[20 60]; % number of observations in the first and second sample group. For example, can be considered as 20 control (healthy) subjects and 80 patients.
    % create data chunks and merge them
    X1=mu(1)+sigma(1).*randn(sum(N),NGold(1));
    X2=[mvnrnd(mu(1).*ones(N(1),sum(NGold(2:3))),sigma(1).*eye(sum(NGold(2:3))));[mvnrnd(mu(2)*ones(N(2),NGold(2)),sigma(2).*eye(NGold(2))),mvnrnd(mu(3)*ones(N(2),NGold(3)),sigma(3).*eye(NGold(3)))]]; % build non-null data
    % Define data and group labels.
    X=[X1,X2]; % define the whole dataset
    g=[zeros(N(1),1);ones(N(2),1)];% define group labels
    rng('shuffle');
    % Analyse
    tic;O=ebi(X,g);toc % Run ebi on data
    [~,~,~,ttest2stat]=ttest2(X(g==0,:),X(g>0,:),'Vartype','unequal'); % run ttest on data to prepare the results for the R implementation
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
    dx=[0.13 0.1.*ones(1,3) 0.02];
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
            set(ah(ploti+(figii-1)*8),'position',mysubplotpos{ploti+(figii-1)*8},'XGrid','off','XMinorGrid','off','FontSize', MyFontSize,'XTick',[1 2],'YTick',0:0.25:1,'XTickLabel',{'Normal','Beta'});
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
    X2=[0+1.*randn(N2a,1);-1+1.*randn(N2b,1)];
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
%% VI. Real EEG Data
if ismember(6,test_n)
    %% Load/Analyse
    rng(FixedRandomSeed);
    SR=2048;
    FreqBands={2:4,5:7,8:10,11:13,14:20,21:30,31:47,53:97,103:147,153:197,203:247,253:297,303:347,353:397}; % Revised
    NfreqBands=length(FreqBands); % number of frequency bands
    EpochTime=1; % in seconds
    EpochLength=SR*EpochTime;
    fiMask=zeros(NfreqBands,EpochLength); % The matrix mask to define the frequency band for bcoherence function
    for fi=1:NfreqBands
        fiMask(fi,FreqBands{fi}+1)=1;
    end 
    
    rawdatafilename1='MS_EBI_F.mat';
    rawdatafilename0='MS_EBI_RS.mat';
    datafilename='MS_EBI_ExpSpectral.mat';
    ThisFreqCompare=1:7;
    if exist(datafilename,'file')==2
        load(datafilename,'O','Pvalues','X0','X1'); %#ok<NASGU>
        disp(['Loaded MAT results file: ' datafilename]);
    else
        % X1: Force Data
        load(rawdatafilename1,'SubjectData');
        disp(['Loaded MAT raw data file: ' rawdatafilename1]);
        GoodIDs=all(SubjectData.Exp2s2EXGmask(:,:)==1,1);
        [~,Ntime,~]=size(SubjectData.Exp2s2EXG);
        NGepochs=sum(GoodIDs);
        
        %%%
        X1=zeros(NGepochs,128,7);
        g1=ones(NGepochs,1);
        %Pvalues=nan(1,128);
        for ch1=1:128
            X1(:,ch1,:)=fft(bsxfun(@times,detrend(squeeze(SubjectData.Exp2s2EXG(ch1,:,GoodIDs))).',hann(Ntime).'./sqrt(sum((hann(Ntime).^2)))).').'*(fiMask(1:7,:).');
        end
        clearvars SubjectData
        
        % X0: Resting-State Data
        load(rawdatafilename0,'SubjectData');
        disp(['Loaded MAT raw data file: ' rawdatafilename0]);
        GoodIDs0=all(SubjectData.Exp0EXGmask(:,:)==1,1);
        [~,Ntime0,~]=size(SubjectData.Exp0EXG);
        NGepochs0=sum(GoodIDs0);
        
        %%%
        X0=zeros(NGepochs0,128,7);
        g0=zeros(NGepochs0,1);
        Pvalues=nan(length(ThisFreqCompare),128);
        for ch1=1:128
            X0(:,ch1,:)=fft(bsxfun(@times,detrend(squeeze(SubjectData.Exp0EXG(ch1,:,GoodIDs0))).',hann(Ntime0).'./sqrt(sum((hann(Ntime0).^2)))).').'*(fiMask(1:7,:).');
            for fi1=ThisFreqCompare
                Pvalues(fi1,ch1)=ranksum(abs(X0(:,ch1,fi1)),abs(X1(:,ch1,fi1)));
            end
        end
              
        % Test Spectral Power:
        O=ebi([reshape(abs(squeeze(X0(:,:,ThisFreqCompare))),NGepochs0,[]);reshape(abs(squeeze(X1(:,:,ThisFreqCompare))),NGepochs,[])],[g0;g1],'Test','AUCZ','Thresholds',[0.9 0.05 0.1 0.05]);
        ebi_report(O,'MS_fig_5_expSpectralPower_report')
        
        save(datafilename,'O','Pvalues','X0','X1');
        disp(['Saved MAT data and Results file: ' datafilename]);
    end
    %% Plot
    close all
    warning('off','MATLAB:handle_graphics:exceptions:SceneNode')
    bcolormapart4=createColormap(0);
    DiffColorMap=bcolormapart4;%[linspace(0.01,0.99,32).' linspace(0.25,0.99,32).' linspace(0.5,0.99,32).';linspace(0.99,0.5,32).' linspace(0.99,0.01,32).' linspace(0.99,0.01,32).'];
    fiLabelLIST={'\delta_{ } (2-4Hz)','\theta_{ } (5-7Hz)','\alpha_l (8-10Hz)','\alpha_h (11-13Hz)','\beta_l (14-20Hz)','\beta_h (21-30Hz)','\gamma_l (31-47Hz)','\gamma_h (53-97Hz)'};
    ChX=[0 20 39 57 66 51 60 76 80 81 79 73 86 92 95 94 89 78 72 84 93 99 100 97 91 86 92 ...
        95 94 89 78 66 6 28 51 60 76 80 81 79 73 57 59 58 55 31 30 29 32 28 28 0 0 0 0 0 0 0 -31 ...
        -30 -29 -32 -28 -28 -16 -28 -51 -60 -55 -58 -59 -81 -80 -76 -49 -66 -78 -89 -94 -95 -100 ...
        -99 -93 -84 -72 -57 -39 -49 -66 -78 -89 -94 -95 -81 -80 -76 -16 -28 -51 -60 -55 -58 -59 ...
        -31 -30 -29 -32 -28 -28 0 6 28 28 0 0 0 0 0 0 31 30 29 32 28 55 58 59 57 ];
    ChY=[0 0 0 0 -28 -51 -60 -55 -58 -59 -57 -53 -28 -30 -31 -30 -29 -32 0 0 0 0 0 0 0 28 30 31 ...
        30 29 32 28 19 28 51 60 55 58 59 57 53 79 81 80 76 95 94 89 78 66 49 39 57 72 84 93 99 100 ...
        95 94 89 78 66 49 12 28 51 60 76 80 81 59 58 55 28 28 32 29 30 31 0 0 0 0 0 0 0 -28 -28 -32 ...
        -29 -30 -31 -59 -58 -55 -12 -28 -51 -60 -76 -80 -81 -95 -94 -89 -78 -66 -49 -39 -19 -28 -49 ...
        -57 -72 -84 -93 -99 -100 -95 -94 -89 -78 -66 -76 -80 -81 -79];
    ChZ=[100 98 92 82 69 69 54 36 17 -3 -23 -42 -42 -23 -3 17 36 54 69 54 36 17 -3 -23 -42 -42 -23 ...
        -3 17 36 54 69 98 92 69 54 36 17 -3 -23 -42 -23 -3 17 36 -3 17 36 54 69 82 92 82 69 54 36 17 ...
        -3 -3 17 36 54 69 82 98 92 69 54 36 17 -3 -3 17 36 82 69 54 36 17 -3 -3 17 36 54 69 82 92 82 ...
        69 54 36 17 -3 -3 17 36 98 92 69 54 36 17 -3 -3 17 36 54 69 82 92 98 92 82 82 69 54 36 17 -3 -3 ...
        17 36 54 69 36 17 -3 -23 ];

    % Find 2D-ready coordinates (for top view plots)
    alpha=ChX./sqrt(ChX.^2+ChY.^2);
    beta=ChY./sqrt(ChX.^2+ChY.^2);
    alpha(1)=1;
    beta(1)=1;
    ChXX=ChY+beta.*abs(100-ChZ);
    ChYY=-alpha.*abs(100-ChZ)-ChX;
    ChXX=ChXX./max(abs(ChXX)).*100;
    ChYY=ChYY./max(abs(ChYY)).*100;
    ChR=sqrt(ChXX.^2+ChYY.^2);
    ChXX(ChR>100)=ChXX(ChR>100)./ChR(ChR>100).*100;
    ChYY(ChR>100)=ChYY(ChR>100)./ChR(ChR>100).*100;
    clearvars alpha beta charii chari numi ngi
    
    MyFontSize=8;
    MyTitleFontSize=10;
    % Row 1: EBI full,
    % Row 2: Contours of EBI criteria.
    % Row 3: Contours of p<0.01, contours of aFDR, and Bonferroni
    plotX=24; plotY=11;
    DX=[ones(1,7) 0.1 0.1 1.6];
    dx=[0.2 0.01.*ones(1,6) 0.4 0.4 0.4 0.15];
    TempSum=sum([DX dx]);
    DX=DX./TempSum; % normalise
    dx=dx./TempSum; % normalsie
    DY=DX(1).*1*ones(1,3).*plotX./plotY;
    dy=[0.15 0.12 0.1 0.03];
    dy(end)=dy(end)+(1-sum([DY dy])); % trim to 1.
    ah=zeros(length(DX),length(DY));
    mysubplotpos=bsubplotpos(DX,DY,dx,dy,1);
    
    fh1=figure('Visible','off');
    set(gcf,'PaperPositionMode', 'manual','PaperUnits', 'centimeters','PaperPosition', [0 0 plotX plotY],'PaperSize', [plotX plotY]);
    
    MaxFreqCond=length(ThisFreqCompare);
    MyVList1=reshape(O.AUC,[128,MaxFreqCond])-0.5;
    MyVList2=Pvalues.'; % reshape(Pvalues,[128,MaxFreqCond]);
    myClim1=[-1 1]*max(abs(MyVList1(:)));
    myClim2=[0 ceil(max(-log10(MyVList2(:))))];
    for rowi=1:3
        for coli=1:MaxFreqCond
            % fdr
            MyV1=MyVList1(:,coli);
            MyV2=MyVList2(:,coli);
            Hafdr=bafdr(MyV2,0.05);
            % ebi
            Hebi=O.H(4,(1:128)+(coli-1).*128);
            MyVebiplot=MyV1;
            MyVebiplot(~Hebi)=NaN;
            
            ah(coli,rowi)=axes('position',mysubplotpos{coli,rowi},'XGrid','off','XMinorGrid','off');
            switch rowi
                case 1
                    btopoplot(MyVebiplot,ChXX,ChYY,ChXX,ChYY,myClim1,0.5,1,1);
                    title(['{\boldmath $' cell2mat(bnth(strsplit(fiLabelLIST{coli},' ('),1)) '$} \textbf{\textsf{(' cell2mat(bnth(strsplit(fiLabelLIST{coli},' ('),2)) '}}'],'FontSize',MyTitleFontSize,'Interpreter','latex');
                    set(gca,'FontSize',MyFontSize,'Box','off');%
                    colormap(ah(coli,rowi),DiffColorMap);
                case 2
                    PlotContour=(O.H(3,(1:128)+(coli-1).*128)).'.*1+(O.H(4,(1:128)+(coli-1).*128)).'.*1+(O.H(2,(1:128)+(coli-1).*128)).'.*1+(O.H(1,(1:128)+(coli-1).*128)).'.*1;
                    btopoplot(PlotContour,ChXX,ChYY,ChXX,ChYY,[0 4],0.5,1,1);
                    [~,sortedmeasures]=sort(O.ThresholdMatrix(1:4,1));
                    TempColorMatrix=[0 0.5 1;1 0 0;0 0.5 0;0 0 0]; % P1,FDR,beta,p1
                    colormap(ah(coli,rowi),([1 1 1;TempColorMatrix(sortedmeasures,:)]));
                case 3
                    btopoplot((MyV2<0.05)+Hafdr*1.0,ChXX,ChYY,ChXX,ChYY,[0 2],0.5,1,1);
                    colormap(ah(coli,rowi),[1;0.75;0.5]*[1 1 1])
            end
            set(ah(coli,rowi),'position',mysubplotpos{coli,rowi},'XGrid','off','XMinorGrid','off');
        end
        
        % colorbar 3
        if rowi==1
            coli=8;
            ah(coli,rowi)=axes('position',mysubplotpos{coli,rowi},'XGrid','off','XMinorGrid','off');
            pcolor([1 2],linspace(myClim1(1)+0.5,myClim1(2)+0.5,100).',linspace(myClim1(1)+0.5,myClim1(2)+0.5,100).'*[1 1]);
            colormap(ah(coli,rowi),DiffColorMap);
            shading flat;
            ylim([myClim1(1) myClim1(2)]+0.5)
            text(-1.65,0.5,0,'AUROC','FontSize',MyFontSize,'Rotation',90,'HorizontalAlignment','center');
            set(gca,'XTick',[],'YTick',[0.05 0.5 0.95],'YTickLabels',{num2str(0.05) num2str(0.5) num2str(0.95)},'FontSize',MyFontSize);
            set(ah(coli,rowi),'position',mysubplotpos{coli,rowi},'XGrid','off','XMinorGrid','off','layer','top');
            
            coli=9;
            ah(coli,rowi)=axes('position',mysubplotpos{coli,rowi},'XGrid','off','XMinorGrid','off');
            pcolor([1 2],[linspace(-myClim2(2),myClim2(1),100).';linspace(myClim2(1),myClim2(2),100).'],[linspace(-myClim2(2),myClim2(1),100).';linspace(myClim2(1),myClim2(2),100).']*[1 1]);
            colormap(ah(coli,rowi),[linspace(0,1,512),linspace(1,0,512)].'*[1 1 1]);
            shading flat;
            text(-1,0,0,'-log_{10}(p)','FontSize',MyFontSize,'Rotation',90,'HorizontalAlignment','center')
            ylim([-myClim2(2) myClim2(2)])
            set(gca,'XTick',[],'YTick',[-myClim2(2) 0 myClim2(2)],'YTickLabels',{num2str(myClim2(2)) num2str(0) num2str(myClim2(2))},'FontSize',MyFontSize);
            set(ah(coli,rowi),'position',mysubplotpos{coli,rowi},'XGrid','off','XMinorGrid','off','layer','top');
            
            coli=10;
            LegMarkerSize=15;
            ah(coli,rowi)=axes('position',mysubplotpos{coli,rowi},'XGrid','off','XMinorGrid','off');
            plot(O.Prlist,O.EstimatedFDR,'Color',[0.5 0 0],'LineWidth',1);
            hold on
            plot(O.Prlist,O.EstimatedBeta,'Color',[0 0.5 0],'LineWidth',1);
            plot(O.ThresholdMatrix(4,1).*[1 1],[0 1],'Color','k','LineWidth',2)
            ball1ah=plot(O.ThresholdMatrix(4,1),O.ThresholdMatrix(4,2),'Color',[0.5 0 0],'Marker','.','MarkerSize',LegMarkerSize,'LineStyle','none');
            ball2ah=plot(O.ThresholdMatrix(4,1),O.ThresholdMatrix(4,3),'Color',[0 0.5 0],'Marker','.','MarkerSize',LegMarkerSize,'LineStyle','none');
            xlim([0.5 1])
            ylim([0 0.5])
            title('EBI','FontSize',MyFontSize,'FontWeight','normal'); % ['p_1: ' num2str(O{1,7}.ThresholdMatrix(4,4),'%.2f')]
            if O.ThresholdMatrix(4,1)>=0.5
                set(gca,'YTick',[0 0.5],'XTick',[0.5 O.ThresholdMatrix(4,1) 1],'XTickLabel',{'0.5',['P_1 = ' num2str(O.ThresholdMatrix(4,1),2)],'1'})
            else
                set(gca,'YTick',[0 0.5],'XTick',[O.ThresholdMatrix(4,1) 0.5 1],'XTickLabel',{['P_1 = ' num2str(O.ThresholdMatrix(4,1),2)],'0.5','1'})
            end
            legend([ball1ah ball2ah],{['FDR = ' num2str(O.ThresholdMatrix(4,2),'%.2f')],['\beta = ' num2str(O.ThresholdMatrix(4,3),'%.2f')]},'Location','best')
            text(-0.07,0.5,0,'FDR','HorizontalAlignment','center','FontWeight','normal','FontSize',MyFontSize, 'Color', [0.5 0 0],'Rotation',90,'Units','normalized');
            text(1.05,0.5,0,'\beta','HorizontalAlignment','center','FontWeight','normal','FontSize',MyFontSize, 'Color', [0 0.5 0],'Rotation',-90,'Units','normalized');
            set(ah(coli,rowi),'position',mysubplotpos{coli,rowi},'XGrid','off','XMinorGrid','off');
        end
        
        if rowi==2          
            ThisTableColor=flipud([0 0.5 0;0 0 0;1 0 0;0 0.5 1]);
            coli=10;
            ah(coli,rowi)=axes('position',[mysubplotpos{coli-2,rowi}(1:2) mysubplotpos{coli,rowi}(3)*1.5 mysubplotpos{coli-2,rowi}(4)],'XGrid','off','XMinorGrid','off');
            myorder=[1 2 4 3];
            MyFontSizeS=6;
            tabd = O.ThresholdMatrix(myorder,myorder);
            cnames = {'{\boldmath$P_1$}','{\boldmath$FDR$}','{\boldmath$p_1$}','{\boldmath$\beta$}'};
            rnames = cellfun(@(x) ['{\boldmath$@$} ' x],cnames,'UniformOutput',0);
            del0=0.3;
            del=0.22;
            for ci=2:5
                text(del0+del*(ci-2),1,0,cnames(ci-1), 'HorizontalAlignment','left','FontWeight','bold','FontSize',MyFontSize, 'Color', [0 0 0],'Interpreter','latex');
                hold on
            end
            
            for rowic=2:5
                text(0,1-(rowic-2).*del-del0,0,rnames(rowic-1), 'HorizontalAlignment','left','FontWeight','bold','FontSize',MyFontSize, 'Color', ThisTableColor(rowic-1,:),'Interpreter','latex');
                for ci=2:5
                    if ci==rowic
                        myweight='bold';
                    else
                        myweight='normal';
                    end
                    text((ci-2)*del+del0,1-(rowic-2).*del-del0,0,num2str(tabd(rowic-1,ci-1),2), 'HorizontalAlignment','left','FontWeight',myweight,'FontSize',MyFontSizeS, 'Color', ThisTableColor(rowic-1,:));
                end
            end
            axis off
            set(ah(coli,rowi),'position',[mysubplotpos{coli-2,rowi}(1:2) mysubplotpos{coli,rowi}(3)*1.5 mysubplotpos{coli-2,rowi}(4)],'XGrid','off','XMinorGrid','off');
        end
        
        if rowi==3
            coli=10;
            ah(coli,rowi)=axes('position',mysubplotpos{coli,rowi},'XGrid','off','XMinorGrid','off');
            grah=zeros(1,2);
            grah(1)=plot(-2,-2,'Color',[1 1 1]*0.75,'LineWidth',4);
            hold on
            grah(2)=plot(-2,-2,'Color',[1 1 1]*0.5,'LineWidth',4);
            xlim([0 1])
            ylim([0 1])
            legend(grah,{'p < 0.5','q < 0.05 (aFDR)'},'Location','east')
            axis off
            set(ah(coli,rowi),'position',mysubplotpos{coli,rowi},'XGrid','off','XMinorGrid','off');
        end
        
    end
    
    titleah = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    set(titleah, 'YTick', [], 'XTick', [], 'YTickLabel', [], 'XTickLabel', [],'Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5,0.99,0,'Empirical Bayesian Inference of Motor-Related Changes in Experimntal Spectral EEG Power','HorizontalAlignment','center','VerticalAlignment','top','FontWeight','bold','FontSize',MyTitleFontSize+2, 'Color', [0 0 0]);
    
    tifffilename='MS_fig5_ExpSpectralPowerTopoplots';
    print(fh1,tifffilename,'-dpng', '-r600');
    disp([tifffilename ' was generated.'])
    close(fh1);
    warning('on','MATLAB:handle_graphics:exceptions:SceneNode')
end
end


%% Functions
function [ Outputs, cmdOut ] = locfdrWrap(zz, bre, df, pct, pct0, nulltype, type, plot, mult, mlests, main,sw) %#ok<*INUSD>
% Calls R and runs locfdr on the data and brings back the results to MATLAB

% debug syntax
% [ Outputs, cmdOut ] = locfdrWrap( [randn(1,180) 1+randn(1,10) -0.7+randn(1,10)], 10, 10, 0, [0.47 53], 1,  0, 0, [0.5 1 1.5],[0 1], 'label',2)

save('locfdrdatain.mat','zz', 'bre', 'df', 'pct', 'pct0', 'nulltype', 'type', 'plot', 'mult', 'mlests','main','sw','-v4');
[cmdStatus,cmdOut]=system('"C:\Program Files\R\R-3.4.4\bin\Rscript.exe" locfdrscript.r'); % needs to be adjusted to the R executable on the system
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

function y=createColormap(x)
bcolormapart4=nan(7,3);
bcolormapart4(1,:)=[0,0,0.5]; % Blue
bcolormapart4(2,:)=[0,0,1];
bcolormapart4(3,:)=[0,0.5,1];
bcolormapart4(4,:)=[0.5,1,1];
bcolormapart4(5,:)=[1,1,1]; % White
bcolormapart4(6,:)=[1,1,0.5];
bcolormapart4(7,:)=[1,0.5,0];
bcolormapart4(8,:)=[1,0,0];
bcolormapart4(9,:)=[0.5,0,0]; % Red
SetX=1:9;
XX=linspace(1,9,1024);
y=[interp1(SetX,bcolormapart4(:,1),XX).' interp1(SetX,bcolormapart4(:,2),XX).' interp1(SetX,bcolormapart4(:,3),XX).'];
end

function y=bnth(x,n)
y=x(n);
end

function [ positions ] = bsubplotpos(DX,DY,dx,dy,varargin)
% [ positions ] = subplot_pos(DX,DY,dx,dy)
if nargin==5
    plotwidth=1;
    plotheight=1;
else
    plotwidth=sum(DX)+sum(dx);
    plotheight=sum(DY)+sum(dy);
end
DY=DY(end:-1:1);
dy=dy(end:-1:1);
nx=length(DX);
ny=length(DY);
positions=cell(nx,ny);
for xi=1:nx
    for yi=1:ny
        positions{xi,ny-yi+1}=[(sum(DX(1:(xi-1)))+sum(dx(1:xi)))/plotwidth (sum(DY(1:(yi-1)))+sum(dy(1:yi)))/plotheight DX(xi)/plotwidth DY(yi)/plotheight];
    end
end
end

function []=btopoplot(V,ChXX,ChYY,ChXXbp,ChYYbp,myClim,varargin)
R=100;
guidewidth=0.025;
if nargin>6
    MapLineWidth=varargin{1};
else
    MapLineWidth=0.5;
end
if nargin>7
    BGcolour=false;
else
    BGcolour=true;
end
if nargin>8
    NaNPlot=false;
else
    NaNPlot=true;
end

% wider plot
if nargin>9
    plotRcut=105;
else
    plotRcut=100;
end
ti=(-plotRcut):1:plotRcut;
[qx,qy] = meshgrid(ti,ti);

Vn=V;%(V-nanmean(V))./nanmean(V);
VnNAN=Vn;
VnNAN(isnan(Vn))=0;
if NaNPlot
Vn(isnan(Vn))=0;
end
F = scatteredInterpolant(ChXXbp.', ChYYbp.',VnNAN);
Fnan = scatteredInterpolant(ChXXbp.', ChYYbp.',Vn,'nearest');
qz = F(qx,qy);
qznan = Fnan(qx,qy);
qz=qz+0.*qznan;

plotR=100;
qz((qx.^2+qy.^2)>(plotRcut.^2))=NaN;
if ~all(isnan(V)) || BGcolour
pcolor(ti,ti,qz); shading flat;caxis(myClim);
end
hold on
plottheta=0:0.001:(2.01*pi);
Yoffset=0;
xlim(115.*[-1 1]);
ylim(115.*[-1 1]);
plot(ChXX(1:128),ChYY(1:128),'Color',[1 1 1]-1,'Marker','.','MarkerSize',1,'LineStyle','none');
set(gca, 'YTick', [], 'XTick', [], 'YTickLabel', [], 'XTickLabel', []);
% Head Cuircle
hold on
plot(plotR.*cos(plottheta),plotR.*sin(plottheta)-Yoffset,'Color','k','LineWidth',MapLineWidth);
eartheta=(0.4*pi):0.001:(1.6*pi);
plot(-plotR+cos(eartheta)*10,25.*sin(eartheta)-Yoffset,'Color','k','LineWidth',MapLineWidth);
plot(plotR-cos(eartheta)*10,25.*sin(eartheta)-Yoffset,'Color','k','LineWidth',MapLineWidth);
plot(plotR.*cos(2*pi*(0.25+[-0.015 0 0.015])),plotR.*sin(2*pi*(0.25+[-0.015 0 0.015]))-Yoffset+[0 15 0],'Color','k','LineWidth',MapLineWidth);
% dots:
plot(ChXX(1:128),ChYY(1:128),'Color',[0 0 0],'LineStyle','none','Marker','.','MarkerSize',1);
plot(ChXX([1 19 54 85 115]),ChYY([1 19 54 85 115]),'Color',[0 0 0],'LineStyle','none','Marker','.','MarkerSize',3);
plot([-1 1].*R,[0 0],'LineWidth',guidewidth,'Color',[0 0 0]);
plot([0 0],[-1 1].*R,'LineWidth',guidewidth,'Color',[0 0 0]);
axis off
end

% The contents of the "locfdrscript.r" file (remove the %'s)

% library(R.matlab)
% library(locfdr)
% data <- readMat('locfdrdatain.mat')
% # data$zz
% # data$bre
% # data$df
% # data$pct
% # data$pct0
% # data$nulltype
% # data$type
% # data$plot
% # data$mult
% # data$mlests
% # data$main
% # data$sw
% #dataout <- locfdr(data$zz, bre = data$bre, df = data$df, pct = data$pct, pct0 = data$pct, nulltype = data$nulltype , type = data$type, plot = data$plot, data$mult, data$mlests, main = "label", sw=data$sw)
% dataout <- locfdr(data$zz)
% # dataout$mat
% writeMat("locfdrdataout.mat",fdr=dataout$fdr,fp0=dataout$fp0,Efdr=dataout$Efdr,cdf1=dataout$cdf1,mat=dataout$mat,z2=dataout$z.2,mult=dataout$mult)

%% 
% Note: the R software and the locfdr package need to be obtained
% separately and are not part of this toolbox. They are separately licensed
% under their terms and conditions. 

%%
% Copyright (c) 2018-2019 Bahman Nasseroleslami, All rights reserved.
%
% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
%
%     Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
%     Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
%     Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.