function [O]=ebi(p_X,p_g,varargin)
% Empirical Bayesian Inference for Multiple Comparisons
%
% Usage:
% O = ebi(X,g,[options])
%
% X: Input Data: Nobs x Nvar Data Matrix (Nobs: Number of observations, Nvar: number of varibles)
% g: Group Labels: Nobs x 1 column vector of 0's for control and 1's for patient/treatment groups
% varargin: options: Name - Value pairs
%                    Thresholds : a 1 x 4 vector of threshold values [P1, FDR, Beta, Alpha] for Posterior Probability
%                                of class 1 (P1), False discovery Rate (FDR), Type II Error Beta and Type I Error Alpha, each between 0 and 1,
%                                at which hypotheses are tested and the corresponding values of other measures
%                                reported. Default: [0.9 0.05 0.2 0.05]
%                          Test: 'AUCZ' (default): AUC as test statitic, transformed automatically to Z-scores
%                                'T': independent sample t, as test statistic
%                                'R': Spearman's rank correlation, transformed to Z-scores,
%                                'WZ': Wilcoxon's (paired) Sign-Rank W-statitic, transformed to Z
%                                'CohZ': Classic Spectral Coherence (Sxy^2/(Sxx.Syy)), the p-value transformed to Z space.
%                     Estimator: 'GMD': (default, recommended) Gaussian
%                                        Mixture Distribution (GMD)-based estimation of
%                                        probability density functions (PDFs), using Akaike
%                                        Information Criterion (AIK), 
%                                        (the only reliable option kept) 
%                         'pct': The percentage range of the null
%                                distribution for calculation of propr probability
%                                p_0. Default: [0.5 0.5] (gives the median)
%
% O : Output:  Structure with multiple fileds listed below
% O.AUC : 1 x Nvar, Area Under the Receiving Operator Curve for each variable
% O.AUCZ : 1 x Nvar, AUC for each variable, transformed to Z space
% O.Posterior : 2 x Nvar, Posterior Probability P(effect|Data). Row 1 is P_0(z), row 2 is P_1(z) for each variable
% O.Zlist, 1 x [] vector, Z-variable values in which PDFs are defined. 
% O.EstimatedPriors  , 1 x 2, [p0 p1], prior probabilities of the null and effect estimated from data
% O.EstimatedDensities, 3 x [] , PDFs for the the null PDF, non-null PDF, and mixed data PDF [f0;f1;f]
% O.EstimatedPosteriors, 2 x [], PDFs for the Posterior Probabilities [P0(z);P1(z)]
% O.AUC0b, Nboot0 x Nvar, Null bootstrapped test statistics values
% O.P0AUC0b
% O.Prlist, 1 x 1001, Pr_1 Threshold variable values in which alpha, beta, fdr are deined
% O.EstimatedFDR, 1 x 1001 False Discovery Rate (FDR) as a function of Threshold Pr
% O.EstimatedAlpha, 1 x 1001 Type I error (alpha) as a function of Threshold Pr
% O.EstimatedBeta, 1 x 1001 Type II error (Beta) as a function of Threshold Pr
% O.GlobalAlpha, 1 x 1 Scalar, Gloabal Type I error regarless of threshold (a measure of overlap of the null pdf f0(z) onto Posterior P1(z))
% O.GlobalBeta, 1 x 1 Scalar, Gloabal Type II error regarless of threshold (a measure of overlap of the non-null pdf f1(z) onto Posterior P0(z))
% O.AUC1b, [not used, reserved for future], Nboot1 x Nvar, Non-Null bootstrapped test statistics values

% Example:
% O = ebi(X,g);
% O = ebi(X,g,'Test','WZ','Thresholds',[0.8 0.1 0.2])
%
% sample test/debug command:
% O = ebi([randn(50,200);randn(50,100),randn(50,100)+1.5],[zeros(50,1);ones(50,1)]);
%
% Written by:
% Bahman Nasseroleslami, Trinity College Dublin, the University of Dublin, 27/05/2016, nasserob@tcd.ie, bahman@neuromotor.org
% based on (Efron, Tibshirani, Storey, Tusher, 2001, Journal of American Statistical Association) and (Efron 2007, The Annals of Statistics).
% Part of the Emprical Bayesian Inference (EBI) toolbox for MATLAB
% Revision: 13/6/2018. Fixed a bug in T-Test option. Added support for
% coherence analysis on raw complex-valued data

% Parse the inputs, specify the defaults, and assert correct formats
p = inputParser;
addRequired(p,'p_X',@(X) isnumeric(X) && all(isfinite(X(:))) && all(size(X)>1));
addRequired(p,'p_g',@(X) all(isfinite(X(:))) && size(X,1)>1 && size(X,2)==1 && all((X==0) | (X==1)) && sum(X==0)>0 && sum(X==1)>0);
addParameter(p,'Thresholds',[0.9 0.05 0.2 0.05],@(X) isnumeric(X) && all(isfinite(X(:))) && all(size(X)==[1 4]) && all((X>0) & (X<1))); % [P1Z, FDR, Beta, Alpha]
addParameter(p,'Test','AUCZ',@(X) ischar(X) && ismember(X,{'AUC','AUCZ','T','R','WZ','CohZ'}));
addParameter(p,'Estimator','GMD',@(X) ischar(X) && ismember(X,{'GMD','Kernel'}));
addParameter(p,'pct',[0.5 0.5],@(X) isnumeric(X) && all(isfinite(X)) && all(size(X)==[1 2]) && all((X>0) & (X<1)) && X(2)>=X(1));
% addParameter(p,'SigmaW2',1,@(X) isnumeric(X) && all(isfinite(X)) &&
% size(X,1)==1 && all(X>0) && all(X<1)); % reserved for advanced coherence anlysis
parse(p,p_X,p_g,varargin{:});
X=p.Results.p_X;
g=p.Results.p_g;
Thresholds=p.Results.Thresholds;
Test=p.Results.Test;
pct=p.Results.pct;
% SigmaW2=p.Results.SigmaW2; % reserved for advanced coherence analysis 
assert(size(X,1)==size(g,1),'X and g must have the same number of rows.')
assert( (~strcmp(Test,'R')) || (strcmp(Test,'R') && (mod(length(g),2)==0) && all(g(1:(length(g)/2))==0) && all(g((1:(length(g)/2))+(length(g)/2))==1)),'For Correlation, the first half of g must be 0s and second half of it must be 1s, in the corresponding order.')
assert( (~strcmp(Test,'WZ')) || (strcmp(Test,'WZ') && (mod(length(g),2)==0) && all(g(1:(length(g)/2))==0) && all(g((1:(length(g)/2))+(length(g)/2))==1)),'For pair-wise comparison using Wilcoxon''s Signed Rank, the first half of g must be 0s and second half of it must be 1s, in the corresponding order.')

rng('shuffle'); % Initialise random number generation

[Nobs,Nvar]=size(X); % Get number of observations and number of variables from input data.

Prec=1e-3; % Increment/Precision
KernelSupport=[-10 10]; % The range for test statitic
Zlist=-20:Prec:20; % the range for PDF and integration
epsilon=2.061e-9; % The epsilon value for modification of the Fisher-Type Transform so that the 0:1 or -1:1rnage is mapped to the -10:10 rather than infinity. 4.55e-5;

Ngroups=[sum(g==0) sum(g>0)]; % Number of observations in the group A (control) and group B (patients)

% Take AUC as Zi variable
% Calculate the test statistic and bring it to the z-domain if needed (also
% suppress the outliers
switch Test
    case 'T'
        [~,~,~,ttest2stat]=ttest2(X(g==0,:),X(g>0,:),'Vartype','unequal');
        AUC=icdf('Normal',cdf('T',ttest2stat.tstat,ttest2stat.df),0,1);
        AUCZ=SmashT(AUC);
    case 'AUCZ'
        AUC=0.5.*ones(1,Nvar); % Initialise as neutral values
        for vi=1:Nvar,     AUC(1,vi)  = bAUROC(g,X(:,vi));     end % calculate
        AUCZ=A2Z(AUC,epsilon); % bring to Z domain
    case 'WZ'
        AUC=0.5.*ones(1,Nvar); % Initialise as neutral values
        for vi=1:Nvar,     AUC(1,vi)  = Wstat(g,X(:,vi));     end
        AUCZ=A2Z(AUC,epsilon);
    case 'R'
        AUC=zeros(1,Nvar); % Initialise as neutral values
        for vi=1:Nvar,     AUC(1,vi)  = 0.5+0.5*corr(X(g==0,vi),X(g>0,vi),'type','Spearman','rows','pairwise'); end
        AUCZ=A2Z(AUC,epsilon);
    case 'CohZ' % under development 
        % use algebraic approximation
        % % C2=(abs(mean(X(g==0,:).*conj(X(g>0,:)),1)).^2)./(mean(abs(X(g==0,:)).^2,1).*mean(abs(X(g>0,:)).^2));
        % % TempP=(1-C2).^(((Nobs/2)-1).*SigmaW2); 
        % % use Hotelling's one-sample T^2
        %TempP=ones(1,Nvar);
        %for vi=1:Nvar, TempP(1,vi)=bhotellingT(complex(X(g==0,vi).*conj(X(g>0,vi)))); end
        %AUC=norminv(1-TempP); % norminv(2*max(1-(TempP),TempP)); 
        %AUCZ=SmashT(AUC);
end

ModelMixed=bgmdfit(AUCZ(:)); % fit GMD model to mixed data f(z)
ModelMixedY=reshape(pdf(ModelMixed,Zlist.'),size(Zlist)); % get the numerical PDF of the f(z)

% Estimate Null Density by Boot
Nboot0=max(20,fix(10000/Nvar)); % Number of bootstraps. For 500+ variables 20 bootsraps are enough (see Efron's paper), for 500- variables more bootsraps are done, for 100 variables, 100 bootsraps, for 10 variables 1000 bootsraps.
AUC0b=zeros(Nboot0,Nvar); % allocate the memory for bootstrap test statitics

% Follow the appropriate null-bootstrapping algorithm for each type of test
switch Test
    case 'T'
        for bi=1:Nboot0
            g0=brandom(g,Nobs);
            for vi=1:Nvar
                [~,~,~,ttest2stat]=ttest2(X(g0==0,vi),X(g0>0,vi),'Vartype','unequal');
                AUC0b(bi,vi)=icdf('Normal',cdf('T',ttest2stat.tstat,ttest2stat.df),0,1);
            end
        end
        AUC0b=SmashT(AUC0b);
    case 'AUCZ'
        for bi=1:Nboot0
            g0=brandom(g,Nobs);
            for vi=1:Nvar
                AUC0b(bi,vi) = bAUROC(g0,X(:,vi)); % Pay attention: The labels should be permuted, not the data.
            end
        end
        AUC0b=A2Z(AUC0b,epsilon);
    case 'WZ'
        for bi=1:Nboot0
            g0=2*randi(2,Ngroups(1),1)-3;
            for vi=1:Nvar
                AUC0b(bi,vi) = Wstat([zeros(Ngroups(1),1);ones(Ngroups(1),1)],X(:,vi).*[g0;g0]); % Pay attention: The labels should be permuted, not the data.
            end
        end
        AUC0b=A2Z(AUC0b,epsilon);
    case 'R'
        for bi=1:Nboot0
            g0=randi(Ngroups(1),Ngroups(1),1);
            g1=randi(Ngroups(1),Ngroups(1),1)+Ngroups(1);
            for vi=1:Nvar
                AUC0b(bi,vi) = 0.5+0.5*corr(X(g0,vi),X(g1,vi),'type','Spearman','rows','pairwise');
            end
        end
        AUC0b=A2Z(AUC0b,epsilon);
        %     case 'CohZ'
        %         for bi=1:Nboot0
        %             g0=randi(Ngroups(1),Ngroups(1),1);
        %             g1=randi(Ngroups(1),Ngroups(1),1)+Ngroups(1);
        %             %C2=(abs(mean(X(g0,:).*conj(X(g1,:)),1)).^2)./(mean(abs(X(g0,:)).^2,1).*mean(abs(X(g1,:)).^2));
        %             %TempP=(1-C2).^(((Nobs/2)-1).*SigmaW2);
        %             TempP=ones(1,Nvar);
        %             for vi=1:Nvar, TempP(1,vi)=bhotellingT(complex(X(g0,vi).*conj(X(g1,vi)))); end
        %             AUC0b(bi,:)=norminv(1-TempP);
        %         end
        %         AUC0b=SmashT(AUC0b);
end

% Sparse Null Data to avoid unnecessaru computation and increase the speed
if numel(AUC0b)>length(Zlist)
    sk=floor(numel(AUC0b)./(length(Zlist)-1));
    AUC0bs=bnth(sort(AUC0b(:)),(1+floor(0.5*sk)):sk:numel(AUC0b));
else
    AUC0bs=AUC0b;
end

% Find the null PDF f_0(z) and its numerical versio
Model0=bgmdfit(AUC0bs(:)); % Still the best, but slow.
Model0Y=reshape(pdf(Model0,Zlist.'),size(Zlist));

% p0: calculate the prior probability
[~,idp]=min(abs(bsxfun(@minus,[1; 1]*cumsum(Model0Y,2)./sum(Model0Y,2),pct.')),[],2);
p0max=min(sum(ModelMixedY(idp(1):idp(2)))./sum(Model0Y(idp(1):idp(2))),1);

% Calculate the Posterior probabilities for each test-statistic value of the given data and variables using the Bayesian formula
Posterior(1,:)=min(1,interp1(Zlist,Model0Y,AUCZ,'linear')./interp1(Zlist,ModelMixedY,AUCZ,'linear').*p0max);
Posterior(2,:)=1-Posterior(1,:);

% Store the prior values
EstimatedPriors=[p0max,1-p0max];

% Calculate the numerical PDF estimates of the Posterior probabilities using the Bayesian formula
EstimatedPosteriors(1,:)=min(1,Model0Y./ModelMixedY.*p0max);
EstimatedPosteriors(2,:)=1-EstimatedPosteriors(1,:);

% Calculate the numerical non-null PDF f_1(z)
switch EstimatedPriors(1)
    case 1
        f1dist=Model0Y;
    otherwise
        f1dist=EstimatedPosteriors(2,:).*ModelMixedY./EstimatedPriors(2);
end
f1dist=f1dist./sum(f1dist).*length(f1dist)./(Zlist(end)-Zlist(1)); % re-normalise for more accuracy.
EstimatedDensities=[Model0Y;f1dist;ModelMixedY];

% Estimate measures as a function of Threshold Probability
Prlist=0:0.001:1; % Threshold Probability
EstimatedFDR=nan(size(Prlist));
EstimatedAlpha=nan(size(Prlist));
EstimatedBeta=nan(size(Prlist));
for Pri=1:length(Prlist)
    DetectionAUD=EstimatedPosteriors(2,:)>=Prlist(Pri); % Theoretical Detection Region: DecisionPr=Prlist(Pri);
    if all(~DetectionAUD)
        EstimatedFDR(Pri)=0;
    else
        EstimatedFDR(Pri)=mean(EstimatedPosteriors(1,DetectionAUD).*EstimatedDensities(3,DetectionAUD))./mean(EstimatedDensities(3,DetectionAUD)); % Integral Formula over Density Functions for in the corresponding detection region.
    end
    EstimatedAlpha(Pri)=sum(EstimatedDensities(1,DetectionAUD))./sum(EstimatedDensities(1,:)); %OK! % Classic Frequentist alpha at each decision criterion
    EstimatedBeta(Pri)=1-sum(EstimatedDensities(2,DetectionAUD))./sum(EstimatedDensities(2,:)); % [consistent with counted version] Classic Frequentist beta at each decision criterion
end

GlobalBeta=mean(EstimatedPosteriors(1,:).*EstimatedDensities(2,:))./mean(EstimatedDensities(2,:));
GlobalAlpha=mean(EstimatedPosteriors(2,:).*EstimatedDensities(1,:))./mean(EstimatedDensities(1,:));

% Calculate Thresholds and test Hypoptheses
%               P1Z (Posterior), FDR, Beta, p1, Alpha                        (Estimated Equivalent, Output)
% P1Z
% FDR
% Beta
% p1
% Alpha
% (Given Criterion, Input)

% Calculate the measure at given thresholds
[~,tid]=min(abs(Prlist-Thresholds(1)));
H(1,:)=Posterior(2,:)>=Thresholds(1);
ThresholdMatrix(1,:)=[Thresholds(1), EstimatedFDR(tid), EstimatedBeta(tid), sum(H(1,:))./Nvar, EstimatedAlpha(tid)];

[~,tid]=min(abs(EstimatedFDR-Thresholds(2)));
H(2,:)=Posterior(2,:)>=Prlist(tid);
ThresholdMatrix(2,:)=[Prlist(tid),Thresholds(2), EstimatedBeta(tid), sum(H(2,:))./Nvar, EstimatedAlpha(tid)];

[~,tid]=min(abs(EstimatedBeta-Thresholds(3)));
H(3,:)=Posterior(2,:)>=Prlist(tid);
ThresholdMatrix(3,:)=[Prlist(tid), EstimatedFDR(tid),Thresholds(3), sum(H(3,:))./Nvar, EstimatedAlpha(tid)];

MyP=prctile(Posterior(2,:),EstimatedPriors(1).*100);
 [~,tid]=min(abs(Prlist-MyP));
switch EstimatedPriors(1)
    case 0
        H(4,:)=1;
        ThresholdMatrix(4,:)=[MyP, EstimatedFDR(tid),EstimatedBeta(tid), EstimatedPriors(2), EstimatedAlpha(tid)];
    case 1
        H(4,:)=0;
        ThresholdMatrix(4,:)=[MyP,0,1, EstimatedPriors(2), 0];
    otherwise
        H(4,:)=Posterior(2,:)>=MyP;
        ThresholdMatrix(4,:)=[MyP, EstimatedFDR(tid),EstimatedBeta(tid), EstimatedPriors(2), EstimatedAlpha(tid)];
end

[~,tid]=min(abs(EstimatedAlpha-Thresholds(4)));
H(5,:)=Posterior(2,:)>=Prlist(tid);
ThresholdMatrix(5,:)=[Prlist(tid), EstimatedFDR(tid),EstimatedBeta(tid), sum(H(5,:))./Nvar, Thresholds(4)];

% GenerateOutput
O.AUC=AUC;
O.AUCZ=AUCZ;
O.Zlist=Zlist;
O.EstimatedDensities=EstimatedDensities;
O.AUC0b=AUC0b;
O.Posterior=Posterior;
O.EstimatedPriors=EstimatedPriors;
O.EstimatedPosteriors=EstimatedPosteriors;
O.Prlist=Prlist;
O.EstimatedFDR=EstimatedFDR;
O.EstimatedAlpha=EstimatedAlpha;
O.EstimatedBeta=EstimatedBeta;
O.GlobalAlpha=GlobalAlpha;
O.GlobalBeta=GlobalBeta;
O.ThresholdMatrix=ThresholdMatrix;
O.H=H;
O.ThresholdLabels={'Posterior P1','FDR','Beta','prior p1','Alpha'};
O.Model0=Model0;
O.ModelMixed=ModelMixed;
O.KernelSupport=KernelSupport;
O.Nvar=Nvar;
O.Ngroups=Ngroups;
end
function Y=A2Z(X,epsilon)
Y=0.5*(log(epsilon+X)-log(epsilon+1-X));
Yb=abs(Y)>9;
for nni=1:size(Y,1)
    Y(nni,Yb(nni,:))=9.*sign(Y(nni,Yb(nni,:)))+METNs(sum(Yb(nni,:),2));
end
end
% Truncate the t-statistics between -9 and 9
function Y=SmashT(X)
X(X==-Inf)=-10;
X(X==Inf)=10;
Yb=abs(X)>9;
Y=X;
for nni=1:size(X,1)
    Y(nni,Yb(nni,:))=9.*sign(X(nni,Yb(nni,:)))+METNs(sum(Yb(nni,:),2));
end
end
% Fit GMDs using AIC 
function [Model]=bgmdfit(W)
Nrep=5;
Nmodels=(floor(length(W)/3)-1);
GMM=cell(1,Nmodels);
GMMAIC=nan(1,Nmodels);
warning('off','stats:gmdistribution:FailedToConvergeReps');
DoLoop=true;
k=0;
while DoLoop
    k=k+1;
    GMM{k} = fitgmdist(W,k,'CovarianceType','diagonal','Options',statset('Display','off','MaxIter',500,'TolFun',1e-6),'RegularizationValue',0.001,'Replicates',Nrep,'SharedCovariance',false,'Start','randSample');
    %     GMM{k} = fitgmdist(W,k,'CovarianceType','diagonal','Options',statset('Display','off'),'RegularizationValue',0.001,'Replicates',20,'SharedCovariance',false,'Start','randSample');
    GMMAIC(k)=GMM{k}.AIC;
    if k>5 && ((k>Nmodels) || all(0<diff(GMMAIC(max(1,k-3):k)))) % Stop when adding more mixtures worsens AIC for 3 consequetive times.
        DoLoop=false;
    end
end
warning('on','stats:gmdistribution:FailedToConvergeReps');
[~,Best]=nanmin(GMMAIC);
Model=GMM{Best};
end
function g0=brandom(g,Nobs)
k=0;
while k==0
    g0=g(randi(Nobs,Nobs,1),1);
    k=sum(g0==0).*sum(g>0);
end
end
% Wilcoxon's Sign-Rank W-statitic
function W=Wstat(g,X)
% Calculated the Normalised Wilcoxon's Signed Rank Test Statistic W. 
% W =Wstat(g,X)
% g: column vector N x 1, labels 0 or 1 for group 0 and group 1,
% X: column vector N x 1, Data
% W: W test statistic (normalised
XX=X(g>0)-X(g==0);
XX2=XX(XX~=0);
N=length(XX2);
[~,O]=sort(abs(XX2));
Y=XX2(O);
W=sum((1:N).'.*(Y>0))./sum(1:N);
end

function [A,varargout]=bAUROC(g,X)
% Calculated the Area Under the Receiver Operating Curve (AUROC) and it's variance using De Long's Method (Zhou 2009, Statistical Methods in Diagnostic Medicine)
% [A [, VA]]=bAUROC(g,X)
% g: column vector N x 1, labels 0 or 1 for group 0 and group 1,
% X: column vector N x 1, Data

% A: Area under the curve
% VA: Estimated Variance
n0=sum(g==0);
n1=sum(g==1);
T1=repmat(X(g==1),1,n0);
T0=repmat(X(g==0).',n1,1);
Phi=(T1>T0)+0.5.*(T1==T0);
V10=mean(Phi,2);
V01=mean(Phi,1);
A=mean(V10,1);
if nargout>1
    varargout{1}=(var(V10,0,1)/n1)+(var(V01,0,2)/n0);
end
end
% Function to create N deterministic data points with truncated normal
% distribution, bounded between [-1 1]
function X=METNs(N)
d1=normcdf(-4);
d2=normcdf(4);
X = norminv(d1+((d2-d1)/(2*N)).*(1:2:(2*N)),0,0.25);
end
% Take the n'th elements
function y=bnth(x,n)
y=x(n);
end

% Copyright (c) 2018-2019 Bahman Nasseroleslami, All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 
%     Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
%     Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
%     Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission. 
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.