function [O]=ebi_uv(p_X,p_g,varargin)
% Power Calculation for Single Comparisons, 
% using emprirical bootstrapping-based inferece (ebi)
%
% usage:
% O = ebi_uv( X, g, ['Nboot0',2000,'Nboot1',2000,'Thresholds',[0.05 0.20]] )
%
% X: Input Data: Nobs x 1 Data Matrix
% g: Group Labels: Nobs x1 column vector of 0's for control and 1's for patient/treatment groups
% varargout: options: a structure with following fields
%                    options.Thresholds : a 1 x 2 vector of threshold values [Alpha, Beta] for Type I Error Alpha and Type II Error Beta, each between 0 and 1, at which hypotheses are tested and the corresponding values of other measures reported. Default: [0.05 0.2]
%                    options.Nboot0 : Number of null bootsrapping 
%                    options.Nboot1 : Number of non-null bootsrapping 
% O : Output:  Structure with multiple fileds listed below
% O.AUC : Area Under the Receiving O-perator Curve for each variable
% Posterior,Zlist,EstimatedPriors,EstimatedDensities,EstimatedPosteriors,AUC0b,P0AUC0b, Prlist, EstimatedFDR, EstimatedAlpha, EstimatedBeta,GlobalAlpha,GlobalBeta,AUC1b
% Posterior: Ngroup x Nvar Probability P(group | Data)

% O : Output:  Structure with multiple fileds listed below
% O.AUC : 1 x 1, Area Under the Receiving Operator Curve
% O.AUCZ : 1 x 1, AUC for the variable, transformed to Z space
% O.AUCV : 1 x 1, variance of AUC
% O.AUCZV : 1 x 1, variance of AUC, transformed to Z space
% O.Zlist, 1 x [] vector, Z-variable values in which PDFs are defined. 
% O.EstimatedDensities, 3 x [] , PDFs for the the null PDF and non-null PDF [f0;f1]
% O.AUC0b, Nboot0 x 1, Null bootstrapped test statistics values
% O.AUC1b, Nboot0 x 1, Non-Null bootstrapped test statistics values
% O.EstimatedFDR, 1 x 1001 False Discovery Rate (FDR) as a function of Threshold Pr
% O.EstimatedAlpha, 1 x 100001 Type I error (alpha), equally spaced values
% O.EstimatedBeta, 1 x 100001 Type II error (Beta) as a function of O.EstimatedAlpha
% O.ThresholdMatrix. Rows: @alpha, @beta, @AUC(l)=AUC(r)=AUC,@alpha=beta
%                    Columns: Alpha, Beta, AUCleft AUCright
% O.pvalue: p-value of the given test data
% O.Ngroups: number of data in each group
% O.Model0: null pdf object
% O.Model1: non-null pdf object

% Example:
% O = ebi_uv(X,g);
% O = ebi_uv(X,g,'Thresholds',[0.1 0.05],'Nboot0',1000,'Nboot1',1000)
%
% sample test/debug command:
% O=ebi_uv([randn(50,1);randn(50,1)+1.5],[zeros(50,1);ones(50,1)],'Thresholds',[0.005 0.05],'Nboot0',10000,'Nboot1',10000);
%
% Written by:
% Bahman Nasseroleslami, Trinity College Dublin, the University of Dublin, 27/05/2016, nasserob@tcd.ie, bahman@neuromotor.org
% Part of the Emprical Bayesian Inference (EBI) toolbox for MATLAB

p = inputParser;
addRequired(p,'p_X',@(X) isnumeric(X) && all(isfinite(X(:))) && size(X,1)>1 && size(X,2)==1);
addRequired(p,'p_g',@(X) all(isfinite(X(:))) && size(X,1)>1 && size(X,2)==1 && all((X==0) | (X==1)) && sum(X==0)>0 && sum(X==1)>0);
addParameter(p,'Thresholds',[0.05 0.2],@(X) isnumeric(X) && all(isfinite(X(:))) && size(X,1)==1 && size(X,2)==2 && all((X>0) & (X<1)));
addParameter(p,'Nboot1',2000,@(X) isnumeric(X) && isfinite(X) && all(size(X)==1) && (X>=0));
addParameter(p,'Nboot0',2000,@(X) isnumeric(X) && isfinite(X) && all(size(X)==1) && (X>=1));
parse(p,p_X,p_g,varargin{:});
X=p.Results.p_X;
g=p.Results.p_g;
Thresholds=p.Results.Thresholds;
Nboot1=p.Results.Nboot1;
Nboot0=p.Results.Nboot0;
assert(all(size(X)==size(g)),'X and g must have the same size.')

rng('shuffle');
[Nobs,~]=size(X);
Ngroups=[sum(g==0) sum(g>0)];
Prec=1e-3;
KernelSupport=[-20 20];
epsilon=2.061e-9;
Zlist=KernelSupport(1):Prec:KernelSupport(2);

% Take AUC as Zi variable
[AUC, AUCV]  = bAUROC(g,X(:,1));
AUCZ=A2Z(AUC,epsilon);
AUCZV=AUCV./(((AUC+epsilon).*(epsilon+1-AUC)).^2);

% Estimate Alternative Density by Boot or Model

switch Nboot1
    case 0
        Model1=makedist('Normal','mu',AUCZ,'sigma',sqrt(AUCZV));
    otherwise
        AUC1b=zeros(Nboot1,1);
        for bi=1:Nboot1
            %     n0=sum(g(randi(Nobs,1,Nobs))==0); Bad ! Biased!
            n0=min(max(sum(g==0),1),Nobs-1); % Good !
            PermID=[randi(Ngroups(1),1,n0) Ngroups(1)+randi(Ngroups(2),1,Nobs-n0)].';
            AUC1b(bi,1) = bAUROC(g,X(PermID,1));
        end
        AUC1b=A2Z(AUC1b,epsilon);
        Model1=bgmdfit(AUC1b(:));
end

% Estimate Null Density by Boot
AUC0b=zeros(Nboot0,1);
for bi=1:Nboot0
    g0=brandom(g,Nobs);
    AUC0b(bi,1) = bAUROC(g0,X(:,1)); % Pay attention: The labels should be permuted, not the data.
end
AUC0b=A2Z(AUC0b,epsilon);
Model0=bgmdfit(AUC0b(:));

% Fit pdfs to null and non-null distribution
    function y=f0(x)
        y=reshape(pdf(Model0,x(:)),size(x));
    end
    function y=f1(x)
        y=reshape(pdf(Model1,x(:)),size(x));
    end

f0dist=f0(Zlist);
f1dist=f1(Zlist);
EstimatedDensities=[f0dist;f1dist];

EstimatedAlpha=0:0.00001:1;
EstimatedBeta=nan(size(EstimatedAlpha));

MyCDF0=cdf(Model0,Zlist.');
MyCDF1=cdf(Model1,Zlist.');
for Pri=1:length(EstimatedAlpha)
    DetectionAUD=((1-MyCDF0)<=EstimatedAlpha(Pri)*0.5) | ((MyCDF0)<=EstimatedAlpha(Pri)*0.5);
    EstimatedBeta(Pri)=1-sum(EstimatedDensities(2,DetectionAUD))./sum(EstimatedDensities(2,:));
end

% Calculate Thresholds and test Hypoptheses
%               Alpha, Beta, [AUCleft AUCright]                      (Estimated Equivalent, Output)
% Alpha
% Beta
% AUC
% (Given Criterion, Input)
ThresholdMatrix=zeros(4,4);
[~,tid]=min(abs(EstimatedAlpha-Thresholds(1)));
[~,tid2]=min(abs(MyCDF0-Thresholds(1)*0.5));
[~,tid3]=min(abs((1-MyCDF0)-Thresholds(1)*0.5));
ThresholdMatrix(1,:)=[Thresholds(1), EstimatedBeta(tid), Z2A(Zlist(tid2),epsilon), Z2A(Zlist(tid3),epsilon)];

[~,tid]=min(abs(EstimatedBeta-Thresholds(2)));
[~,tid2]=min(abs(MyCDF1-Thresholds(2)*0.5));
[~,tid3]=min(abs((1-MyCDF1)-Thresholds(2)*0.5));
ThresholdMatrix(2,:)=[EstimatedAlpha(tid), Thresholds(2), Z2A(Zlist(tid2),epsilon), Z2A(Zlist(tid3),epsilon)];

prep=cdf(Model0,AUCZ);
pvalue=2.*min(prep,1-prep);
[~,tid]=min(abs(EstimatedAlpha-pvalue));
ThresholdMatrix(3,:)=[pvalue, EstimatedBeta(tid), AUC, AUC];

[~,eq]=min(abs(EstimatedAlpha-EstimatedBeta));
[~,tid2]=min(abs(MyCDF0-EstimatedAlpha(eq)*0.5));
[~,tid3]=min(abs((1-MyCDF0)-EstimatedAlpha(eq)*0.5));
ThresholdMatrix(4,:)=[EstimatedAlpha(eq), EstimatedBeta(eq), Z2A(Zlist(tid2),epsilon), Z2A(Zlist(tid3),epsilon)];
% GenerateOutput
O.AUC=AUC;
O.AUCZ=AUCZ;
O.AUCV=AUCV;
O.AUCZV=AUCZV;
O.Zlist=Zlist;
O.EstimatedDensities=EstimatedDensities;
O.AUC0b=AUC0b;
if Nboot1>0,  O.AUC1b=AUC1b; end
O.EstimatedAlpha=EstimatedAlpha;
O.EstimatedBeta=EstimatedBeta;
O.ThresholdMatrix=ThresholdMatrix;
O.pvalue=pvalue;
O.Ngroups=Ngroups;
O.Nvar=1;
O.Model0=Model0;
O.Model1=Model1;
end
function Y=A2Z(X,epsilon)
Y=0.5*(log(epsilon+X)-log(epsilon+1-X));
Yb=abs(Y)>9;
%Y=Y.*~Yb+(9.*sign(Y)+min(max(0.25.*randn(size(Y)),-0.5),0.5)).*Yb;
for nni=1:size(Y,1)
    Y(nni,Yb(nni,:))=9.*sign(Y(nni,Yb(nni,:)))+METNs(sum(Yb(nni,:),2));
end
end
function Y=Z2A(X,epsilon)
Y=((epsilon+1).*(1+exp(2*X)-(2*epsilon+1)))./(1+exp(2*X));
end
function [Model]=bgmdfit(W)
Nmodels=(floor(length(W)/3)-1);
GMM=cell(1,Nmodels);
GMMAIC=nan(1,Nmodels);
warning('off','stats:gmdistribution:FailedToConvergeReps');
DoLoop=true;
k=0;
while DoLoop
    k=k+1;
    GMM{k} = fitgmdist(W,k,'CovarianceType','diagonal','Options',statset('Display','off','MaxIter',500,'TolFun',1e-6),'RegularizationValue',0.001,'Replicates',5,'SharedCovariance',false,'Start','randSample');
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
function X=METNs(N)
d1=normcdf(-4);
d2=normcdf(4);
X = norminv(d1+((d2-d1)/(2*N)).*(1:2:(2*N)),0,0.25);
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

% Copyright (c) 2018-2019 Bahman Nasseroleslami, All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 
%     Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
%     Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
%     Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission. 
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
