function O=ebi_empiricalstat(x,X0,X1,varargin)
% Emprical estimation of p-values, and alpha-beta curve from null and non-null bootstrapping data (Single or Multiple Comparisons)
%
% Description: Takes the Test Statitic value(s) x in, as well as the null sample X0 data and
% non-null sample X1 data, then it calculates the empirical p-value, as
% well as apha-beta function and the values at x
% 
% Syntax:
% O=ebi_empiricalstat(x,X0,X1,['Thresholds',Thresholds],['Precision',Precision],['OutputPrecision',OutputPrecision])
%
% x : Input Data: 1 x N vector, Test statitic values to be tested (against null and non-null distribution)
% X0: Input Data: 1 x N vector, Data obtained from bootstrapping under Null hypothesis
% X1: Input Data: 1 x N vector, Data obtained from bootstrapping under Non-Null hypothesis
% varargin: options: Name - Value pairs
%                    Thresholds : a 1 x 2 vector of threshold values [Alpha, Beta] for  Type I Error Alpha and Type II Error Beta, each between 0 and 1,
%                                at which hypotheses are tested and the corresponding values of other measures reported. Default: [0.05 0.2]
%                     Precision : a 1 x 1 scalar, precision for p-value estimation (default: 1e-6)
%               OutputPrecision : a 1 x 1 scalar, precision for estimating alpha-beta curve (defualt: 1e-3)
%
% Examples:
% O=ebi_empricalstat(x,X0,X1)
% O=ebi_empricalstat(x,X0,X1,'Thresholds',[0.005 0.05])
%
% sample test/debug command:
% O=ebi_empricalstat([2 3 4],randn(1,2000),2+randn(1,2000))
%
% Written by:
% Bahman Nasseroleslami, Trinity College Dublin, the University of Dublin, 2016, nasseroleslami@gmail.com
% Revision: 4/8/2017.
% To Do: Add empirical FDR
% Part of the Emprical Bayesian Inference (EBI) toolbox for MATLAB

pip = inputParser;
addRequired(pip,'x',@(X) isnumeric(X) && all(isfinite(X(:))) && size(X,2)>0 && size(X,1)==1);
addRequired(pip,'X0',@(X) isnumeric(X) && all(isfinite(X(:))) && size(X,2)>1 && size(X,1)==1);
addRequired(pip,'X1',@(X) isnumeric(X) && all(isfinite(X(:))) && size(X,2)>1 && size(X,1)==1);
addParameter(pip,'Thresholds',[0.05 0.2],@(X) isnumeric(X) && all(isfinite(X(:))) && size(X,1)==1 && size(X,2)==2 && all((X>0) & (X<1)));
addParameter(pip,'Precision',1e-6,@(X) isnumeric(X) && isfinite(X) && all(size(X)==1) && (X>=0));
addParameter(pip,'OutputPrecision',1e-3,@(X) isnumeric(X) && isfinite(X) && all(size(X)==1) && (X>=0));
parse(pip,x,X0,X1,varargin{:});
x=pip.Results.x;
X0=pip.Results.X0;
X1=pip.Results.X1;
Thresholds=pip.Results.Thresholds;
qd=pip.Results.Precision;
qd2=pip.Results.OutputPrecision;
%assert(all(size(X)==size(g)),'X and g must have the same size.')

X0=X0(:).'; % Column
X1=X1(:).'; % Column
x=x(:); % Row
Q=0:qd:1;
[~,Pid]=min(abs(quantile(X0,Q)-x),[],2);
P0=(Pid-1).*qd;
P=min([P0,1-P0],[],2)*2; % 2-sided P

% Power at Alpha
Alpha=0:qd2:1;

TempV1=quantile(X0,Alpha/2);
TempV2=quantile(X0,1-(Alpha/2));
TempVV=sort([TempV1; TempV2],1);
TempVV1=TempVV(1,:);
TempVV2=TempVV(2,:);
[~,Pidb1]=min(abs(quantile(X1,Q)-TempVV1.'),[],2);
[~,Pidb2]=min(abs(quantile(X1,Q)-TempVV2.'),[],2);
Beta=(Pidb2-Pidb1).*qd;

% Matrix
% Calculate Thresholds and test Hypoptheses
%               Alpha, Beta, [statleft Statright]                      (Estimated Equivalent, Output)
% Alpha
% Beta
% TestStat X_x
% (Given Criterion, Input)
ThresholdMatrix=zeros(4,4);
[~,TempID]=min(abs(Alpha-Thresholds(1)));
ThresholdMatrix(1,:)=[Thresholds(1) Beta(TempID) TempVV1(TempID) TempVV2(TempID)];

[~,TempID]=min(abs(Beta-Thresholds(2)));
ThresholdMatrix(2,:)=[Thresholds(2) Alpha(TempID) TempVV1(TempID) TempVV2(TempID)];

% [~,P3id]=min(abs(quantile(X0,Q)-Thresholds(3)),[],2);
% P30=(P3id-1).*qd;
% P3=min([P30,1-P30],[],2)*2; % 2-sided P
% [~,TempID]=min(abs(Alpha-P3));
% ThresholdMatrix(3,:)=[P3 Beta(TempID) TempVV1(TempID) TempVV2(TempID)];

%%%
O.Alpha=Alpha;
O.Beta=Beta;
O.Power=1-Beta;
O.P=P;
O.ThresholdMatrix=ThresholdMatrix;

% Copyright (c) 2018 Bahman Nasseroleslami, All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 
%     Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
%     Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
%     Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission. 
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.