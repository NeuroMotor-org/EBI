function [H,varargout] = bafdr(P, q)
% Adaptive False Discovery Rate, based on the algorithm in: 
% Benjamini, Krieger, Yekutieli, 2006, Biometrika, DOI: 10.1093/biomet/93.3.491
%
% Usage:
% H = bfdr(P, q)
% [H, k] = bfdr(P, q)
% [H, k, Pcr] = bfdr(P, q)

% P: m x n x ... (arbitrary dimension and shape) vector/matrix of p-values
% q: False Discovery Rate, (e.g. 0.05 or a value comparable to sig. level ?)
% 
% H: m x n x ... array (vector/matrix) of 0's and 1's.
% k: The number of significant detections (out of m x n x .. elements of P)
% Pcr: Critical P-value threshold used for accepting/rejecting P's hypotheses. 
% 
% Description
% The function corrects the p-values in array P (any dimensions and shapes),
% according to the FDR value q. The elements of H show the significant (1) 
% or non significant difference (0) for rejecting null hypotheses H_0. 
% the number of significant detections (k) and the 
% criterion p-value Pcr=k*q/m are optional outputs. 

% Written by Bahman Nasseroleslami, nasserob@tcd.ie, bahman@neuromotor.org

qs=q./(1+q);
m=sum(~isnan(P(:)));
[H, r1, pthresh]=LinearSetup(P,qs);

if r1~=0 && r1~=m
    m0=m-r1;
    [H,r1,pthresh]=LinearSetup(P,qs.*m./m0);
end

if nargout>1
    varargout{1}=r1;
end
if nargout>2
    varargout{2}=pthresh;
end

end % Function End

function [H,r,pthresh]=LinearSetup(P,q)
Pvec=P(:);
Ps=sort(Pvec(~isnan(Pvec)),'ascend');
m=length(Ps);
r=find(Ps<=((1:m).'./m*q),1,'last');

if isempty(r)
    H=false(size(P));
    r=0;
    pthresh=0;
else
    pthresh=(r./m.*q);
    H=(P<=pthresh);
end
end

% Copyright© 2015-2016 Bahman Nasseroleslami
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.