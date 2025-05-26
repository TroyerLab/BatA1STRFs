function spkhist = TTLNcreate_spkhist(len,ncoef,overlap)
% spkhist = TTLNcreate_spkhist(spkobs,nSpikehist)
% make spike history structure

% INPUTS: 
%     len - max number of lags in spike history
%     ncoef - number of coefs
%     overlap - overlap fraction in hanning basis
% OUTPUTS: 
%     spkhist: struct of spkhistory parameters
%       above values plus basismat
% If len = ncoef, then assume individual parameters, basis = []
%     alternatively, use overlapping hanning windows on a log scale
%     if overlap is a matrix then set basis = overlap

if nargin < 2 || isempty(ncoef)
    ncoef = len; 
end
if nargin < 3 || isempty(overlap)
    overlap = .5; %defualt overlap is 50%
elseif prod(size(overlap))>1
    spkhist.basis = overlap;
end
        

spkhist.len = len;
spkhist.ncoef = ncoef;
spkhist.coef = zeros(ncoef,1);

if len == ncoef
    spkhist.basis = [];
else
    spkhist.basis = TTLNcreate_spkhistbasis(len,ncoef,overlap);
end

