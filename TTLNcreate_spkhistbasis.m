function [spkhist_basis] = TTLNcreate_spkhistbasis(len,ncoef,overlap)
% [spkhist_basis] = TTLNcreate_spikehistbasis(len,ncoef,overlap)
% create matrix of spkhist basis vectors
% assumed to be hanning functions on a log scale
% len is range of lags to cover (in bins)
% spike impulse reponse f = spkhist_basis*nim.spk_hist.coefs;
% H = Xspkhist*f
% For gradient calculation dHdf = spkhist_basis

tlag = 1:len;
tlaglog = log(tlag);
spkhist_basis = zeros(length(tlag),ncoef);

winwidth = (log(len+.5)-log(.5))/(1+(1-overlap)*(ncoef-1));
leftedges = log(.5)+(0:(ncoef-1))*(1-overlap)*winwidth;

for ii=1:ncoef
    idx = find(tlaglog>=leftedges(ii) & tlaglog<=leftedges(ii)+winwidth);
    spkhist_basis(idx,ii) = 1/2+(1/2)*cos(2*pi*((tlaglog(idx)-leftedges(ii))/winwidth-1/2));
end