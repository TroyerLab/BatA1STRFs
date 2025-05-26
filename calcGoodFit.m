% calculate goodness of fit values

summaryfilename = 'summary.mat';

% LL - Log-likelihood per spike of data
% STACC - cross correlation of STA from 1st vs 2nd half
% STRFSTACC - cross correlation of STRF and STA for whole data (model
%      consistency measure)
% Gstd - standard deviation of input from STRF
% GvHfrac - fraction of variance of G vs. H = varG/(VarG+varH)
% GvTotfrac - fraction of variance of G vs. total: G+H = varG/Var(G+H)% NLfit - explained variance in binned firing probability based on NLfit
% H50, H90,  - 50th/90th percentile of cumulative spike
%       filter (absolute value) (msec)
% NLvbg - fraction of spikes 'explained' by constant bg firing rate vs.
%    rectified linear component
% NLfit_rms, NLfit_r2,  - fit of firing rate output to binned estimate of
%   firing rate, quantified as rMSE (Hz), or r2


initdata = 1;
savegoodfit = 1;

load('summary.mat','STRFSTACCall','STRFioall');

STRFinfo = readSTRFinfo('STRFinfo.xlsx');
N = length(STRFinfo.bat);

%LL
nLL = fitnLLall;
% input magnitudes
Gstd = STRFioall.Gstd;
GvHfrac = (STRFioall.Gstd.^2)./((STRFioall.Gstd.^2)+(STRFioall.Hstd.^2));
GvTotfrac = (STRFioall.Gstd.^2)./(STRFioall.GHstd.^2);
% spike history timing
dt = 2;
% NL fits
tmpNL.type = 'linthresh_bg';
for ii =1:N
    bins = STRFioall.bincent(ii,:);
    rout = STRFioall.rout(ii,:);
    tmpNL.P = STRFioall.spkNLP(ii,:);
    rpred = TTLNtransfer(bins,tmpNL)/(dt/1000);
    rbg = (tmpNL.P(3)/(dt/1000))*ones(size(rpred));
    NLvbg(ii) = sum(rpred-rbg)/sum(rpred);
    NLerr = rpred-rout;
    NLfit_rms(ii) = std(NLerr);
    NLfit_r2(ii) = 1-(var(NLerr)/var(rpred));
end

save('goodfit.mat','nLL','STRFSTACC','GvTotfrac','GvHfrac','Gstd','NLfit_r2','NLfit_rms','NLvbg')
