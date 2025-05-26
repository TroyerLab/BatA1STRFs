function [nLL dnLLdp] = TTLNcalcnLL(LNmod,spkobs,pout)
% [nLL dnLLdp] = TTcalcLL(LNmod,spkobs,pout)

if strcmp(LNmod.opt.LLtype,'Poiss')
    nLL = -(spkobs.* log(pout) - pout); %up to an overall constant
    dnLLdp = -(spkobs./pout-1);
else % binomial
    nLL = -(spkobs.*log(pout)+(1-spkobs).*log(1-pout));
    dnLLdp = -(spkobs./pout-(1-spkobs)./(1-pout));
end

%% NORMALIZE BY number of spikes
Nspks = sum(spkobs);
nLL = nLL/Nspks;
dnLLdp = dnLLdp/Nspks;
