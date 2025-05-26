function [nLL nLLgrad] = TTLNgrad_spkNL(params,LNmod,GH,spkobs)
% [nLL nLLgrad] = TTLNgrad_spkNL(params,LNmod,GH,spkobs)
LNmod.spkNL.P = params;
[pout dpdparams dpdG] = TTLNtransfer(GH,LNmod.spkNL);
[nLL dnLLdp] = TTLNcalcnLL(LNmod,spkobs,pout);
nLL = sum(nLL);
nLLgrad = sum((dnLLdp*ones(1,size(dpdparams,2))).*dpdparams);

