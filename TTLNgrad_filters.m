function [nLL, nLLgrad] = TTLNgrad_allfilt(params, LNmod, spkobs, Xstim,Xspkhist)
%
% [nLL, nLLgrad] = TTLNgrad_allfilt(params, LNmod, spkobs, Xstim,Xspkhist)
%
% nLL and nLLgrad with respect to both the stimulus and spike history filters, and
% bias
%
% params list is [filtK; spkhist.coef; bgrate; bias] 

nLL = 0;
nLLgrad = zeros(length(params),1);

filtlen = length(LNmod.filtK);
LNmod.filtK = params(1:filtlen);
LNmod.spkNL.P(1) = params(end);
LNmod.spkNL.P(end) = params(end-1);
if LNmod.spkhist.ncoef>0
    LNmod.spkhist.coef = params(filtlen + (1:LNmod.spkhist.ncoef));
end

if length(params) ~= (filtlen+LNmod.spkhist.ncoef+2)
    error('paramater length does not match filter specifications.')
end

% contribution from filter and spike hisgtory
G = Xstim*LNmod.filtK;
if LNmod.spkhist.ncoef > 0 
    if isempty(LNmod.spkhist.basis)
        H = Xspkhist*LNmod.spkhist.coef;
    else
        H = Xspkhist*LNmod.spkhist.basis*LNmod.spkhist.coef;
    end
else
    H= 0;
end

%calculate dnLLdG 
[pout dpdparams dpdG] = TTLNtransfer(G+H,LNmod.spkNL);
% nLL
[nLL dnLLdp] = TTLNcalcnLL(LNmod,spkobs,pout);
nLL = sum(nLL);
dnLLdG = dnLLdp.*dpdG;
dnLLdpbg = dnLLdp.*dpdparams(:,end);

% assemble gradient
if LNmod.spkhist.len > 0 
    if isempty(LNmod.spkhist.basis)
        nLLgrad = [dnLLdG'*Xstim dnLLdG'*Xspkhist sum(dnLLdpbg) sum(dnLLdG)]';
    else
        nLLgrad = [dnLLdG'*Xstim dnLLdG'*Xspkhist*LNmod.spkhist.basis sum(dnLLdpbg) sum(dnLLdG)]';
    end
else
    nLLgrad = [dnLLdG'*Xstim sum(dnLLdpbg) sum(dnLLdG)]';
end

% add any L2 penalties
if length(LNmod.regP.L2List)>0
    [pen grad_pen] = TTLNcalcL2pen(LNmod);
    nLL = nLL+pen;
    nLLgrad(1:filtlen) = nLLgrad(1:filtlen)+grad_pen;
end

end
