function [modevalP, G, H] = TTLNmodeval(LNmod, spkobs, Xstim,Xspkhist)
%
% [modevalP, G, H] = TTLNmodeval(LNmod, spkobs, Xstim,Xspkhist)
%
% evaluate the model nLL is the sum of the nLL plus L2pen and L1pen
% if nargin<4 assume there is no spike history
% if size(Xstrim,2) == 1 assume it is equal to G
% if size(Xspkhist,2) == 1 assume it is equal to H


filtLen = length(LNmod.filtK);
nLL = 0;

% calculate filter inputs
if size(Xstim,2) == 1
    G = Xstim;
else
    G = Xstim*LNmod.filtK;
end
% calculate spk history input
if nargin<4 || isempty(Xspkhist)
    H = 0;
elseif size(Xspkhist,2) == 1
    H = Xspkhist;
else
    if LNmod.spkhist.ncoef > 0 
        if isempty(LNmod.spkhist.basis)
            H = Xspkhist*LNmod.spkhist.coef;
        else
            H = Xspkhist*LNmod.spkhist.basis*LNmod.spkhist.coef;
        end
    else
        H= 0;
    end
end


%calculate dnLLdG 
[pout dpdparams dpdG] = TTLNtransfer(G+H,LNmod.spkNL);
% nLL
[nLL dnLLdp] = TTLNcalcnLL(LNmod,spkobs,pout);
nLL = sum(nLL);
L2pen = TTLNcalcL2pen(LNmod);
L1pen  = LNmod.regP.lam_L1*sum(abs(LNmod.filtK));

modevalP.nLL = nLL;
modevalP.L2pen = L2pen;
modevalP.L1pen = L1pen;
modevalP.Gstd = std(G);
modevalP.Hstd = std(H);
modevalP.filtKmn = mean(abs(LNmod.filtK)); 
modevalP.spkhistmn = mean(abs(LNmod.spkhist.basis*LNmod.spkhist.coef)); 
modevalP.bgrate = LNmod.spkNL.P(end);

end
