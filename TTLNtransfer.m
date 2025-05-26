function [pout dpdspkNL dpdG] = TTLNtransfer(G,spkNL)
% [pout dpdspkNL dpdG] = TTLNtransfer(G,spkNL)
% compute transfer function and gradient
% spkNL is structure with
% spkNL.type  
%    logexp: pout = spkNL.P(3)*log(1+exp(spkNL.P(2)*(G+spkNL.P(1))))
%    logexp_bg: pout = spkNL.P(4)+spkNL.P(3)*log(1+exp(spkNL.P(2)*(G+spkNL.P(1))))
%    exp: pout = spkNL.P(3)*exp(spkNL.P(2)*(G+spkNL.P(1)))
%    exp_bg: pout = spkNL.P(4)+spkNL.P(3)*exp(spkNL.P(2)*(G+spkNL.P(1)))
%    linthresh: pout = spkNL.P(2)*max(G+spkNL.P(1),0)

inp = spkNL.P(2)*(G+spkNL.P(1));
max_inp = 50; %to prevent numerical overflow
too_large = (inp > max_inp);
pout = zeros(size(G));
switch spkNL.type
    case 'logexp'
        pout(~too_large) = spkNL.P(3)*log(1+exp(inp(~too_large)));
        pout(too_large) = spkNL.P(3)*inp(too_large);
    case 'logexp_bg'
        pout(~too_large) = spkNL.P(4)+spkNL.P(3)*log(1+exp(inp(~too_large)));
        pout(too_large) = spkNL.P(4)+spkNL.P(3)*inp(too_large);
    case 'exp'
        pout(~too_large) = spkNL.P(3)*exp(inp(~too_large));
        pout(too_large) = spkNL.P(3)*exp(max_inp);
    case 'exp_bg'
        pout(~too_large) = spkNL.P(4)+spkNL.P(3)*exp(inp(~too_large));
        pout(too_large) = spkNL.P(4)+spkNL.P(3)*exp(max_inp);
    case 'linthresh'
        pout = spkNL.P(2)*max(G+spkNL.P(1),0);
    case 'linthresh_bg'
        pout = spkNL.P(3)+spkNL.P(2)*max(G+spkNL.P(1),0);
    otherwise
        error('Unknown output type in TTtransfer');
        pout = [];
        dpdspkNL = [];
        dpdG = [];
        return
end

if nargout>1 % calculate derivatives needed for gradient
    switch spkNL.type
        case 'logexp'
            dlogexp = zeros(size(G));
            dlogexp(~too_large) = spkNL.P(3)*exp(inp(~too_large))./(1+exp(inp(~too_large)));
            dlogexp(too_large) = spkNL.P(3);
            dpdspkNL = [dlogexp*spkNL.P(2) dlogexp.*(G+spkNL.P(1)) pout/spkNL.P(3)];
            dpdG = dlogexp*spkNL.P(2);
        case 'logexp_bg'
            dlogexp = zeros(size(G));
            dlogexp(~too_large) = spkNL.P(3)*exp(inp(~too_large))./(1+exp(inp(~too_large)));
            dlogexp(too_large) = spkNL.P(3);
            dpdspkNL = [dlogexp*spkNL.P(2) dlogexp.*(G+spkNL.P(1)) (pout-spkNL.P(4))/spkNL.P(3) ones(size(pout))];
            dpdG = dlogexp*spkNL.P(2);
        case 'exp'
            dpdspkNL = [pout*spkNL.P(2) pout.*(G+spkNL.P(1)) pout/spkNL.P(3)]; 
            dpdG = pout*spkNL.P(2);
        case 'exp_bg'
            dpdspkNL = [pout*spkNL.P(2) pout.*(G+spkNL.P(1)) (pout-spkNL.P(4))/spkNL.P(3) ones(size(pout))]; 
            dpdG = (pout-spkNL.P(4))*spkNL.P(2);
        case 'linthresh'
            posvals = G+spkNL.P(1)>0;
            dpdspkNL = [spkNL.P(2)*posvals pout/spkNL.P(2)];
            dpdG = spkNL.P(2)*posvals;
        case 'linthresh_bg'
            posvals = G+spkNL.P(1)>0;
            dpdspkNL = [spkNL.P(2)*posvals pout/spkNL.P(2) ones(size(pout))];
            dpdG = spkNL.P(2)*posvals;
    end
end

%enforce minimum predicted firing rate to avoid nan LLs
pout =  max(pout,1e-20); %minimum predicted probability
pout =  min(pout,1-1e-5); %maximum predicted probability

