function spkhist = TTLNcreate_Xspkhist(spkobs,nSpikehist)
% spkhist = TTLNcreate_Xspkhist(spkobs,nSpikehist)
% make spike history matrix so contribution from spike history is 
% spkhist*nim.spk_hist.coefs

NT = length(spkobs);

spkhist = sparse(NT+nSpikehist,nSpikehist);
for ii = 1:nSpikehist
    spkhist(ii+(1:NT),ii) = spkobs;
end
spkhist(NT+1:end,:) = [];

