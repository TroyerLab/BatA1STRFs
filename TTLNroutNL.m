function [rout bin_centers bin_edges Gn Gx] = TTLNroutNL(LNmod,spkobs,G)
% [rout bin_centers bin_edges Gn Gx] = TTLNroutNL(LNmod,spkobs,G) 
% calculate rout nonlinearity fit and data for single unit model

Gdist_bins = 200;
[Gn,Gx] = hist(G,Gdist_bins);
rout_bins = 100;
% bin_edges = my_prctile(G,linspace(0.05,99.95,rout_bins+1));
bin_edges = my_prctile(G,1:99);
rout_bins = 99-1;
bin_centers = 0.5*bin_edges(1:end-1) + 0.5*bin_edges(2:end);
rout = nan(rout_bins,1);
for i = 1:rout_bins
    cur_set = find(G >= bin_edges(i) & G < bin_edges(i+1));
    if ~isempty(cur_set)
        rout(i) = mean(spkobs(cur_set)==1)/LNmod.stimP.dt;
    end
end
pred = TTLNtransfer(bin_centers,LNmod.spkNL);
rpred = pred/LNmod.stimP.dt;
% pred = LNmod.NLparams(3)*log(1 + exp(LNmod.NLparams(2)*(bin_centers + LNmod.NLparams(1))));
