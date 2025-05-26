function [pen grad_pen] = TTLNcalcL2pen(LNmod)
% [pen pengrad] = TTLNcalcpen(LNmod,L2List,L2_mats)
% compute L2 penalties and gradient

%% COMPUTE L2 PENALTIES AND ASSOCIATED CONTRIBUTIONS TO THE LL GRADIENT
pen = 0;
grad_pen = zeros(size(LNmod.filtK));

for ii = 1:length(LNmod.regP.L2List)
    lambda = getfield(LNmod.regP,['lam_' LNmod.regP.L2List{ii}]);
    mat = getfield(LNmod.regP.L2mats,['L2_' LNmod.regP.L2List{ii}]);
    pen = pen + lambda*sum((mat * LNmod.filtK).^2);
    grad_pen = grad_pen+2*lambda*(mat' * mat * LNmod.filtK);
end
