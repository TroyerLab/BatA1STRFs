function [regP] = TTLNcreate_L2mats(regP,stimP)
%
% [regP] = TTLNcreate_L2mats(regP,stimP)
%
% Creates a set of matrices A specifying different L2 penalties as ||Ak||^2
% INPUTS:
%     LNmod: Input model to which we're adding regularization
% OUTPUTS:
%     L2mats: struct containing all the needed regularization matrices
% 
% The method of computing sparse differencing matrices used here is adapted from 
% Bryan C. Smith's and Andrew V. Knyazev's function "laplacian", available
% here: http://www.mathworks.com/matlabcentral/fileexchange/27279-laplacian-in-1d-2d-or-3d


%% EXTRACT ALL REG PARAMS


%% CREATE L2 REG MATRICES
L2mats = [];
L2List = {};

% set derivative boundaries
x_bound = strcmp(regP.spatial_boundaries,'free');
t_bound = strcmp(regP.temporal_boundaries,'free');

%check for temporal derivative
if regP.lam_dT > 0
    L2List{end+1} = 'dT';
    et = ones(stimP.nLags,1);
    if t_bound == 1 %if free boundary
        et([end]) = 0;
    end
    D1t = spdiags([et -et], [0 1], stimP.nLags, stimP.nLags)';
    
    Ix = speye(stimP.nFreq);
    L2mats.L2_dT = kron(Ix,D1t);
end

%check for spatial derivative
if regP.lam_dX > 0
    L2List{end+1} = 'dX';
    ex = ones(stimP.nFreq,1);
    if x_bound == 1 %if free boundary
        ex([end]) = 0;
    end
    D1x = spdiags([ex -ex], [0 1], stimP.nFreq, stimP.nFreq)';
    
    It = speye(stimP.nLags);
    L2mats.L2_dX = kron(D1x,It);
end

%check for temporal laplacian
if regP.lam_d2T > 0
    L2List{end+1} = 'd2T';
    et = ones(stimP.nLags,1);
    if t_bound == 1 %if free boundary
        et([1 end]) = 0;
    end
    D1t = spdiags([et -2*et et], [-1 0 1], stimP.nLags, stimP.nLags)';
    
    Ix = speye(stimP.nFreq);
    L2mats.L2_d2T = kron(Ix,D1t);
end

%check for spatial laplacian
if regP.lam_d2X > 0
    L2List{end+1} = 'd2X';
    ex = ones(stimP.nFreq,1);
    if x_bound == 1 %if free boundary
        ex([1 end]) = 0;
    end
    D1x = spdiags([ex -2*ex ex], [-1 0 1], stimP.nFreq, stimP.nFreq)';
    
    It = speye(stimP.nLags);
    L2mats.L2_d2X = kron(D1x,It);
end
    
%check for spatio-temporal laplacian
if regP.lam_d2XT > 0
    L2List{end+1} = 'd2XT';
    et = ones(stimP.nLags,1);
    if t_bound == 1 %if free boundary
        et([1 end]) = 0;
    end
    ex = ones(stimP.nFreq,1);
    if x_bound == 1 %if free boundary
        ex([1 end]) = 0;
    end
    D1t = spdiags([et -2*et et], [-1 0 1], stimP.nLags, stimP.nLags)';
    D1x = spdiags([ex -2*ex ex], [-1 0 1], stimP.nFreq, stimP.nFreq)';
    
    It = speye(stimP.nLags);
    Ix = speye(stimP.nFreq);
    L2mats.L2_d2XT = kron(Ix,D1t) + kron(D1x,It);
end

regP.L2mats = L2mats;
regP.L2List = L2List;