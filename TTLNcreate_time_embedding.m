function Xmat = TTLNcreate_time_embedding(stim,stimP)
%
% Xmat = TTLNcreate_time_embedding(stim,stimP)
%
% Takes a Txd stimulus matrix and creates a time-embedded matrix of size Tx(d*L),
% where L is the desired number of time lags. If stim is a 3d array 
% the 'spatial dimensions are folded into the 2nd dimension. Assumes zeros-padding.
% Optional up-sampling of stimulus and tent-basis representation for filter
% estimation. Note that Xmatrix is formatted so that adjacent time lags 
% are adjacent within a time-slice of the Xmatrix. Thus X(t,1:nLags) gives
% all time lags of the first spatial pixel at time t.
%
% INPUTS:
%       stim: stimulus matrix (time must be in the first dim.)
%       params: struct of stimulus params (see NIM_create_stim_params)
% OUTPUTS: 
%       Xmat: Time embedded stim matrix

%%
sz = size(stim);

%if there are two spatial dims, fold them into one
if length(sz) > 2
    stim = reshape(stim,sz(1),prod(sz(2:end)));
end

%check that the size of stim matches with the specified stim_params
%structure
[NT,Npix] = size(stim);
if stimP.nFreq ~= Npix
    error('Stimulus dimension mismatch');
end

%for temporal only stimuli 
if Npix == 1 
    Xmat = toeplitz(stim,[stim(1) zeros(1,stimP.dims(1) - 1)]);
else
    %otherwise loop over lags and manually shift the stim matrix
    Xmat = zeros( NT,stimP.nLags*stimP.nFreq);
    for n = 1:stimP.nLags
        Xmat(:,n-1+(1:stimP.nLags:(Npix*stimP.nLags))) = TTLNshift_mat_zpad( stim, n-1, 1);
    end    
end
