function stimP = TTLNcreate_stimP(dt,nLags,nFreq,Freqrange,baseFreq)
%
% stimP = TTLNcreate_stimP(dt,nLags,nFreq,Freqrange,baseFreq)
%
% INPUTS: 
%     dt: time bin for STRF (seconds) (usually .002 = 2 msec)
%     nLags: number of time lags (usually 26)
%     nFreq: number of frequency bins (equally spaced on log scale (udually 81)
%     Freqrange: range of frequencies in octaves (usually 0 4)
%     basefreq: base frequency in Hz (usually 5000 Hz)


%%  ADD TO stimP structure
stimP.dt = dt;
stimP.nLags = nLags;
stimP.Lags = (stimP.dt*1000)*(0:stimP.nLags-1)'; % msec
stimP.nFreq = nFreq;
stimP.Freqs = linspace(Freqrange(1),Freqrange(2),stimP.nFreq); % octaves
stimP.baseFreq = baseFreq;

