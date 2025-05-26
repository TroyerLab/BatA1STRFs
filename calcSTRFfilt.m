% make STRF smoothing filter as translation from central template

display = 1; % flags display of filter

cutoffLagmod = 40; % modulation cutoff in time (Hz)
cutoffFreqmod = 4; % modulation cutoff in frequancy (cycles/octave)
nLags = stimP.nLags;
dt = stimP.dt*1000; % time bin in msec
centLagidx = ceil(nLags/2);
nFreq = stimP.nFreq;
df = 4/(nFreq-1); % frequency bin in octaves
centFreqidx = ceil(nFreq/2);

% create half sinusoid centered atidx
lagfilt = zeros(1,nLags);
lagperiod = (1000/cutoffLagmod)/dt; % cut off period in bins
lagsidx = -floor(lagperiod/4):floor(lagperiod/4);
lagfilt(centLagidx+lagsidx) = cos(2*pi*lagsidx/lagperiod);

freqfilt = zeros(nFreq,1);
freqperiod = (1/cutoffFreqmod)/df; % cut off period in octave bins
freqidx = -floor(freqperiod/4):floor(freqperiod/4);
freqfilt(centFreqidx+freqidx) = cos(2*pi*freqidx/freqperiod);

% centfilt = sqrt(freqfilt*lagfilt);
centfilt = (freqfilt*lagfilt);

if display == 1
    clf;
    subplot(2,2,1)
    imagesc(stimP.Lags,stimP.Freqs,centfilt);
            set(gca,'xdir','reverse','ydir','normal','clim',max(abs(centfilt(:)))*[-1 1]);
            set(gca,'ytick',[0:4]);
            ylabel('Octaves')
    xlabel('Lags (msec()')
    
    subplot(2,2,3)
    hold on;
    % plot(stimP.Lags,DMRCCfilt(centFreqidx,:));
    plot(stimP.Lags,lagfilt);
    set(gca,'xdir','reverse','ytick',[]);
    xlabel('Lags (msec()')
    
    subplot(2,2,2)
    hold on;
    plot(freqfilt,stimP.Freqs);
    set(gca,'xtick',[],'yaxisLocation','right')
            ylabel('Octaves')
end

% make convolution matrix
STRFsmfilt = zeros(nFreq*nLags,nFreq*nLags);
% make correlation for each pixel
for jj = 1:nLags
    tmplags = zeros(1,nLags);
    lagoffset = jj-centLagidx;
    if lagoffset>=0
        tmplags(1+lagoffset:nLags) = lagfilt(1:nLags-lagoffset);
    else
        tmplags(1:nLags+lagoffset) = lagfilt(-lagoffset+1:nLags);
    end
    for ii = 1:nFreq
        tmpfreq = zeros(nFreq,1);
        freqoffset = ii-centFreqidx;
        if freqoffset>=0
            tmpfreq(1+freqoffset:nFreq) = freqfilt(1:nFreq-freqoffset);
        else
            tmpfreq(1:nFreq+freqoffset) = freqfilt(-freqoffset+1:nFreq);
        end
        centfilt = (tmpfreq*tmplags)';
        STRFsmfilt(:,nLags*(ii-1)+jj) = centfilt(:);
    end
end

save('STRFsmfilt.mat','freqfilt','lagfilt','STRFsmfilt');
