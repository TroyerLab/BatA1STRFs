# BatA1STRFs
Measure Spectro-temporal Receptive Fields based on fitting linear-nonlinear cascade model

The model has two inputs:
•	Input from the stimulus obtained by convolving STRF filter with the stimulus (G)
•	Input from spike history, obtained by taking a spike history filter, pasting it on after each spike, and adding over spikes (H)
These two inputs are summed and the neuron uses a rectified linear model plus baseline to generate the probability of having a spike pout = P3+P2*min(0,G+H+P1), where P1 is a bias term (negative of threshold), P2 is the gain, and P3 is the background probability.

Optimization is based on maximizing the log-likelihood. At each time step the likelihood (probability) calculated of the model matching the data (producing a spike or no spike) and these are multiplied across all bins. It’s more convenient to take the log likelihood so that we’re optimizing the sum of the log of the probability across bins. Since most numerical algorithms are based on minimizing error, the model minimizes the negative log likelihood (nLL). 

Any numerical fitting routine is subject to overfitting, in which the parameterized are adjusted precisely to match the particular data obtained. This data includes ‘noise’ in addition to the underlying relationship that the model is intended to extract. One strategy to prevent overfitting is to ‘regularize’ the parameter values, adding penalty terms that restrict the model’s ability to contort itself to fit particular values. We used second order local smoothing in both the time and frequency directions.


A model in contained in an LNmod structure with fields: filtK, spkhist, spkNL, stimP, regP, and opt
•	filtK – the input filter (STRF) as a column vector. reshape(filtK, nLags, nFreqs) returns the STRF matrix
•	spkhist – structure with the following fields
o	len – overall length of the spike history filter (50*2 msec bins = 100 msec history)
o	ncoef – number of coefficients used to specify the filter
o	basis – basis vectors corresponding to each coefficient
o	coef – coefficient values
o	the spike history filter is calculated as spkhist.basis*spkhist.coef
o	currently the basis functions are a series of localized bumps evenly spaced on a log scale so that the resolution is greatest for shortest lags
•	spkNL – having to do with the spike nonlinearity
o	P – list of parameters
o	type – name of spike nonlinearity used (= ‘linthresh_bg’ for this analysis)
•	stimP – stimulus parameters
o	Fields related to STRF time lags: dt, nLags, Lags (=dt*(0:nLags-1)), 
o	Fields related to frequencies: nFreqs, Freqs (in octaves above baseFreq), baseFreq (freq in Hz = baseFreqs*2^Freqs)
•	regP – regularization parameters
o	lam_d2X,lam_d2T (weighting on smoothing penalty in time and freq (X) directions
o	lam_L1 weighting of L1 smoothing penalty = sum(lam_L1*abs(LNmod.filtK))
•	opt – related to optimization – keeps track of nLL and L2 and L1 penalty terms throughout optimization process 
![image](https://github.com/user-attachments/assets/699a872e-0d09-4f7e-8b24-f73d9ebc8e53)
