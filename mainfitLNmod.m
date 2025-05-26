% mainfitLNmod - main data flow for fitting LN cascade models to data 

% initial files and directories needed
%   

% flags
doinitmodel = 1; % initialize model
dofitLNmod = 0; % fit nonlinear model to data, one unit at a time
docompileData = 1; % compile data into summary structures
dogoodfit = 1; % calculate goodness of fit
docalcSTRFfilt = 1; % calculate smoothing filter
doSTRFsmooth = 1; % apply smoothing filter to data

%% initialize model parameters save to file 'initmod.mat'
% stimulus parameters
dt = .002; % width of time bins (sec)
nLags = 26; % number of time bins
nFreq = 81; % number of frequency bins
Freqrange = [0 4]; % frequency range in octaves
baseFreq = 5000; % base frequency (Hz)
% regularization parameters
lam_L1 = .01;
lam_d2X = 10;
lam_d2T = 10;
% optimization parameters

if doinitmodel == 1
    disp('initializing model')
    stimP = TTLNcreate_stimP(dt,nLags,nFreq,Freqrange,baseFreq);
    regP = TTLNcreate_regP(stimP,'lam_L1',lam_L1,'lam_d2X',lam_d2X,'lam_d2T',lam_d2T); % regularization parameters
    initmod = TTLNcreate_model(stimP,regP,[],[],'linthresh_bg');
    % optimization meta parameters
    opts.style = ''; % use default style for different functions
    opts.MaxIter = 2000;
    opts.MaxIterations = 2000;
    opts.MaxFunEvals = 2000;
    opts.verbose = 0;
    opts.TolFun = 1e-4; opts.TolX = 1e-8; % matlab
    opts.optTol = 1e-4; opts.progTol = 1e-8; % Mark Schmidt
    save('initmod.mat','initmod','opts');
else
    load('initmod.mat');
end

% load DMR spectrum as resolution of STRF 
load('DMRspec.mat');

% use nonlinear optimization to fit models to data, also calculate spike
% triggered average (STA) - roughly 24 hour run time

if dofitLNmod == 1
    disp('fiting models to data')
    fitLNmod;
end

% compile data into structures for batch processing - save to 'summary.mat'
if docompileData == 1
    disp('compiling data fits')
    compileData;
end

% calculate goodness of fit - append to 'summary.mat'
if dogoodfit == 1
    disp('calculating goodness of fit parameters')
    calcGoodFit;
end

% smooth STRF 
if docalcSTRFfilt == 1
    calcSTRFfilt;
else
    load('STRFsmfilt.mat','STRFsmfilt');
end
if doSTRFsmooth == 1
    disp('smoothing STRFs')
    load('./summary.mat','STRFall');
    STRFsmall = (STRFsmfilt*STRFall')';
    save('./summary.mat','STRFsmfilt','STRFsmall','-append');
end