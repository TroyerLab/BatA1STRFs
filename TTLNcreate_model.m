function LNmod = TTLNcreate_model(stimP,regP,spkhist,init_filt,spkNL_type)
%
% LNmod = TTLNcreate_model(stim_params,<regP>,<init_filt>,<spkNL_type>)
%
% Creates a new TTLN model
%     above input plus
%     stimP: Struct of stimulus parameters (see TTLNcreate_stimP for details).
%     spkNL: Struct for output nonlinearity
%       spkNL.type: String specifying the type of spiking NL (default is 'logexp_bg')
%       spkNL.P: Vector of parameters for the spiking NL function
%     filtK: Stimulus filter
%     regP: Struct of regularization parameters (see TTcreate_regP)
%     spkhist: Struct containing parameters of the spike history term
%         len: len of spike-history filter 
%         ncoef: number of spike-history coefficients
%         coefs: Coefficients of the spike-history filter
%         basis: matrix of basis functions
%     opt: structure containing the optimization history
%       opt.nLL_seq: Sequence of negative log-likelihood values (recorded after each stage of optimization)
%       opt.L2pen_seq: Same as LL_seq but for the penalized LL 
%       opt.L1pen_seq: Same as LL_seq but for the penalized LL 
%       opt.history: String array that specifies the sequence of parameter optimization stages applied
%       opt.flag: stop criteria flag
%       opt.opt: String array that specifies the sequence of parameter optimization stages applied


if nargin < 2 || isempty(regP)
    regP = TTLNcreate_regP(); %use default regP
end
if nargin < 3 || isempty(spkhist)
    len = ceil(.1/stimP.dt);
    ncoef = 10;
    spkhist = TTLNcreate_spkhist(len,ncoef); % default spkhist
end
if nargin < 4 || isempty(init_filt)
    %use a random vector (this is mostly arbitrary, although initializing to zeros sometimes puts you at a local miLNmodum)
    init_filt = 0.1*randn(stimP.nLags*stimP.nFreq,1)/(stimP.nLags*stimP.nFreq); 
end
if nargin < 5
    spkNL_type = 'logexp_bg'; %default spiking NL
end


filtK = init_filt(:); %store current filter coefs
filtlen = length(filtK); % length of filter

spkNL.type = spkNL_type;
switch spkNL_type
    case 'logexp'
        spkNL.P = [0 1 .01]; 
    case 'logexp_bg'
        spkNL.P = [0 1 .01 .01]; 
    case 'linthresh'
        spkNL.P = [0 .01]; 
    case 'linthresh_bg'
        spkNL.P = [0 .01 .01]; 
    case 'exp'
        spkNL.P = [0 .01]; 
    case 'exp_bg'
        spkNL.P = [0 .01 .01]; 
    otherwise
        error('Unsupported spiking NL type')
end

opt.LLtype = 'Poiss'; % assumption of output distribution for calculating LL
opt.nLL_seq = [];
opt.L2pen_seq = [];
opt.L1pen_seq = [];
opt.history = [];
opt.style = {};
opt.flag = [];
opt.output = {};

% % add L2mats - now done in create_regP
% [L2List L2mats]  = TTLNcreate_L2mats(LNmod);
% LNmod.regP.L2List = L2List;
% LNmod.regP.L2mats = L2mats;

%create model STRUCT
LNmod = struct('opt',opt,'stimP',stimP,'regP',regP,'filtK',filtK,'spkNL',spkNL,'spkhist',spkhist);


