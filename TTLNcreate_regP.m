function regP = TTLNcreate_regP(stimP,varargin)
%
% regP = TTLNcreate_regP(stimP,'lam_type1',lam1_val,'lam_type2',lam2_val',...)
%
% INPUTS: 
%     Lam_type: String that specifies type of regularization
%     Lam_val: Value of regularization hyperparameters
%     CAN SET ANY NUMBER OF PARAMETERS BY INPUTTING PAIRS OF THESE INPUTS
% OUTPUTS: 
%     regP: struct of regularization parameters

%% INITIALIZE REG_PARAMS WITH DEFAULT VALUES
regP.lam_NLd2 = 0; %second derivative of tent basis coefs
regP.lam_dX = 0; %first spatial deriv
regP.lam_dT = 0; %first temporal deriv
regP.lam_d2XT = 0; %spatiotemporal laplacian
regP.lam_d2X = 0; %2nd spatial deriv
regP.lam_d2T = 0; %2nd temporal deriv
regP.lam_L2 = 0; %L2 on filter coefs
regP.lam_L1 = 0; %L1 on filter coefs
regP.L2List = {};
regP.L2mats = [];

%boundary conditions
regP.spatial_boundaries = 'zero'; %default is assume zeros-boundaries in spatial dims (alt. 'free')
regP.temporal_boundaries = 'zero'; %assume 'zero' boundaries for temporal dims (alt. 'free')

%% PARSE INPUTS AND ADD TO REG_PARAMS
if mod(length(varargin),2) ~= 0
    error('Inputs must be in the form: "lam_name",lam_val"');
end
n_inputs = length(varargin)/2;
for ii = 1:n_inputs
   input_name = varargin{(ii-1)*2+1};
   input_val = varargin{(ii-1)*2 + 2};
   if ~isfield(regP,input_name)
       error(['Invalid regularization type ' input_name]);
   else
      regP = setfield(regP,input_name,input_val); 
   end
end

% add L2mats 
[regP]  = TTLNcreate_L2mats(regP,stimP);

