function LNmod_out = TTNLfit_filters(LNmod,spkobs,Xstim,Xspkhist,opt_params)
%
% LNmod_out = TTNLfit_filters(LNmod,spkobs,Xstim,Xspkhist,silent,opttype)
%
% Fit parameters of the input and spike hist filters. 

% INPUTS:
%     LNmod: Model structure containing the current estimates of the spiking nonlinearity parameters (spk_alpha, spk_beta, and spk_theta)
%     spkobs: vector of binned spike counts
%     Xstim: stimulus matrix
%     Xspkhist: spikehistory matrix
%     <silent>: (0 or 1) to determine whether to monitor convergence 
%     <hold_const>: vector specifying which of the parameters to hold constant
%
% OUTPUTS:
%     LNmod_out: updated model structure with fit spk NL parameters

if isempty(opt_params.style) % 'MS'- Mark Schmidt; MAT - Matlab; GD - gradient descent
    opt_params.style = 'MS';
    % opttype = 'GD';
end

if opt_params.verbose == 1
    opts.Display = 'iter';
else
    opts.Display = 'final';
end

init_params = [LNmod.filtK; LNmod.spkhist.coef; LNmod.spkNL.P(end);  LNmod.spkNL.P(1)];

LNmod_out = LNmod;

if LNmod.regP.lam_L1>0 % use L1 'Lasso' regularization
    lamL1 = zeros(size(init_params));
    lamL1(1:length(LNmod.filtK)) = LNmod.regP.lam_L1;
    switch opt_params.style
        case 'MS'
            [fit_params, val, flag, output] = TTLN_L1General2_PSSas(@(K) TTLNgrad_filters(K,LNmod,spkobs,Xstim,Xspkhist), init_params,lamL1,opt_params);
        case 'GD'
            if isempty(opts.lrate) opts.lrate = 1e-4; end
            if isempty(opts.momentum) opts.momentum = .9; end
            [fit_params, val, flag, output] = TTLNgraddescL1(@(K) TTLNgrad_filters(K,LNmod,spkobs,Xstim,Xspkhist), init_params,lamL1,opt_params);
        case 'MAT'
            error('Matlab does not have an optimizer with L1 regularization')
            return
    end
else  % no L1 so must have differentiable objective function
    switch opt_params.style 
        case 'MS'
            %if using Mark Schmidt's optimization, some differences in option parameters
            % opts.optTol = 1e-4; opts.progTol = 1e-8; opts.Method = 'lbfgs';
            % opts = addStructFields(opts,opt_params);
            [fit_params, val, flag, output] = minFunc( @(K) TTLNgrad_filters(K,LNmod,spkobs,Xstim,Xspkhist), init_params,opt_params);
        case 'GD'
            % opts.lrate = 1e-4; opts.momentum = .9; opts.progtol = 1e-8;
            [fit_params, val, flag, output] = TTLNgraddesc(@(K) TTLNgrad_filters(K,LNmod,spkobs,Xstim,Xspkhist), init_params,opt_params);
        case 'MAT'
            % opts.HessUpdate = 'bfgs'; opts.GradObj = 'on';      
            % opts = addStructFields(opts,opt_params);
            [fit_params, val, flag, output] = fminunc(@(K) TTLNgrad_filters(K,LNmod,spkobs,Xstim,Xspkhist), init_params,opt_params);
    end
end
% set output model structure
LNmod_out = LNmod;
filtlen = length(LNmod_out.filtK);
LNmod_out.filtK = fit_params(1:filtlen);
LNmod_out.spkhist.coef = fit_params(filtlen+(1:LNmod_out.spkhist.ncoef));
LNmod_out.spkNL.P(end) = fit_params(end-1);
LNmod_out.spkNL.P(1) = fit_params(end);

% append optimization history
[modeval] = TTLNmodeval(LNmod_out, spkobs, Xstim,Xspkhist);
LNmod_out.opt.nLL_seq = cat(1,LNmod_out.opt.nLL_seq,modeval.nLL);
LNmod_out.opt.L2pen_seq = cat(1,LNmod_out.opt.L2pen_seq,modeval.L2pen);
LNmod_out.opt.L1pen_seq = cat(1,LNmod_out.opt.L1pen_seq,modeval.L1pen);
LNmod_out.opt.history = cat(1,LNmod_out.opt.history,{'allfilt'});
LNmod_out.opt.style{end+1} = opt_params.style;
LNmod_out.opt.flag = cat(1,LNmod_out.opt.flag,flag);
LNmod_out.opt.output{end+1} = output;

