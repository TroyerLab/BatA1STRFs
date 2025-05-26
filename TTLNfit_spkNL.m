function LNmod_out = TTNLfit_spkNL(LNmod,spkobs,GH,opt_params)
%
% LNmod_out = TTNLfit_spkNL(LNmod,spkobs,GH,opt_params)
%
% Fit parameters of the spiking nonlinearity.
% INPUTS:
%     LNmod: Model structure containing the current estimates of the spiking nonlinearity parameters (spk_alpha, spk_beta, and spk_theta)
%     spkobs: vector of binned spike counts
%     GH: total input from filter and spike history
%     <silent>: (0 or 1) to determine whether to monitor convergence 
%     <hold_const>: vector specifying which of the parameters to hold constant
%
% OUTPUTS:
%     LNmod_out: updated model structure with fit spk NL parameters

init_params = LNmod.spkNL.P;

if isempty(opt_params.style)
    opt_params.style = 'MAT';
end

%INITIALIZE CONSTRAINTS 
LB = -1e3 *ones(size(init_params));
UB = 1e3 *ones(size(init_params));
LB(2:4) = 0;

switch opt_params.style
    case 'GD'
        if isempty(opts.lrate) opts.lrate = 1e-4; end
        if isempty(opts.momentum) opts.momentum = .9; end
        [fit_params, val, flag, output] = TTLNgraddescCon(@(K) TTLNgrad_spkNL(K,LNmod,GH,spkobs), init_params,LB,UB,opt_params);
    case 'MAT'
        Aeq = [];Beq = [];
        % opt_params.GradObj = 'off';opt_params.Algorithm = 'active-set';
        % [fit_params, val, flag, output] = fmincon(@(K) TTLNgrad_spkNL(K,LNmod,GH,spkobs), init_params,[],[],Aeq,Beq,LB,UB,[],opt_params);
        opts.Method = 'lbfgs'; opts.HessUpdate = 'bfgs'; opts.GradObj = 'on';
        [fit_params, val, flag, output] = fminunc(@(K) TTLNgrad_spkNL(K,LNmod,GH,spkobs), init_params,opt_params);
end

LNmod_out = LNmod;
LNmod_out.spkNL.P = fit_params;

[modeval] = TTLNmodeval(LNmod_out, spkobs, GH, 0);
LNmod_out.opt.nLL_seq = cat(1,LNmod_out.opt.nLL_seq,modeval.nLL);
LNmod_out.opt.L2pen_seq = cat(1,LNmod_out.opt.L2pen_seq,modeval.L2pen);
LNmod_out.opt.L1pen_seq = cat(1,LNmod_out.opt.L1pen_seq,modeval.L1pen);
LNmod_out.opt.history = cat(1,LNmod_out.opt.history,{'spkNL'});
LNmod_out.opt.style{end+1} = opt_params.style;
LNmod_out.opt.flag = cat(1,LNmod_out.opt.flag,flag);
LNmod_out.opt.output{end+1} = output;


end

