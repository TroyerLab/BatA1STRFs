% find LNmodel to data 

datadir = './spkdata';
savedir = './STRFdata';

STRFinfo = readClassification('STRFclassification.xlsx');
N = length(STRFinfo.bat);

loadXstim = 1;
loadmodel = 1;

if loadXstim == 1
    'loading stimulus'
    load('DMRspec.mat')
    Xstim = TTLNcreate_time_embedding(DMRspec,stimP); 
    Nhalf = floor(size(Xstim,1)/2);
end


%% make base model
if loadmodel == 1
    load('initmod.mat');
end

% calculate fitted STRFs
for ii = 1:N
    % create filename
    filename = createFilename(STRFinfo,ii);
    if ~exist(fullfile(savedir,[filename '_STRFinfo.mat']))
        % load spikes
        disp([num2str(ii) ': loading spikes'])
        load(fullfile(datadir,[filename  '_info.mat']),'spk');
        spk = spk/1000; % convert to seconds
        spkidx = ceil(spk/stimP.dt);
        % Bin spikes into correct time resolution
        spkobs = zeros(size(Xstim,1),1);
        spkobs(spkidx)=1;
        % make spkhist
        Xspkhist = TTLNcreate_Xspkhist(spkobs,initmod.spkhist.len); 
        % make STA mod and check for 1st/2nd half data consistency
        spkidx1 = spkidx(spkidx<=Nhalf);
        Nspkidx1 = length(spkidx1);
        STA1 =  mean(Xstim(spkidx1,:))';
        spkidx2 = spkidx(spkidx>Nhalf);
        Nspkidx2 = length(spkidx1);
        STA2 =  mean(Xstim(spkidx2,:))';
        STA = (Nspkidx1*STA1+Nspkidx2*STA2)/(Nspkidx1+Nspkidx2);
        tmp = corrcoef(STA1,STA2);
        STACC = tmp(1,2);
        % initial spkNL based on STA
        STAmod = initmod;
        STAmod.filtK = STA;
        G = Xstim*STAmod.filtK;
        STAmod.filtK = STAmod.filtK/std(G);
        G = G/std(G);
        % fit spike nonlinearity based on STA
        disp(['fitting STA spkNL- ' filename]);
        STAmod = TTLNfit_spkNL(STAmod, spkobs, G,opts);
        STANLopt.spkNLP = STAmod.spkNL.P;
        STANLopt.flag = STAmod.opt.flag(end);
        STANLopt.nLL = STAmod.opt.nLL_seq(end);
        STANLopt.output = STAmod.opt.output;
        % fit model
        disp([num2str(ii) ': fitting LNmod - ' filename]);
        fitmod = initmod;
        fitmod.spkNL = STAmod.spkNL; % initialize spkNL to STA fitted
        % scale filtK so that mean size of random filter coefficient match STA
        fitmod.filtK = (mean(abs(STAmod.filtK))/mean(abs(fitmod.filtK)))*fitmod.filtK;
        % fit model
        fitmod = TTLNfit_filters(fitmod,spkobs,Xstim,Xspkhist,opts);
        % store aspects of fit
        STRFopt.spkNLP = fitmod.spkNL.P;
        STRFopt.flag = fitmod.opt.flag(end);
        STRFopt.nLL = fitmod.opt.nLL_seq(end);
        STRFopt.L2pen = fitmod.opt.L2pen_seq(end);
        STRFopt.L1pen = fitmod.opt.L1pen_seq(end);
        STRFopt.output = fitmod.opt.output;
        % refit spkNL
        disp([num2str(ii) 'fitting LNmod spkNL - ' filename])
        STRFG = Xstim*fitmod.filtK;
        STRFH = Xspkhist*fitmod.spkhist.basis*fitmod.spkhist.coef;
        fitmod = TTLNfit_spkNL(fitmod,spkobs,STRFG+STRFH,opts);
        % store relevat values
        STRFNLopt.spkNLP = STAmod.spkNL.P;
        STRFNLopt.flag = fitmod.opt.flag(end);
        STRFNLopt.nLL = fitmod.opt.nLL_seq(end);
        STRFNLopt.output = fitmod.opt.output;
        % find CC between STRF and STA
        STRF = fitmod.filtK;
        tmp = corrcoef(STA,STRF);
        STRFSTACC = tmp(1,2);
        % save input output statistics 
        STRFG = Xstim*fitmod.filtK;
        STRFH = Xspkhist*fitmod.spkhist.basis*fitmod.spkhist.coef;
        % write into saveable data structures
        STRFio.spkNLP = fitmod.spkNL.P;
        STRFio.Gstd = std(STRFG); 
        STRFio.Hstd = std(STRFH);  
        STRFio.GHstd = std(STRFG+STRFH);  
        % find match to data
        [rout bin_centers bin_edges GHn GHx] = TTLNroutNL(fitmod,spkobs,STRFG+STRFH);
        STRFio.rout = rout;
        STRFio.bincent = bin_centers;
        STRFio.GHn = GHn;
        STRFio.GHx = GHx;
        % save data
        save(fullfile(savedir,[filename '_STRFs.mat']),'STRF','STRFSTACC','fitmod','STAmod','STRFio','STRFNLopt','STANLopt','STRFopt');
        datestr(now)
    end
end



