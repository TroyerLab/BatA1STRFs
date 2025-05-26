% compile group data from STRFdata directory


datadir = './spkdata';
savedir = './STRFdata';
summaryfilename = 'summary.mat'; % stores summary data over units

initdata = 1;
savesummary = 1;

STRFinfo = readSTRFinfo('STRFinfo.xlsx');
N = length(STRFinfo.bat);

if initdata == 1
    % get model from first file
    filename = createFilename(STRFinfo,1);
    foundSTRF = 0;
    if exist(fullfile(savedir,[filename '_STRFs.mat']))==2
        % get info from individual file
        load(fullfile(savedir,[filename '_STRFs.mat']));
        stimP = fitmod.stimP;
        routN = length(STRFio.rout);
        GhistN = length(STRFio.GHx);
    end
    STAall = zeros(N,stimP.nLags*stimP.nFreq);
    STRFall = zeros(N,stimP.nLags*stimP.nFreq);
    STRFSTACCall = zeros(N,1);
    % nLL
    fitnLLall = zeros(N,1); 
    % input/output statistics
    STAioall.spkNLP = zeros(N,length(fitmod.spkNL.P));
    STRFioall.spkNLP = zeros(N,length(fitmod.spkNL.P));
    STRFioall.Gstd = zeros(N,1); 
    STRFioall.Hstd = zeros(N,1); 
    STRFioall.GHstd = zeros(N,1); 
    STRFioall.filtKmn = zeros(N,1); 
    STRFioall.STRFinfopkhist = zeros(N,fitmod.spkhist.len); 
    STRFioall.STRFinfopkhistmn = zeros(N,1); 
    STRFioall.rout = zeros(N,routN);
    STRFioall.bincent = zeros(N,routN);
    STRFioall.GHx = zeros(N,GhistN);
    STRFioall.GHn = zeros(N,GhistN);
end
    
% compile 
for ii = 1:N
    % disp(num2str(ii));
    % create filename
    filename = createFilename(STRFinfo,ii);
    foundSTRF = 0;
    if exist(fullfile(savedir,[filename '_STRFs.mat']))==2
        % get info from individual file
        load(fullfile(savedir,[filename '_STRFs.mat']),'STRF','STRFSTACC',...
            'STAmod','fitmod','STRFio');
        if exist('STRF')
            foundSTRF = 1;
        else
            clear('STRF');
        end
    end
    if foundSTRF == 1
        % compile raw STRFinfo/STAs
        STRFall(ii,:) = STRF;
        STRFSTACCall(ii) = STRFSTACC;
        STAall(ii,:) = STAmod.filtK;
        % nLL
        fitnLLall(ii) = fitmod.opt.nLL_seq(end);
        % % find G and H statistics        
        STRFioall.spkNLP(ii,:) = STRFio.spkNLP;
        STRFioall.rout(ii,:) = STRFio.rout;
        STRFioall.bincent(ii,:) = STRFio.bincent;
        STRFioall.GHx(ii,:) = STRFio.GHx; % combined G+H
        STRFioall.GHn(ii,:) = STRFio.GHn; % combined G+H
        STRFioall.Gstd(ii) = STRFio.Gstd; 
        STRFioall.Hstd(ii) = STRFio.Hstd; 
        STRFioall.spkhist(ii,:) = fitmod.spkhist.basis*fitmod.spkhist.coef;
        STRFioall.spkhistmn(ii) = mean(abs(STRFioall.spkhist(ii,:))); 
    end
    regP = fitmod.regP;
    stimP = fitmod.stimP;
end

if savesummary == 1
    if exist(summaryfilename)
        save(summaryfilename,'STRFinfo','stimP','regP',...
            'STRFall','STRFSTACCall','STRFioall','STAall','STAioall','fitnLLall','-append');
    else
        save(summaryfilename,'STRFinfo','stimP','regP',...
            'STRFall','STRFSTACCall','STRFioall','STAall','STAioall','fitnLLall');
    end
end



