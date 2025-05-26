function filename = creatFilename(STRFinfo,idx)
% filename = creatFilename(STRFinfo,idx)
% create filename from bat,site,chan,unit info on the idx row of the STRFinfo
% structure

filename = ['Tb' num2str(STRFinfo.bat(idx)) '_Site'];
if STRFinfo.site(idx)<10 filename(end+1) = '0'; end
filename = [filename num2str(STRFinfo.site(idx)) '_Chn'];
if STRFinfo.chn(idx)<10 filename(end+1) = '0'; end
filename = [filename num2str(STRFinfo.chn(idx)) '_Unit'];
filename = [filename num2str(STRFinfo.unit(idx))];
