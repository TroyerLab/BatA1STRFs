function STRFinfo = readSTRFinfo(filename)
% STRFinfo = readSTRFinfo(filename)
% read STRFR metadata and return a STRFinfo structure with fields
% STRFs.bat, STRFs.site, STRFs.chn, STRFs.unit 
% data for each field is an Nx1 length numerica vaector
%
% the filename stem (string) for a unit with index = idx can be found with
% filestr = createFilename(STRFinfo,idx);


[num,txt] = xlsread(filename);
N = length(txt)-1;
STRFinfo.bat = zeros(N,1);
STRFinfo.site = zeros(N,1);
STRFinfo.chn = zeros(N,1);
STRFinfo.unit = zeros(N,1);

for ii = 1:N
    if txt{ii+1}(end) ~= ' '
        txt{ii+1}(end+1) = ' ';
    end
    STRFinfo.bat(ii) = str2num(txt{ii+1}(3:5));
    sp = find(txt{ii+1}==' ');
    STRFinfo.site(ii) = str2num(txt{ii+1}(sp(1)+1:sp(2)-1));
    STRFinfo.chn(ii) = str2num(txt{ii+1}(sp(2)+1:sp(3)-1));
    STRFinfo.unit(ii) = str2num(txt{ii+1}(sp(3)+1:sp(4)-1));
end

STRFinfo.screen = num;