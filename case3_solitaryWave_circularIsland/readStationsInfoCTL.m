function [StaInfo, nStations] = readStationsInfoCTL(fn)
f = fopen(fn,'r');
nStations = 0;
while ~feof(f)
    s = fgetl(f);
    nStations = nStations + 1;
    [s1,s] = strtok(s); StaInfo(nStations).lon = str2double(s1);
    [s1,s] = strtok(s); StaInfo(nStations).lat = str2double(s1);
    StaInfo(nStations).name = strtrim(s);
end
fclose(f);
