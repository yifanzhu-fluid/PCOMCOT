function dat = COMCOT_readBinaryDataStation(fn)
% read data from binary file to memory
% dat = COMCOT_readBinaryDataStation(fn)

fid = fopen(fn, 'rb');
StartTag = fread(fid, 1, 'integer*4');
ComputeGreen = fread(fid, 1, 'integer*4');
NFaults = fread(fid, 1, 'integer*4');
NDataLength = fread(fid, 1, 'integer*4');
EndTag = fread(fid, 1, 'integer*4');

if(ComputeGreen ~= 1)
    nDatCol = NFaults+1;
else
    nDatCol = 2;
end
dat = zeros(NDataLength, nDatCol);

StartTag = fread(fid, 1, 'integer*4');
t = fread(fid, NDataLength, 'real*8');
dat(:,1) = t(:);
EndTag = fread(fid, 1, 'integer*4');

for iCol = 1:nDatCol-1
    StartTag = fread(fid, 1, 'integer*4');
    h = fread(fid, NDataLength, 'real*8');
    dat(:,1+iCol) = h(:);
    EndTag = fread(fid, 1, 'integer*4');
end

fclose(fid);