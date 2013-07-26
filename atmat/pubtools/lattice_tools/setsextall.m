function newring=setsextall(ring,sextfile)
%newring=setsextall(ring,fam,val);
newring=ring;

%filenameset=[qemres.datadir '/settings.mat'];

fid=fopen(sextfile,'rt');
if fid < 0, error('setsextall:NoFile','Cannot open file %s',sextfile); end

%use textscan() instead ?
[sfam sval]=textread(sextfile,'%s%f');
%[sfam sval]=textscan(sextfile,'%s%f');
for l=1:length(sfam)
sfam{l};
sval(l);
newring=setsext(newring,sfam{l},sval(l));
end