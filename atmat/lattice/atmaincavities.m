function idcav = atmaincavities(ring)
%ATMAINCAVITIES Get the fundamental mode cavities
%
%IDCAV=ATMAINCAVITIES(RING)
%   Return a refpts pointing to fundamental mode cavities
%   (cavities with the lowest frequency)

idcav=atgetcells(ring,'Frequency');
if any(idcav)
    fall=atgetfieldvalues(ring(idcav),'Frequency');
    [~,~,ic]=unique(fall);
    idcav(ic~=1)=false;
end
idcav=find(idcav);
end