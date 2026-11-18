function pyring=atwritepy(ring,varargin)
%ATWRITEPY Creates pyAT lattice from a Matlab lattice
%
%PYRING=ATWRITEPY(RING)
%   return the python lattice object
%
%PYRING=ATWRITEPY(RING,'file',FILENAME)
%   In addition, store the python lattice object in FILENAME. The ring can
%   be retrieved in python with the commands:
%
%   >>> with open(<filename>,'rb') as fid:
%   ...  ring=pickle.load(fid)
%
%PYRING=ATWRITEPY(...,quiet=true)
%   Do not emit warnings for undefined elements (when Class is not specified
%
%PYRING=ATWRITEPY(...,name='name')
%PYRING=ATWRITEPY(...,energy=ENERGY)
%PYRING=ATWRITEPY(...,energy=PERIODICITY)
%   Set the corresponding attribute on the python lattice. Useful when the
%   lattice has no RingParam element

[filename,varargs]=getoption(varargin,'file','');
if ~isempty(filename)
    fid=py.open(filename,'wb');
end

pyring=py.at.load_var(ring',pyargs(varargs{:}));

if ~isempty(filename)
    py.pickle.dump(pyring, fid, 2);
    fid.close();
end
end
