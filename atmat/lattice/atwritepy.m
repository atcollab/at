function pyring=atwritepy(ring,filename)
%ATWRITEPY Creates pyAT lattice from a Matlab lattice
%
%PYRING=ATWRITEPY(RING)
%   return the python lattice object
%
%PYRING=ATWRITEPY(RING,FILENAME)
%   In addition, store the python lattice object in FILENAME. The ring can
%   be retrieved in python with the commands:
%
%   >>> with open(<filename>,'rb') as fid:
%   ...  ring=pickle.load(fid)

if nargin>=2
    fid=py.open(filename,'wb');
end

[energy,periods]=atenergy(ring);
mlist=cellfun(@(elem) at2py(elem), ring, 'UniformOutput', false);
pyring=py.at.Lattice(reshape(mlist,1,[]),pyargs('energy',energy,'periodicity',periods));

if nargin>=2
    py.pickle.dump(pyring, fid, 2);
    fid.close();
end
end
