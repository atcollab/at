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

[filename,varargs]=getoption(varargin,'file','');
if ~isempty(filename)
    fid=py.open(filename,'wb');
end

% [keep_all,varargs]=getoption(varargs,'keep_all',false);
% [name,energy,periods,ring_h,particle]=atGetRingProperties(ring,'FamName',...
%     'Energy','Periodicity','HarmNumber','Particle');
% if ~keep_all
%     rp=atgetcells(ring,'Class','RingParam');
%     ring=ring(~rp);
% end
% mlist=cellfun(@(elem) at2py(elem), ring, 'UniformOutput', false);
% pyring=py.at.Lattice(reshape(mlist,1,[]),pyargs('name',name,'energy',energy,...
%     'periodicity',int32(periods),'harmonic_number',ring_h,'particle',particle.name));

pyring=py.at.load_var(ring',pyargs(varargs{:}));

if ~isempty(filename)
    py.pickle.dump(pyring, fid, 2);
    fid.close();
end
end
