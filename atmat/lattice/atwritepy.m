function pyring=atwritepy(ring,filename)
%ATWRITEPY Creates a .m file to store an AT structure
%
%ATWRITEPY(RING)
%   Prints the result in the command window
%ATWRITEPY(RING,FILENAME)
%   Prints the result in a file

if nargin>=2
    fid=py.open(filename,'wb');
    if fid==-1
        error('AT:FileErr','Cannot Create file %s\n%s',fn,mess);
    end
end

pyring=py.list(reshape(cellfun(@(elem) at2py(elem), ring, 'UniformOutput', false),1,[]));

if nargin>=2
    py.pickle.dump(pyring, fid, 2);
    fid.close();
end
end

