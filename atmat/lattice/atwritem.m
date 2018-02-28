function atwritem(ring,filename)
%ATWRITEM Creates a .m file to store an AT structure
%
%  INPUTS
%  1. ring    - Lattice structure (.mat file)
%  2. filname - output file where to write the lattice as a line by line file using
%  element constructors
%
%  EXAMPLES
%  1. atwritem(ring) % Prints the result in the command window
%  2. atwritem(ring,filename) % Prints the result in a file
%
%  See also at2str

if nargin>=2
    %get filename information
    [pname,fname,ext]=fileparts(filename);
    
    %Check file extension
    if isempty(ext), ext='.m'; end
    
    % Make fullname
    fn=fullfile(pname,[fname ext]);
    
    % Open file to be written
    [fid,mess]=fopen(fullfile(pname,[fname ext]),'wt');
    
    if fid==-1
        error('AT:FileErr','Cannot Create file %s\n%s',fn,mess);
    end
    
    fprintf(fid,'function ring=%s()\n',fname);
else
    fid=1;
end

fprintf(fid,'ring={...\n');
ok=cellfun(@(elem) fprintf(fid,'%s;...\n',at2str(elem)),ring); %#ok<NASGU>
fprintf(fid,'};\n');

if nargin>=2
    %close file
    fprintf(fid,'end\n');
    fclose(fid);
end
end

