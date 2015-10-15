function atwritem(ring,filename)
%ATWRITEM Creates a .m file to store an AT structure
%
%ATWRITEM(RING)
%   Prints the result in the command window
%ATWRITEM(RING,FILENAME)
%   Prints the result in a file

if nargin>=2
    [pname,fname,ext]=fileparts(filename);
    if isempty(ext), ext='.m'; end
    fn=fullfile(pname,[fname ext]);
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
    fprintf(fid,'end\n');
    fclose(fid);
end
end

