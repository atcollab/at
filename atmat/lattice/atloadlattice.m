function lattice = atloadlattice(fspec,varargin)
%ATLOADLATTICE  load a lattice from a file
%
%LATTICE=ATLOADLATTICE(FILEPATH)
%   Create an AT lattice (cell array of AT elements) from FILEPATH
%
%LATTICE=ATLOADLATTICE(FILEPATH,'matkey',VARIABLE_NAME,...)
%   Give the name of a variable containing the lattice in .mat files
%   containing several variables. By default, if the .mat file contains a
%   single variable, it will be loaded. Otherwise ATLOADLATTICE will look
%   in order for 'ring', 'lattice' and 'RING'.
%
%LATTICE=ATLOADLATTICE(FILEPATH,KEY,VALUE...)
%   (Key,value) pairs are added to the RingProperties of the lattice or
%   overload existing ones.
%   Standard keys include: FamName, Energy, Periodicity, HarmNumber, Particle
%   Custom properties are allowed
%
%Allowed file types:
%   .mat    Binary Matlab file. If the file contains several variables,
%           the variable name must be specified using the 'matkey' keyword.
%
%   .m      Matlab function. The function must output a valid AT structure.

persistent link_table

if isempty(link_table)
    link_table.mat=@load_mat;
    link_table.m=@load_m;
end

[~,~,fext]=fileparts(fspec);
fext=char(fext);

if isempty(fext), fext='.mat'; end
try
    load_func=link_table.(fext(2:end));
catch
    error('AT:load','AT cannot load %s files', fext)
end
[lattice,varargs]=load_func(fspec,varargin);
lattice=atSetRingProperties(lattice,varargs{:});

    function [lattice,opts]=load_m(fspec,opts)
        [fpath,fname,~]=fileparts(fspec);
        if isempty(fpath)
            lattice=feval(fname);
        else
            here=cd;
            cd(fpath);
            lattice=feval(fname);
            cd(here);
        end
    end

    function [lattice,opts]=load_mat(fpath, opts)
        dt=load(fpath);
        vnames=fieldnames(dt);
        key='RING';
        if length(vnames) == 1
            key=vnames{1};
        else
            for v={'ring','lattice'}
                vn=v{1};
                if isfield(dt,vn)
                    key=vn;
                    break;
                end
            end
        end
        [key,opts]=getoption(opts,'matkey',key);
        try
            lattice=dt.(key);
        catch
            error('AT:load','Cannot find variable %s\nmatkey must be in: %s',...
                key, strjoin(vnames,', '));
        end

    end

end