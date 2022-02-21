function lattice = atloadlattice(fspec,varargin)
%ATLOADLATTICE  load a lattice from a file
%
%LATTICE=ATLOADLATTICE(FILEPATH)
%
%LATTICE=ATLOADLATTICE(FILEPATH,KEY,VALUE...)
%   (Key,value) pairs are added to the RingProperties of the lattice or
%   overload existing ones.
%   Standard keys include: FamName, Energy, Periodicity, HarmNumber, Particle
%   Custom properties are allowed

persistent link_table

if isempty(link_table)
    link_table.mat=@load_mat;
    link_table.m=@load_m;
end

[~,~,fext]=fileparts(fspec);

if isempty(fext), fext='.mat'; end
load_func=link_table.(fext(2:end));
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
        if length(vnames) == 1
            key=vnames{1};
        else
            key='RING';
        end
        [key,opts]=getoption(opts,'matkey',key);
%       if any(cellfun(@(x) strcmp(key,x), vnames))
        try
            lattice=dt.(key);
        catch
            error('AT:load','Cannot find variable %s\nmatkey must be in: %s',...
                key, strjoin(vnames,', '));
        end

    end

end