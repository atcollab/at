function varargout=atwritejson(ring, varargin)
%ATWRITEJSON Create a JSON file to store an AT lattice
%
%JS=ATWRITEJSON(RING)
%   Return the JSON representation of RING as a character array
%
%ATWRITEJSON(RING, FILENAME)
%   Write the JSON representation of RING to the file FILENAME
%
%ATWRITEJSON(RING, ..., 'compact', true)
%   If compact is true, write a compact JSON file (no linefeeds)
%
%see also atloadlattice

[compact, varargs]=getoption(varargin, 'compact', false);
[filename, ~]=getargs(varargs,[]);

if ~isempty(filename)
    %get filename information
    [pname,fname,ext]=fileparts(filename);

    %Check file extension
    if isempty(ext), ext='.json'; end

    % Open file to be written
    [fid,mess]=fopen(fullfile(pname,[fname ext]),'wt');

    if fid==-1
        error('AT:FileErr','Cannot Create file %s\n%s',fn,mess);
    else
        fprintf(fid, sjson(ring));
        fclose(fid);
    end
    varargout={};
else
    varargout={sjson(ring)};
end

    function jsondata=sjson(ring)
        ok=~atgetcells(ring, 'Class', 'RingParam');
        data.atjson= 1;
        data.elements=ring(ok);
        data.properties=get_params(ring);
        jsondata=jsonencode(data, 'PrettyPrint', ~compact);
    end

    function prms=get_params(ring)
        % Get "standard" properties
        [name, energy, part, periodicity, harmonic_number]=...
            atGetRingProperties(ring,'FamName', 'Energy', 'Particle',...
            'Periodicity', 'HarmNumber');
        prms=struct('name', name, 'energy', energy, 'periodicity', periodicity,...
            'particle', saveobj(part), 'harmonic_number', harmonic_number);
        % Add user-defined properties
        idx=atlocateparam(ring);
        if ~isempty(idx)
            flist={'FamName','PassMethod','Length','Class',...
                'Energy', 'Particle','Periodicity','cell_harmnumber'};
            present=isfield(ring{idx}, flist);
            p2=rmfield(ring{idx},flist(present));
            for nm=fieldnames(p2)'
                na=nm{1};
                prms.(na)=p2.(na);
            end
        end
    end

end