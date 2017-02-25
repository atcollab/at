function [ring_output,cavitiesIndex]=atcavityoff(ring_input,varargin)
%ATCAVITYOFF	switches cavities off
%
%  [RING2, CAVITIESINDEX] = ATCAVITYOFF(RING,CAVIPASS)
%    Changes passmethods to turn off radiation
%   damping.
%
%  INPUTS:
%    1. RING      initial AT structure
%    2. CAVIPASS  pass method for cavities (default IdentityPass)
%                 '' makes no change,
%
%  OUPUTS
%    1. RING2          output ring with cavities off
%    2. CAVITIESINDEX  indices of radiative elements and cavities
%
%  See also ATCAVITYON, ATRADON, ATRADOFF

%
%% Written by Laurent S. Nadolski

[cavipass] = parseargs({'IdentityPass','auto',''},varargin);

ring_output = ring_input;

if ~isempty(cavipass)
    cavitiesIndex=atgetcells(ring_output,'Frequency');
    if ~any(cavitiesIndex)
        warning('AT:atradon:NoCavity', 'No cavity found in the structure');
    end
    ring_output(cavitiesIndex)=changepass(ring_output(cavitiesIndex),cavipass);
else
    cavitiesIndex=false(size(ring_output));
end

if any(cavitiesIndex)
    atdisplay(1,['Cavities located at position ' num2str(find(cavitiesIndex)')]);
else
    atdisplay(1,'No cavity');
end

    function newline=changepass(line,newpass)
        if strcmp(newpass,'auto')
            passlist=strrep(atgetfieldvalues(line,'PassMethod'),'RadPass','Pass');
        else
            passlist=repmat({newpass},size(line));
        end
        newline=cellfun(@newelem,line,passlist,'UniformOutput',false);
        
        function elem=newelem(elem,newpass)
            elem.PassMethod=newpass;
            %elem=rmfield(elem,'Energy');
        end
    end

end
