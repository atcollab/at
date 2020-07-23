function [ring_output,cavitiesIndex,energy]=atcavityon(ring_input,varargin)
%ATCAVITYON switches RF cavities on
%
%  [RING2,CAVINDEX,ENERGY]=ATCAVITYON(RING,CAVIPASS)
%    Changes passmethods to get RF cavity acceleration and radiation
%    damping. ATRADON also sets the "Energy" field on the modified elements,
%    looking for the machine energy in:
%       1) 1st 'RingParam' element
%       2) 1st 'RFCavity' element
%       3) field "E0" of the global variable "GLOBVAL"
%
%  INPUTS
%  1. RING		initial AT structure
%  2. CAVIPASS	pass method for cavities (default CavityPass)
%               '' makes no change
%
%  OUPUTS
%  1. RING2          output ring with cavities off
%  2. CAVITIESINDEX  indices of radiative elements and cavities
%  3. ENERGY         energy
%
%  See also ATRADOFF, ATRADON, ATCAVITYOFF

%
%% Written by Laurent S. Nadolski


[cavipass]=parseargs({'CavityPass'},varargin);

ring_output=ring_input;

energy=atenergy(ring_output);
if ~isempty(cavipass)
    cavitiesIndex=atgetcells(ring_output,'Frequency');
    if any(cavitiesIndex)
        ring_output(cavitiesIndex)=changepass(ring_output(cavitiesIndex),cavipass,energy);
    end
else
    cavitiesIndex=false(size(ring_output));
end

if any(cavitiesIndex)
    atdisplay(1,['Cavities located at position ' num2str(find(cavitiesIndex)')]);
else
    atdisplay(1,'No cavity');
end

    function newline=changepass(line,newpass,nrj)
        if strcmp(newpass,'auto')
            passlist=atgetfieldvalues(line,'PassMethod');
            ok=cellfun(@(psm) isempty(strfind(psm,'RadPass')),passlist);
            passlist(ok)=strrep(passlist(ok),'Pass','RadPass');
        else
            passlist=repmat({newpass},size(line));
        end
        newline=cellfun(@newelem,line,passlist,'UniformOutput',false);
        
        function elem=newelem(elem,newpass)
            elem.PassMethod=newpass;
            elem.Energy=nrj;
        end
    end
end
