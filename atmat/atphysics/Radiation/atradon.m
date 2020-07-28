function [ring,radelemIndex,cavitiesIndex,energy]=atradon(ring,varargin)
%ATRADON switches RF and radiation on
%
%  [RING2,RADINDEX,CAVINDEX,ENERGY] = ATRADON(RING,CAVIPASS,BENDPASS,QUADPASS)
%    Changes passmethods to get RF cavity acceleration and radiation
%    damping. ATRADON also sets the "Energy" field on the modified elements,
%    looking for the machine energy in:
%       1) 1st 'RingParam' element
%       2) 1st 'RFCavity' element
%       3) field "E0" of the global variable "GLOBVAL"
%
%  INPUTS
%  1. RING      initial AT structure
%  2. CAVIPASS	pass method for cavities
%                '' makes no change,
%               (default: CavityPass)
%  3. BENDPASS	pass method for bending magnets
%               '' makes no change,
%               'auto' substitutes 'Pass' with 'RadPass' in any method
%               (default: 'auto')
%  4. QUADPASS	pass method for quadrupoles
%               '' makes no change,
%               'auto' substitutes 'Pass' with 'RadPass' in any method
%               (default: '')
%  5. SEXTUPASS pass method for sextupoles
%               '' makes no change,
%               'auto' substitutes 'RadPass' with 'Pass' in any method
%               (default: '')
%  6. OCTUPASS  pass method for octupoles
%               '' makes no change,
%               'auto' substitutes 'RadPass' with 'Pass' in any method
%               (default: '')
%
%  OUPUTS
%  1. RING2     Output ring
%  2. RADINDEX  Indices of elements with radiation
%  3. CAVINDEX  Indices of cavities
%
%  See also ATRADOFF, ATCAVITYON, ATCAVITYOFF


[cavipass,bendpass,quadpass,sextupass,octupass]=getargs(varargin,'CavityPass','auto','','','');

energy=atenergy(ring);

[ring,cavities]=changepass(ring,cavipass,@(rg) atgetcells(rg,'Frequency'), 'Cavity');

[ring,dipoles]=changepass(ring,bendpass,@(rg) atgetcells(ring,'BendingAngle',@(elem,bangle) bangle~=0),'Bend');
 
[ring,quadrupoles]=changepass(ring,quadpass,@(rg) selmpole(rg,2),'Quad');

[ring,sextupoles]=changepass(ring,sextupass,@(rg) selmpole(rg,3),'Sextu');

[ring,octupoles]=changepass(ring,octupass,@(rg) selmpole(rg,4),'Octu');

cavitiesIndex=atgetcells(ring,'PassMethod',@(elem,pass) endsWith(pass,'CavityPass'));
radelemIndex=atgetcells(ring,'PassMethod',@(elem,pass) endsWith(pass,'RadPass'));

if any(cavities)
    atdisplay(1,['Cavities modified at position ' num2str(find(cavities)')]);
end

radnum=sum(dipoles|quadrupoles|sextupoles|octupoles);
if radnum > 0
    atdisplay(1,[num2str(radnum) ' elements switched to include radiation']);
end

    function [ring,mask]=changepass(ring,newpass,selfunc,code)
        if isempty(newpass)
            mask=false(size(ring));
        else
            mask=selfunc(ring);
            if any(mask)
                passlist(1:sum(mask),1)=cellstr(newpass);
                ok=~cellfun(@isempty,passlist);
                mask(mask)=ok;
            end
            if any(mask)
                ring(mask)=cellfun(@newelem,ring(mask),passlist(ok),'UniformOutput',false);
            else
                warning(['AT:atradon:NO' code], ['no ' code ' modified']),
            end
        end
        
        function elem=newelem(elem,newpass)
            % Return the modified element
            if ~endsWith(elem.PassMethod,'RadPass')
                if strcmp(newpass,'auto')
                    newpass=strrep(elem.PassMethod,'Pass','RadPass');
                end
                elem.PassMethod=newpass;
                elem.Energy=energy;
            end
        end
    end

    function mask=selmpole(ring,order)
        % Select multipoles of order "order"
        mask=atgetcells(ring,'PolynomB',@ismpole);
        
        function ok=ismpole(elem,polyb)
            om1=order-1;
            ok = ~(isfield(elem,'BendingAngle') && elem.BendingAngle~=0);
            ok=ok && elem.MaxOrder >= om1 && all(polyb(1:om1) == 0) && polyb(order) ~= 0;
        end
    end

end
