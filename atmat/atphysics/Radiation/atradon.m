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
%   The default is to turn cavities ON and set radiation in dipoles and wigglers.
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
%
%  [...] = ATRADON(...,keyword,value)
%   The following keywords trigger the processing of the following elements:
%
%   'bendpass'      pass method for bending magnets. Default 'auto'
%   'quadpass'      pass method for quadrupoles. Default 'auto'
%   'sextupass'     pass method for sextupoles. Default ''
%   'octupass'      pass method for octupoles. Default ''
%   'wigglerpass'   pass method for wigglers. Default 'auto'
%
%  OUPUTS
%  1. RING2     Output ring
%  2. RADINDEX  Indices of elements with radiation
%  3. CAVINDEX  Indices of cavities
%  4. ENERGY	Ring energy
%
%  EXAMPLES
%
%>> ringrad=atradon(ring);
%   Turns cavities on and sets radiation in bending magnets, quadrupoles and wigglers (default)
%
%>> ringrad=atradon(ring,'CavityPass','auto','auto');
%   Turns cavities on and sets radiation in bending magnets, wigglers and quadrupoles
%
%>> ringrad=atradon(ring,'quadpass','');
%   Turns cavities on and sets radiation in bending magnets, wigglers
%
%  See also ATRADOFF, ATCAVITYON, ATCAVITYOFF

[octupass,varargs]=getoption(varargin,'octupass','');
[sextupass,varargs]=getoption(varargs,'sextupass','');
[quadpass,varargs]=getoption(varargs,'quadpass','auto');
[wigglerpass,varargs]=getoption(varargs,'wigglerpass','auto');
[bendpass,varargs]=getoption(varargs,'bendpass','auto');
[cavipass,varargs]=getoption(varargs,'cavipass','CavityPass');
[cavipass,bendpass,quadpass]=getargs(varargs,cavipass,bendpass,quadpass);

energy=atenergy(ring);

[ring,cavities]=changepass(ring,cavipass,@(rg) atgetcells(rg,'Frequency'), 'Cavity');

[ring,dipoles]=changepass(ring,bendpass,@(rg) atgetcells(ring,'BendingAngle',@(elem,bangle) bangle~=0),'Bend');
 
[ring,quadrupoles]=changepass(ring,quadpass,@(rg) selmpole(rg,2),'Quad');

[ring,sextupoles]=changepass(ring,sextupass,@(rg) selmpole(rg,3),'Sextu');

[ring,octupoles]=changepass(ring,octupass,@(rg) selmpole(rg,4),'Octu');

[ring,wigglers]=changepass(ring,wigglerpass,@(rg) atgetcells(rg,'Class','Wiggler'));

cavitiesIndex=atgetcells(ring,'PassMethod',@(elem,pass) endsWith(pass,'CavityPass'));
radelemIndex=atgetcells(ring,'PassMethod',@(elem,pass) endsWith(pass,'RadPass'));

if any(cavities)
    atdisplay(1,['Cavities modified at position ' num2str(find(cavities)')]);
end

radnum=sum(dipoles|quadrupoles|sextupoles|octupoles|wigglers);
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
            elseif nargin >= 4
                warning(['AT:atradon:NO' code], ['no ' code ' modified']),
            end
        end
        
        function elem=newelem(elem,newpass)
            % Return the modified element
            if ~endsWith(elem.PassMethod,'RadPass')
                if strcmp(newpass,'auto')
                    newpass=strrep(elem.PassMethod,'Pass','RadPass');
                end
                % Check the existence of the new PassMethod
                if ~isempty(which(newpass))
                    elem.PassMethod=newpass;
                    elem.Energy=energy;
                end
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
