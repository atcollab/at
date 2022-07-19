function [ring,radelemIndex,cavitiesIndex]=atradoff(ring,varargin)
%ATRADOFF  switches radiation and cavity off
%
%   [RING2,RADINDEX,CAVINDEX] = ATRADOFF(RING,CAVIPASS,BENDPASS,QUADPASS)
%    Changes passmethods to turn off radiation damping. 
%
%   The default is to turn cavities OFF and remove radiation in all magnets.
%
%  INPUTS
%  1. RING      initial AT structure
%  2. CAVIPASS  pass method for cavities
%               '' makes no change,
%               'auto' sets'IdentityPass' or 'DriftPass' depending of cavity length)
%               (default: 'auto')
%  3. BENDPASS  pass method for bending magnets
%               '' makes no change,
%               'auto' substitutes 'RadPass' with 'Pass' in any method
%               (default: 'auto')
%  4. QUADPASS  pass method for quadrupoles
%               '' makes no change,
%               'auto' substitutes 'RadPass' with 'Pass' in any method
%               (default: 'auto')
%
%  [...] = ATRADOFF(...,keyword,value)
%   The following keywords trigger the processing of the following elements:
%
%   'bendpass'      pass method for bending magnets. Default 'auto'
%   'quadpass'      pass method for quadrupoles. Default 'auto'
%   'sextupass'     pass method for sextupoles. Default 'auto'
%   'octupass'      pass method for bending magnets. Default 'auto'
%   'wigglerpass'	pass method for wigglers. Default 'auto'
%
%   OUPUTS
%   1. RING2     Output ring
%   2. RADINDEX  Indices of elements with radiation
%   3. CAVINDEX  Indices of active cavities
%
%  See also ATRADON, ATCAVITYON, ATCAVITYOFF

[octupass,varargs]=getoption(varargin,'octupass','auto');
[sextupass,varargs]=getoption(varargs,'sextupass','auto');
[quadpass,varargs]=getoption(varargs,'quadpass','auto');
[wigglerpass,varargs]=getoption(varargs,'wigglerpass','auto');
[bendpass,varargs]=getoption(varargs,'bendpass','auto');
[cavipass,varargs]=getoption(varargs,'cavipass','auto');
[cavipass,bendpass,quadpass]=getargs(varargs,cavipass,bendpass,quadpass);

[ring,cavities]=changepass(ring,cavipass,@(rg) atgetcells(rg,'Frequency'),@autoCavityPass,'Cavity');

[ring,dipoles]=changepass(ring,bendpass,@(rg) atgetcells(ring,'BendingAngle',@(elem,bangle) bangle~=0),@autoMultiPolePass,'Bend');

[ring,quadrupoles]=changepass(ring,quadpass,@(rg) selmpole(rg,2),@autoMultiPolePass,'Quad');

[ring,sextupoles]=changepass(ring,sextupass,@(rg) selmpole(rg,3),@autoMultiPolePass,'Sextu');

[ring,octupoles]=changepass(ring,octupass,@(rg) selmpole(rg,4),@autoMultiPolePass,'Octu');

[ring,wigglers]=changepass(ring,wigglerpass,@(rg) atgetcells(rg,'Class','Wiggler'),@autoMultiPolePass,'Wiggler');

cavitiesIndex=atgetcells(ring,'PassMethod',@(elem,pass) endsWith(pass,'CavityPass'));
radelemIndex=atgetcells(ring,'PassMethod',@(elem,pass) endsWith(pass,'RadPass'));

if any(cavities)
    atdisplay(1,['Cavities modified at position ' num2str(find(cavities)')]);
end

radnum=sum(dipoles|quadrupoles|sextupoles|octupoles|wigglers);
if radnum > 0
    atdisplay(1,[num2str(radnum) ' elements switched to include radiation']);
end

    function [ring,mask]=changepass(ring,newpass,selfunc,autopass,code) %#ok<INUSD>
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
            end
        end
        
        function elem=newelem(elem,newpass)
            % Return the modified element
            if strcmp(newpass,'auto')
                newpass=autopass(elem);
            end
            elem.PassMethod=newpass;
            %elem=rmfield(elem,'Energy');
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

    function newpass=autoMultiPolePass(elem)
        % Default PassMethod for multipoles
        newpass=strrep(elem.PassMethod,'RadPass','Pass');
    end

    function newpass=autoCavityPass(elem)
        % Default PassMethod for cavities
        if (elem.Length == 0)
            newpass='IdentityPass';
        else
            newpass='DriftPass';
        end
    end
    
end
