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
%  [...] = ATRADOFF(...[,keyword,value]...)
%   The following keywords trigger the processing of the following elements:
%
%   'allpass'       Defines the default pass method for all elements not
%                   explicitly specified. Replaces the following default
%                   values.
%   'bendpass'      pass method for bending magnets. Default 'auto'
%   'quadpass'      pass method for quadrupoles. Default 'auto'
%   'sextupass'     pass method for sextupoles. Default 'auto'
%   'octupass'      pass method for bending magnets. Default 'auto'
%   'wigglerpass'	pass method for wigglers. Default 'auto'
%   'quantdiffpass' pass method for quantum diffusion. Default 'auto'
%
%   OUPUTS
%   1. RING2     Output ring
%   2. RADINDEX  Indices of elements with radiation
%   3. CAVINDEX  Indices of active cavities
%
%  See also ATRADON, ATCAVITYON, ATCAVITYOFF

[allpass,varargs]=getoption(varargin,'allpass',[]);
[octupass,varargs]=getoption(varargs,'octupass',default_pass('auto'));
[sextupass,varargs]=getoption(varargs,'sextupass',default_pass('auto'));
[quadpass,varargs]=getoption(varargs,'quadpass',default_pass('auto'));
[wigglerpass,varargs]=getoption(varargs,'wigglerpass',default_pass('auto'));
[bendpass,varargs]=getoption(varargs,'bendpass',default_pass('auto'));
[cavipass,varargs]=getoption(varargs,'cavipass',default_pass('auto'));
[quantdiffpass,varargs]=getoption(varargs,'quantdiffpass',default_pass('auto'));
[cavipass,bendpass,quadpass]=getargs(varargs,cavipass,bendpass,quadpass);


[ring,cavities]=changepass(ring,cavipass,@(rg) atgetcells(rg,'Frequency'),...
    autoRFPass(cavipass),'Cavity');

[ring,dipoles]=changepass(ring,bendpass,@(rg) atgetcells(ring,'BendingAngle',@(elem,bangle) bangle~=0),...
    autoMultipolePass(bendpass),'Bend');

[ring,quadrupoles]=changepass(ring,quadpass,@(rg) selmpole(rg,2),...
    autoMultipolePass(quadpass),'Quad');

[ring,sextupoles]=changepass(ring,sextupass,@(rg) selmpole(rg,3),...
    autoMultipolePass(sextupass),'Sextu');

[ring,octupoles]=changepass(ring,octupass,@(rg) selmpole(rg,4),...
    autoMultipolePass(octupass),'Octu');

[ring,wigglers]=changepass(ring,wigglerpass,@(rg) atgetcells(rg,'Class','Wiggler'),...
    autoMultipolePass(wigglerpass),'Wiggler');

[ring,others]=changepass(ring,quantdiffpass,@(rg) atgetcells(rg,'Class','QuantDiff'),...
    autoElemPass(quantdiffpass,'QuantDiffPass'),'QuantDiff');

cavitiesIndex=atgetcells(ring,'PassMethod',@(elem,pass) endsWith(pass,'CavityPass'));
radelemIndex=atgetcells(ring,'PassMethod',@(elem,pass) endsWith(pass,'RadPass'));

if any(cavities)
    atdisplay(1,['Cavities modified at position ' num2str(find(cavities)')]);
end

radnum=sum(dipoles|quadrupoles|sextupoles|octupoles|wigglers|others);
if radnum > 0
    atdisplay(1,[num2str(radnum) ' elements switched to disable longitudinal motion']);
end

    function [ring,mask]=changepass(ring,newpass,selfunc,modfunc,code) %#ok<INUSD>
        if isempty(newpass)
            mask=false(size(ring));
        else
            mask=selfunc(ring);
            if any(mask)
                ring(mask)=cellfun(modfunc,ring(mask),'UniformOutput',false);
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

    function modfun=autoMultipolePass(newpass)
        if strcmp(newpass, 'auto')
            modfun=@varelem;
        else
            modfun=@fixelem;
        end
        
        function elem=varelem(elem)
            % 'auto' multipole modification
            elem.PassMethod=strrep(elem.PassMethod,'RadPass','Pass');
        end
        
        function elem=fixelem(elem)
            % Explicit multipole modification
            elem.PassMethod=newpass;
        end
    end

    function modfun=autoElemPass(newpass,defpass)
        if strcmp(newpass, 'auto')
            newpass=defpass;
        end
        modfun=@modelem;
        
        function elem=modelem(elem)
            % Default element modification
            elem.PassMethod=newpass;
        end
    end

    function modfun=autoRFPass(newpass)
        if strcmp(newpass, 'auto')
            modfun=@varelem;
        else
            modfun=@fixelem;
        end
        
        function elem=varelem(elem)
            % 'auto' RF modification
            if elem.Length > 0
                elem.PassMethod='DriftPass';
            else
                elem.PassMethod='IdentityPass';
            end
        end
        
        function elem=fixelem(elem)
            % Explicit RF modification
            elem.PassMethod=newpass;
        end
    end                

    function defpass=default_pass(defpass)
        if ischar(allpass)
            defpass=allpass;
        end
    end
    
end
