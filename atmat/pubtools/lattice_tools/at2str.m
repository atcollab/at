function elstr=at2str(elem)
%ELSTR=AT2STR(ELEM) String representation of an AT element
%
%AT2STR Creates a string such that EVAL(AT2STR(ELEM)) recreates an
%identical element.

args={};
% if isfield(elem,'Class')
%     atclass=elem.Class;
% else
    atclass=atguessclass(elem, 'UseClass');
% end

switch atclass
    case 'Drift'
        elstr=mout('atdrift',elem.FamName,elem.Length);
    case 'Quadrupole'
        if any(abs(elem.PolynomA)>0)
            args=[{'PolynomA',elem.PolynomA} args];
        end
        anom=abs(elem.PolynomB)~=0;
        anom(2)=false;
        if any(anom)
            args=[{'PolynomB',elem.PolynomB} args];
        end
        elstr=mout('atquadrupole',elem.FamName,elem.Length,...
            elem.PolynomB(2),elem.PassMethod,args{:});
    case 'Sextupole'
        if any(abs(elem.PolynomA)>0)
            args=[{'PolynomA',elem.PolynomA} args];
        end
        anom=abs(elem.PolynomB)~=0;
        anom(3)=false;
        if any(anom)
            args=[{'PolynomB',elem.PolynomB} args];
        end
        elstr=mout('atsextupole',elem.FamName,elem.Length,...
            elem.PolynomB(3),args{:});
    case 'Bend'
        if any(abs(elem.PolynomA)>0)
            args=[{'PolynomA',elem.PolynomA} args];
        end
        anom=abs(elem.PolynomB)~=0;
        anom(2)=false;
        if any(anom)
            args=[{'PolynomB',elem.PolynomB} args];
        end
        angles=[elem.EntranceAngle elem.ExitAngle];
        if all(angles == 0)
            code='atsbend';
        elseif all(angles == 0.5*elem.BendingAngle)
            code='atrbend';
        else
            args=[{'EntranceAngle',elem.EntranceAngle,'ExitAngle',elem.ExitAngle} args];
            code='atrbend';
        end
        elstr=mout(code,elem.FamName,elem.Length,elem.BendingAngle,...
            elem.PolynomB(2),elem.PassMethod,args{:});
    case 'Corrector'
        elstr=mout('atcorrector',elem.FamName,elem.Length,elem.KickAngle,args{:});
    case 'Multipole'
        elstr=mout('atmultipole',elem.FamName,elem.Length,...
            elem.PolynomA,elem.PolynomB);
    case 'ThinMultipole'
        elstr=mout('atmultipole',elem.FamName,elem.PolynomA,elem.PolynomB,args{:});
    case 'RFCavity'
        elstr=mout('atrfcavity',elem.FamName,elem.Length,elem.Voltage,...
            elem.Frequency,elem.HarmNumber,elem.Energy);
    case 'RingParam'
        elstr=mout('atringparam',elem.FamName,elem.Energy,elem.Periodicity,args{:});
    otherwise %'Marker'
        fnames=fieldnames(elem);
        fnames(strcmp('FamName',fnames))=[];
        fvals=cellfun(@(fn) elem.(fn),fnames,'UniformOutput',false);
        flds=[fnames';fvals'];
        elstr=mout('atmarker',elem.FamName,flds{:},args{:});
end

    function elstr=mout(fname,elname,varargin)
        fmt=['%15s(''%s''' repmat(',%s',1,length(varargin)) ')'];
        strargs=cellfun(@mat2str,varargin,'UniformOutput',false);
        elstr=sprintf(fmt,fname,elname,strargs{:});
    end
end
