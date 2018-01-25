function elstr=at2str(elem)
%AT2STR - String representation of an AT element
%
%AT2STR Creates a string such that EVAL(AT2STR(ELEM)) recreates an
%identical element.
%
%  INPUTS
%  1. ELEM - Elem to write
%
%  OUTPUTS
%  1. ELSSTR - String given the AT constructor of the element
%
%  See Also atwritem, atwritepy

atclass=atguessclass(elem, 'UseClass');

switch atclass
    case 'Drift'
        create=@atdrift;
        [options,args]=doptions(elem,create,{'Length'});
    case 'Quadrupole'
        create=@atquadrupole;
        args={elem.Length,elem.PolynomB(2)};
        [options,args]=doptions(elem,create,{'Length','K'},{'PolynomB','PolynomA'},args);
        if ~any(abs(getfield2(elem,'PolynomA',0))>0)
            options=rmfield2(options,{'PolynomA'});
        end
        anom=abs(getfield2(elem,'PolynomB',0))~=0;
        anom(2)=false;
        if ~any(anom)
            options=rmfield2(options,{'PolynomB'});
        end
    case 'Sextupole'
        create=@atsextupole;
        args={elem.Length,elem.PolynomB(3)};
        [options,args]=doptions(elem,create,{'Length'},{'PolynomB','PolynomA'},args);
        if ~any(abs(getfield2(elem,'PolynomA',0))>0)
            options=rmfield2(options,{'PolynomA'});
        end
        anom=abs(getfield2(elem,'PolynomB',0))~=0;
        anom(3)=false;
        if ~any(anom)
            options=rmfield2(options,{'PolynomB'});
        end
    case 'Bend'
        args={elem.Length,elem.BendingAngle,elem.PolynomB(2)};
        [options,args]=doptions(elem,@atsbend,{'Length','BendingAngle','K'},...
            {'EntranceAngle','ExitAngle','PolynomB','PolynomA'},args);
        angle1=getfield2(elem,'EntranceAngle',0);
        angle2=getfield2(elem,'ExitAngle',0);
        if all([angle1 angle2] == 0.5*elem.BendingAngle)
            create=@atrbend;
            options=rmfield2(options,{'EntranceAngle','ExitAngle'});
        else
            create=@atsbend;
            if angle1 == 0
                options=rmfield2(options,{'EntranceAngle'});
            end
            if angle2 == 0
                options=rmfield2(options,{'ExitAngle'});
            end
        end
        if ~any(abs(getfield2(elem,'PolynomA',0))>0)
            options=rmfield2(options,{'PolynomA'});
        end
        anom=abs(getfield2(elem,'PolynomB',0))~=0;
        anom(2)=false;
        if ~any(anom)
            options=rmfield2(options,{'PolynomB'});
        end
    case 'Corrector'
        create=@atcorrector;
        [options,args]=doptions(elem,create,{'Length','KickAngle'});
    case 'Multipole'
        create=@atmultipole;
        [options,args]=doptions(elem,create,{'Length','PolynomA','PolynomB'});
    case 'ThinMultipole'
        create=@atthinmultipole;
        [options,args]=doptions(elem,create,{'PolynomA','PolynomB'});
    case 'RFCavity'
        create=@atrfcavity;
        [options,args]=doptions(elem,create,{'Length','Voltage','Frequency','HarmNumber','Energy'});
    case 'RingParam'
        create=@atringparam;
        [options,args]=doptions(elem,create,{'Energy','Periodicity'});
    case 'Aperture'
        create=@ataperture;
        [options,args]=doptions(elem,create,{'Limits'});
    case 'QuantDiff'
        create=@atQuantDiff;
        [options,args]=doptions(elem,create,{'Lmatp'});
    case 'Matrix66'
        create=@atM66;
        [options,args]=doptions(elem,create,{'M66'});
    case 'MatrixTijkPass'
        create=@atM66Tijk;
        [options,args]=doptions(elem,create,{'M66','Tijk'});
    case 'Monitor'
        create=@atmonitor;
        [options,args]=doptions(elem,create);
    otherwise %'Marker'
        create=@atmarker;
        [options,args]=doptions(elem,create);
        if isfield(options,'Energy')
            options=rmfield(options,'Energy');
        end
end
varg=[args,reshape(cat(2,fieldnames(options),struct2cell(options))',1,[])];
fmt=['%15s(%s' repmat(',%s',1,length(varg)-1) ')'];
%strargs=cellfun(@mat2str,varg,'UniformOutput',false);
strargs=cellfun(@(arg) mat2str(reshape(arg,[],size(arg,2))),varg,'UniformOutput',false);
elstr=sprintf(fmt,func2str(create),strargs{:});

    function s2 = rmfield2(s1,fieldnames)
        ok=cellfun(@(fn) isfield(s1,fn), fieldnames);
        s2=rmfield(s1,fieldnames(ok));
    end

    function result = getfield2(s,varargin)
        try
            result=getfield(s,varargin{1:end-1});
        catch
            result=varargin{end};
        end
    end

    function [opts,args]=doptions(elem,create,argn,dontcheck,args)
        if nargin<4, dontcheck={}; end
        if nargin<3, argn={}; end
        if nargin<5, args=cellfun(@(field) elem.(field),argn,'UniformOutput',false); end
        args=[{elem.FamName},args];
        defelem=create(args{:});                        % Build sample element
        if ~strcmp(elem.PassMethod,defelem.PassMethod)
            args=[args,{elem.PassMethod}];
        end
        argnames=[{'FamName'},argn,{'PassMethod'}];
        defel=rmfield(defelem,[argnames,dontcheck]);	% Keep only optional fields
        fnames=fieldnames(defel)';
        fvalues=struct2cell(defel)';
        ok=isfield(elem,fnames);                        % Check for default values
        ok(ok)=cellfun(@(fname,fval) isequal(elem.(fname),fval),fnames(ok),fvalues(ok));
        opts=rmfield2(elem,[argnames,fnames(ok)]);      % Keep only non-default values
    end
end
