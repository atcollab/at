function elstr=at2py(elem)
%ELSTR=AT2PY(ELEM) convert AT element tp pyat
%
%AT2PY Creates a pyat element

atclass=atguessclass(elem, 'UseClass');

switch atclass
    case 'Drift'
        create=@atdrift;
        pycreate=@py.at.elements.Drift;
        [args,options]=doptions(elem,create,{'Length'});
    case 'Quadrupole'
        create=@atquadrupole;
        pycreate=@py.at.elements.Quadrupole;
        args={elem.Length,elem.PolynomB(2)};
        [args,options]=doptions(elem,create,{'Length','K'},{'PolynomB','PolynomA'},args);
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
        pycreate=@py.at.elements.Sextupole;
        args={elem.Length,elem.PolynomB(3)};
        [args,options]=doptions(elem,create,{'Length'},{'PolynomB','PolynomA'},args);
        if ~any(abs(getfield2(elem,'PolynomA',0))>0)
            options=rmfield2(options,{'PolynomA'});
        end
        anom=abs(getfield2(elem,'PolynomB',0))~=0;
        anom(3)=false;
        if ~any(anom)
            options=rmfield2(options,{'PolynomB'});
        end
    case 'Bend'
        create=@atsbend;
        pycreate=@py.at.elements.Dipole;
        args={elem.Length,elem.BendingAngle,elem.PolynomB(2)};
        [args,options]=doptions(elem,create,{'Length','BendingAngle','K'},...
            {'EntranceAngle','ExitAngle','PolynomB','PolynomA'},args);
        angle1=getfield2(elem,'EntranceAngle',0);
        angle2=getfield2(elem,'ExitAngle',0);
        if angle1 == 0
            options=rmfield2(options,{'EntranceAngle'});
        end
        if angle2 == 0
            options=rmfield2(options,{'ExitAngle'});
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
        pycreate=@py.at.elements.Drift;
        [args,options]=doptions(elem,create,{'Length','KickAngle'});
    case 'Multipole'
        create=@atmultipole;
        pycreate=@py.at.elements.Multipole;
        [args,options]=doptions(elem,create,{'Length','PolynomA','PolynomB'});
    case 'ThinMultipole'
        create=@atthinmultipole;
        pycreate=@py.at.elements.ThinMultipole;
        [args,options]=doptions(elem,create,{'PolynomA','PolynomB'});
    case 'RFCavity'
        create=@atrfcavity;
        pycreate=@py.at.elements.RFCavity;
        [args,options]=doptions(elem,create,{'Length','Voltage','Frequency','HarmNumber','Energy'});
    case 'RingParam'
        create=@atringparam;
        pycreate=@py.at.elements.RingParam;
        [args,options]=doptions(elem,create,{'Energy','Periodicity'});
    case 'Aperture'
        create=@ataperture;
        pycreate=@py.at.elements.Aperture;
        [args,options]=doptions(elem,create,{'Limits'});
    case 'QuantDiff'
        create=@atQuantDiff;
        pycreate=@py.at.elements.Marker;
        [args,options]=doptions(elem,create,{'Lmatp'});
    case 'Matrix66'
        create=@atM66;
        pycreate=@py.at.elements.Marker;
        [args,options]=doptions(elem,create,{'M66'});
    case 'MatrixTijkPass'
        create=@atM66Tijk;
        pycreate=@py.at.elements.Marker;
        [args,options]=doptions(elem,create,{'M66','Tijk'});
    case 'Monitor'
        create=@atmonitor;
        pycreate=@py.at.elements.Monitor;
        [args,options]=doptions(elem,create);
    otherwise %'Marker'
        create=@atmarker;
        pycreate=@py.at.elements.Marker;
        [args,options]=doptions(elem,create);
        if isfield(options,'Energy')
            options=rmfield(options,'Energy');
        end
end
% varg=[args,reshape(cat(2,fieldnames(options),struct2cell(options))',1,[])];
% fmt=['%15s(%s' repmat(',%s',1,length(varg)-1) ')'];
% %strargs=cellfun(@mat2str,varg,'UniformOutput',false);
% strargs=cellfun(@(arg) mat2str(reshape(arg,[],size(arg,2))),varg,'UniformOutput',false);
% elstr=sprintf(fmt,func2str(create),strargs{:});

args=cellfun(@pyarray,args,'UniformOutput',false);
keys=fieldnames(options);
vals=cellfun(@(k) pyarray(options.(k)), keys, 'UniformOutput',false);
v=reshape([keys vals]',1,[]);
try
    elstr=pycreate(args{:},pyargs(v{:}));
catch
    elstr=py.at.elements.Marker(args{1},pyargs(v{:}));
end

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

    function ret=pyarray(v)
        if isnumeric(v)
            if isscalar(v)
                ret=v;
            elseif isvector(v)
                ret=py.numpy.array(v(:)');
            else
                ret=py.numpy.array(num2cell(v,2));
            end
        else
            ret=v;
        end
    end

    function [args,opts]=doptions(elem,create,argn,dontcheck,args)
        if nargin<4, dontcheck={}; end
        if nargin<3, argn={}; end
        if nargin<5, args=cellfun(@(field) elem.(field),argn,'UniformOutput',false); end
        args=[{elem.FamName},args];
        defelem=create(args{:});                        % Build sample element
        argnames=[{'FamName'},argn];
        if strcmp(elem.PassMethod,defelem.PassMethod)
            argnames=[argnames,{'PassMethod'}];
        end
        defel=rmfield(defelem,[argnames,dontcheck]);	% Keep only optional fields
        fnames=fieldnames(defel)';
        fvalues=struct2cell(defel)';
        ok=isfield(elem,fnames);                        % Check for default values
        ok(ok)=cellfun(@(fname,fval) isequal(elem.(fname),fval),fnames(ok),fvalues(ok));
        opts=rmfield2(elem,[argnames,fnames(ok)]);      % Keep only non-default values
    end
end
