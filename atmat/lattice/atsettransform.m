function ring = atsettransform(ring,refpts,varargin)
%ATSETTRANSFORM Set the translations and rotations of selected elements
%
%NEWRING=ATSETTRANSFORM(RING,REFPTS,KEY,VALUE,...)
%   Apply ATTRANSFORMELEM to the elements of RING selected by REFPTS.
%   Each VALUE is either a scalar, applied to all the selected elements,
%   or a vector with one value per selected element.
%
%KEYS:
%   'dx','dy','dz','tilt','pitch','yaw','tilt_frame'
%                   See ATTRANSFORMELEM. Default: no change
%   'reference'     'centre' or 'entrance'. See ATTRANSFORMELEM
%
%NEWRING=ATSETTRANSFORM(...,'relative')
%   The input values are added to the previous ones.
%
%See also attransformelem, atsetshift, atsettilt

names={'dx','dy','dz','pitch','yaw','tilt','tilt_frame'};
[relative,args]=getflag(varargin,'relative');
[reference,args]=getoption(args,'reference',[]);
vals=cell(size(names));
for i=1:length(names)
    [vals{i},args]=getoption(args,names{i},[]);
end
if ~isempty(args)
    error('AT:WrongArgument','Unexpected argument "%s"',num2str(args{1}))
end

if islogical(refpts)
    refpts=find(refpts);
end
nb=length(refpts);
for i=1:length(names)
    v=vals{i};
    if isscalar(v)
        vals{i}=v*ones(1,nb);
    elseif ~isempty(v) && length(v) ~= nb
        error('AT:length','Vector lengths are incompatible: %i/%i.',nb,length(v))
    end
end

opts={};
if relative, opts={'relative'}; end
if ~isempty(reference), opts=[opts {'reference',reference}]; end
for k=1:nb
    kv={};
    for i=1:length(names)
        if ~isempty(vals{i})
            kv=[kv {names{i},vals{i}(k)}]; %#ok<AGROW>
        end
    end
    ring{refpts(k)}=attransformelem(ring{refpts(k)},kv{:},opts{:});
end
end
