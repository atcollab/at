function [s,distance]=atsurvey2spos(r,xycoord,varargin)
% returns closest lattics s coordinates to xycoord points 
% 
% input:
%   r: AT lattice
%   xycoord: 2xN vector of [x,y] cartesian coordinates
%   'slices', value: number of slices to split r 
%                   (more slices = more precision, longer computation time)
%   
% output:
%   s: 1xN vector of s positions in r, closest to xycoord
%   distance: 1xN vector of distances of s to xycoord
% 
%see also: distance2curve

% parse inputs
p = inputParser;
defaultslices = 10^5;

addRequired(p,'r',@iscell);
addRequired(p,'xycoord',@isnumeric);
addOptional(p,'slices',defaultslices,@isnumeric);

parse(p,r,xycoord,varargin{:});
r = p.Results.r;
mapxy = p.Results.xycoord;
npts= p.Results.slices;

% split lattice
rs=splitlattice(r,npts);
G=atgeometry(rs,1:length(rs)+1);

% lattice cartesian coordinates
rx=[G.x];
ry=[G.y];
curvexy=[rx;ry]';

[xy,distance,~] = distance2curve(curvexy,mapxy,'linear');

indmin=arrayfun(@(x)find(curvexy(:,1)>x,1,'first'),xy(:,1));

s=findspos(rs,indmin);

end


function rsplit=splitlattice(ring0,npts)
elmlength=findspos(ring0,1+length(ring0))/npts;
r2=cellfun(@(a)splitelem(a,elmlength),ring0,'UniformOutput',false);
rsplit=cat(1,r2{:});
end

function newelems=splitelem(elem,elmlength)
if isfield(elem,'Length') && elem.Length > 0
    nslices=ceil(elem.Length/elmlength);
    newelems=atdivelem(elem,ones(1,nslices)./nslices);
else
    newelems={elem};
end
end