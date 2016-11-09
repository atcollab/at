function ring=atsetbpmerr(varargin)
%ATSETBPMERR sets the misalignment vectors
%
% RING=ATSETBPMERR(RING, ELEMINDEX, OX,OY, GX,GY, RX,RY, Rot) 
% sets the bpm errors of one element or a group of elements 
% 
% ELEMINDEX contains indexes of elements to be misaligned [1xN]
% OX,OY, are reading offsets [1x1 (all identical) or 1xN]
% GX,GY, are gain errors     [1x1 (all identical) or 1xN]
% RX,RY, are resolution errors  [1x1 (all identical) or 1xN]
% Rot is the rotation of the BPM about the s-axis. [1x1 (all identical) or 1xN]
%
% ATSETBPMERR(ELEMINDEX, OX,OY, GX,GY, RX,RY, Rot) Uses the global variable THERING
%
% See also: ATSETTILT ATSETSHIFT

global THERING
if ~iscell(varargin{1})
    THERING=atsetbpmerr(THERING,varargin{:});
else
    [ring,idx,ox,oy,gx,gy,rx,ry,rot]=deal(varargin{:});
    
    if length(ox) == 1
        ox=ox*ones(size(idx)); %#ok<*NASGU>
    elseif length(ox) ~= length(idx)
        error('AT:length','Vector lengths are incompatible: %i/%i.',length(idx),length(ox));
    end
    if length(oy) == 1
        oy=oy*ones(size(idx));
    elseif length(oy) ~= length(idx)
        error('AT:length','Vector lengths are incompatible: %i/%i.',length(idx),length(oy))
    end

    if length(gx) == 1
        gx=gx*ones(size(idx));
    elseif length(gx) ~= length(idx)
        error('AT:length','Vector lengths are incompatible: %i/%i.',length(idx),length(gx));
    end
    if length(gy) == 1
        gy=gy*ones(size(idx));
    elseif length(gy) ~= length(idx)
        error('AT:length','Vector lengths are incompatible: %i/%i.',length(idx),length(gy))
    end

    if length(rx) == 1
        rx=rx*ones(size(idx));
    elseif length(rx) ~= length(idx)
        error('AT:length','Vector lengths are incompatible: %i/%i.',length(idx),length(rx));
    end
    if length(ry) == 1
        ry=ry*ones(size(idx));
    elseif length(ry) ~= length(idx)
        error('AT:length','Vector lengths are incompatible: %i/%i.',length(idx),length(ry))
    end
    
    if length(rot) == 1
        rot=rot*ones(size(idx));
    elseif length(rot) ~= length(idx)
        error('AT:length','Vector lengths are incompatible: %i/%i.',length(idx),length(rot))
    end
    
        
    ring(idx)=cellfun(@(el,a,b)setfield(el,'Offset',[a,b]),ring(idx),num2cell(ox),num2cell(oy),'un',0);
    ring(idx)=cellfun(@(el,a,b)setfield(el,'Scale',1+[a,b]),ring(idx),num2cell(gx),num2cell(gy),'un',0);
    ring(idx)=cellfun(@(el,a,b)setfield(el,'Reading',[a,b]),ring(idx),num2cell(rx),num2cell(ry),'un',0);
    ring(idx)=cellfun(@(el,a)setfield(el,'Rotation',a),ring(idx),num2cell(rot),'un',0); %#ok<*SFLD>
    
    
end
