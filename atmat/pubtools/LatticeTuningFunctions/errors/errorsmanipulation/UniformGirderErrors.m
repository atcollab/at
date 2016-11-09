function [ring,mag_groups]=UniformGirderErrors(ring,mag_groups)
% function ring=UniformGirderErrors(ring)
%
% makes all error in a section delimited by the markers GS and GE
% identical to the ones on the first element
% of the group. Groups are defined in the lattice by the MagNum Field
%
% BPM are moved with the girder setting the appropriate offset.
% 
%see also: getMagGroupsFromGirderIndex

if nargin<2
mag_groups=getMagGroupsFromGirderIndex(ring);
end

for i=1:length(mag_groups)
    
    % group indexes
    ind=mag_groups{i};
    if length(ind)>1
        % uniform BPM offsets
        % find bpm
        bpmingroup=findcells(ring(ind),'Class','Monitor');
        if ~isempty(bpmingroup)
            % assign first bpm alignment and rotation errors to all other BPM
            offerrx=atgetfieldvalues(ring,ind(bpmingroup(1)),'Offset',{1,1});
            offerry=atgetfieldvalues(ring,ind(bpmingroup(1)),'Offset',{1,2});
            roterr=atgetfieldvalues(ring,ind(bpmingroup(1)),'Rotation',{1,1});
            indALL=ind(bpmingroup);
            ring=atsetfieldvalues(ring,indALL,'Offset'   ,offerrx,1,1);
            ring=atsetfieldvalues(ring,indALL,'Offset'   ,offerry,1,2);
            ring=atsetfieldvalues(ring,indALL,'Rotation' ,roterr,1,1);
        end
        
        % uniform magnet errors
        
        % get first element errors
        if isfield(ring{ind(1)},'T1')
            T1=ring{ind(1)}.T1;
        else
            T1=zeros(6,1);
        end
        
        if isfield(ring{ind(1)},'T2')
            T2=ring{ind(1)}.T2;
        else
            T2=zeros(6,1);
        end
        
        if isfield(ring{ind(1)},'R1')
            R1=ring{ind(1)}.R1;
        else
            R1=eye(6,6);
        end
        
        if isfield(ring{ind(1)},'R2')
            R2=ring{ind(1)}.R2;
        else
            R2=eye(6,6);
        end
        
        if isfield(ring{ind(1)},'DeltaS')
            DeltaS=ring{ind(1)}.DeltaS;
        else
            DeltaS=0;
        end
       
        if isfield(ring{ind(1)},'RotAboutS')
            RotAboutS=ring{ind(1)}.RotAboutS;
        else
            RotAboutS=0;
        end
        
        if isfield(ring{ind(1)},'RotAboutX')
            RotAboutX=ring{ind(1)}.RotAboutX;
        else
            RotAboutX=0;
        end
        
        if isfield(ring{ind(1)},'RotAboutY')
            RotAboutY=ring{ind(1)}.RotAboutY;
        else
            RotAboutY=0;
        end
        
        for isl=1:length(ind)
            
            % set all other slices to this value
            ring{ind(isl)}.T1=T1;
            ring{ind(isl)}.T2=T2;
            ring{ind(isl)}.R1=R1;
            ring{ind(isl)}.R2=R2;
            ring{ind(isl)}.DeltaS=DeltaS;
            ring{ind(isl)}.RotAboutS=RotAboutS;
            ring{ind(isl)}.RotAboutX=RotAboutX;
            ring{ind(isl)}.RotAboutY=RotAboutY;
            % ring{ind(isl)}.Offset=[T1(1) T1(3)];
            % ring{ind(isl)}.Rotation=-Tilt;
        end
    end
    
end

%ring=setBpmOffsetOnDipoleRef(ring);


return
