function ring=UniformMagGroupsErrors(ring)
% function ring=UniformMagGroupsErrors(ring)
%
% makes T1, R1, DeltaS and Tilt error in a mag group 
% identical to the ones on the first element of the group. 
% 
% Groups are defined in the lattice by the MagNum Field
%
%see also: getMagGroupsFromMagNum

mag_groups=getMagGroupsFromMagNum(ring);

for i=1:length(mag_groups)
    
    % group indexes
    ind=mag_groups{i};
    if length(ind)>1
        
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
        
        if isfield(ring{ind(1)},'Tilt')
            Tilt=ring{ind(1)}.Tilt;
        else
            Tilt=0;
        end
        
        if isfield(ring{ind(1)},'DeltaS')
            DeltaS=ring{ind(1)}.DeltaS;
        else
            DeltaS=0;
        end
        
        
        for isl=1:length(ind)
            
            % set all other slices to this value
            ring{ind(isl)}.T1=T1;
            ring{ind(isl)}.T2=T2;
            ring{ind(isl)}.R1=R1;
            ring{ind(isl)}.R2=R2;
            ring{ind(isl)}.Tilt=Tilt;
            ring{ind(isl)}.DeltaS=DeltaS;
            
        end
    end
    
end


