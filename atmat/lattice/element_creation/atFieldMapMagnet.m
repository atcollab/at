function [ elem ] = atFieldMapMagnet( fname, L, Energy, FieldMapFile )
%ATFIELDMAPMAGNET Element creator for FieldMapMagnets
%
%   [ elem ] = atFieldMapMagnet( fname, L, Energy, FieldMapFile )
%
%   FieldMapFile is the name of a matlab file which includes:
%        a vector called X
%        a vector called Y
%        a matrix called LUT (look up table) with
%             size(LUT)=2 x length(X) x length(Y)
%             LUT(1,i,j) is Bx at point X(i),Y(j);        
%             LUT(2,i,j) is By at point X(i),Y(j);        
%
%   see also: atbaselem

load(FieldMapFile);
delta_x=X(2)-X(1);
delta_y=Y(2)-Y(1);
nx=length(X);
ny=length(Y);
LUT_1=zeros([2,size(LUT)]); %#ok<NODEF>
LUT_2=zeros([3,size(LUT)]);

% second order interpolation (Benoit's idea)
LUT_1(1, :, 2:nx-1, :) = (LUT(:, 3:nx, :) - LUT(:, 1:nx-2, :)) / (2*delta_x);
LUT_1(2, :, :, 2:ny-1) = (LUT(:, :, 3:ny) - LUT(:, :, 1:ny-2)) / (2*delta_y);
LUT_2(1, :, 2:nx-1, :)      = 0.5*(LUT(:, 3:nx, :) - 2*LUT(:, 2:nx-1, :) + LUT(:, 1:nx-2, :)) / (delta_x^2);
LUT_2(2, :, 2:nx-1, 2:ny-1) = (LUT(:, 3:nx, 3:ny) - LUT(:, 1:nx-2, 3:ny) - LUT(:, 3:nx, 1:ny-2) + LUT(:, 1:nx-2, 1:ny-2)) / (4*delta_x*delta_y);
LUT_2(3, :, :, 2:ny-1)      = 0.5*(LUT(:, :, 3:ny) - 2*LUT(:, :, 2:ny-1) + LUT(:, :, 1:ny-2)) / (delta_y^2);

cl='FieldMapMagnet';
method='FieldMapPass';
elem=atbaselem(fname,method,'Class',cl,'Length',L,'NumIntSteps',20,...
    'LUT_Bx',squeeze(LUT(1,:,:)),'LUT_By',squeeze(LUT(2,:,:)),...
    'LUT_dBxdx',squeeze(LUT_1(1,1,:,:)),'LUT_dBydx',squeeze(LUT_1(1,2,:,:)),...
    'LUT_dBxdy',squeeze(LUT_1(2,1,:,:)),'LUT_dBydy',squeeze(LUT_1(2,2,:,:)),...
    'LUT_d2Bxdxdx',squeeze(LUT_2(1,1,:,:)),'LUT_d2Bxdxdy',squeeze(LUT_2(2,1,:,:)),'LUT_d2Bxdydy',squeeze(LUT_2(3,1,:,:)),...
    'LUT_d2Bydxdx',squeeze(LUT_2(1,2,:,:)),'LUT_d2Bydxdy',squeeze(LUT_2(2,2,:,:)),'LUT_d2Bydydy',squeeze(LUT_2(3,2,:,:)),...
    'X',X,'Y',Y,'Nx',nx,'Ny',ny,'xmin',min(X),'ymin',min(Y),...
    'delta_x',delta_x,'delta_y',delta_y,'Energy',Energy);

end
