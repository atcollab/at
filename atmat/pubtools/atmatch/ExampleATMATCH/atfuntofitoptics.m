function [val,mi,ma,w]=atfuntofitoptics(r,dpp,opticname,refpts,min,max,weight)
% [val,mi,ma,w]=atfuntofitoptics(r,dpp,opticname,refpts,min,max,weight)
% 
% will output a column array VAL with the value of the parameters
% given in OPTICNAME at the REFPTS positions in the lattice.
% 
% the output provides also an adequate vector of minima maxima and weight
% of the same size of val 
% 
% opticname={
%            'betx',
%            'bety',
%            'alfx',
%            'alfy',
%            'dispx',
%            'dispy',
%            'disppx',
%            'disppy',
%            'gammax',
%            'gammay',
%            'mux',
%            'muy',
%            'xco',
%            'xpco',
%            'xco',
%            'xpco',
%            'Qx',
%            'Qy',
%            'chromX',
%            'chromY',
%            'Mij', (M44 i,j component)
%            'SPos'
%            }
%
%
%
% refpts={[m n k],[j l u],[m u],...}
% min={2,123.34,0,...}
% max={2,123.34,0,...}
% Weigth={2,123.34,0,...}
%
% evaluates atlinopt at the required positions and returns the value in a
% column vector suitable for atmatch constraint
%
% if min max weigth are given, returns vectors usable in atmatch constraint


% get positions
posarray=cell2mat(refpts); % the output will be long as this vector
positions=unique(posarray);


% atlinopt
%[l,t,c]=atlinopt(r,dpp,positions);
[l,t,c]=twissring(r,dpp,positions,'chrom');

val=[];
mi=[];
ma=[];
w=[];

for i=1:length(opticname)
    
    indx=find(positions==refpts{i});
    
    mi=[mi min{i}*ones(size(indx))]; %#ok<*AGROW>
    ma=[ma max{i}*ones(size(indx))];
    w =[w weight{i}*ones(size(indx))];
    
    switch opticname{i}
        case 'SPos'
            b=cat(1,l.SPos);
            val=[val b(indx,1)];
        case 'betx'
            b=cat(1,l.beta);
            val=[val b(indx,1)];
        case 'bety'
            b=cat(1,l.beta);
            val=[val b(indx,2)];
        case 'xco'
            b=cat(1,l.ClosedOrbit);
            val=[val b(indx,1)];
        case 'yco'
            b=cat(1,l.ClosedOrbit);
            val=[val b(indx,3)];
        case 'xpco'
            b=cat(1,l.ClosedOrbit);
            val=[val b(indx,2)];
        case 'ypco'
            b=cat(1,l.ClosedOrbit);
            val=[val b(indx,4)];
        case 'alfx'
            b=cat(1,l.alpha);
            val=[val b(indx,1)];
        case 'alfy'
            b=cat(1,l.alpha);
            val=[val b(indx,2)];
        case 'gammax'
            b=cat(1,l.gamma);
            val=[val b(indx,1)];
        case 'gammay'
            b=cat(1,l.gamma);
            val=[val b(indx,2)];
        case 'dispx'
            b=cat(2,l.Dispersion);
            val=[val b(1,indx)];
        case 'dispy'
            b=cat(2,l.Dispersion);
            val=[val b(3,indx)];
        case 'disppx'
            b=cat(2,l.Dispersion);
            val=[val b(2,indx)];
        case 'disppy'
            b=cat(2,l.Dispersion);
            val=[val b(4,indx)];
        case 'mux'
            b=cat(1,l.mu);
            val=[val b(indx,1)/(2*pi)];
        case 'muy'
            b=cat(1,l.mu);
            val=[val b(indx,2)/(2*pi)];
        case 'Qx'
            val=[val t(1)];
        case 'Qy'
            val=[val t(2)];
        case 'chromX'
            val=[val c(1)];
        case 'chromY'
            val=[val c(2)];
        case 'M11'
            val=[val l.M44(1,1)];
        case 'M12'
            val=[val l.M44(1,2)];
        case 'M13'
            val=[val l.M44(1,3)];
        case 'M14'
            val=[val l.M44(1,4)];
        case 'M21'
            val=[val l.M44(2,1)];
        case 'M22'
            val=[val l.M44(2,2)];
        case 'M23'
            val=[val l.M44(2,3)];
        case 'M24'
            val=[val l.M44(2,4)];
        case 'M31'
            val=[val l.M44(3,1)];
        case 'M32'
            val=[val l.M44(3,2)];
        case 'M33'
            val=[val l.M44(3,3)];
        case 'M34'
            val=[val l.M44(3,4)];
        case 'M41'
            val=[val l.M44(4,1)];
        case 'M42'
            val=[val l.M44(4,2)];
        case 'M43'
            val=[val l.M44(4,3)];
        case 'M44'
            val=[val l.M44(4,4)];
        otherwise
            error([opticname{i} ' : not an allowed variable'])
    end
    
end


return