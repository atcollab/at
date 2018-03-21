function [z, L] = idtable_global(fname, Nslice, filename, Energy, method)
%IDTABLE Creates a RADIA-Map based element
% idtable(fname, Nslice, filename, Energy, method)
%
% FamName	family name
% Nslice	number of slices (1 means the wiggler is represented by a
%           single kick in the center of the device).
% filename	name of file with wiggler tracking tables.
% Energy    Energy of the machine, needed for scaling
% method    name of the function to use for tracking. Use 'WigTablePass'
%
% The tracking table is described in
% P. Elleaume, "A new approach to the electron beam dynamics in undulators
% and wigglers", EPAC92.
%
% returns assigned address in the FAMLIST that is uniquely identifies
% the family

%---------------------------------------------------------------------------
% Modification Log:
% -----------------
% 13-09-2007:  Created by M. Munoz, based in J. Safranek code.
% 17-11-2008:  Modificated by Z. Marti
%---------------------------------------------------------------------------

ElemData.FamName        = fname;  % add check for identical family names

ElemData.Nslice      	= Nslice;
ElemData.MaxOrder		= 3;
ElemData.NumIntSteps 	= 10;
ElemData.R1             = diag(ones(6,1));
ElemData.R2             = diag(ones(6,1));
ElemData.T1             = zeros(1,6);
ElemData.T2             = zeros(1,6);
ElemData.PassMethod 	= method;


factor=1/((Energy/0.299792458)^2);
factor1=-1/((Energy/0.299792458));
    
% Read the file
D=importdata(filename);
if isfield(D,'Kick1x')
    x=(D.x)';
    y=(D.y)';
    xkick1=factor1*D.Kick1x;
    ykick1=factor1*D.Kick1y;
    xkick=factor*D.Kick2x;
    ykick=factor*D.Kick2y;
    L=D.Len;
    nn=size(xkick1);
    Ny=nn(1);
    Nx=nn(2);
%     ElemData.MultiKick= 1;
%     ElemData.nkicks= nn(3);
else
    [~, data_mat]=mhdrload_bis(filename);
    L=data_mat(:,:,1);
    Nx=data_mat(:,:,2);
    Ny=data_mat(:,:,3);
    A = importdata(filename,' ',10);
    x=A.data;
    x=x(1,1:Nx);
    A = importdata(filename,' ',11);
    txkick=A.data;
    y=txkick(1:Ny,1);
    txkick=txkick(:,2:end);
    A=importdata(filename,' ',11+Ny+3);
    tykick=A.data;
    tykick=tykick(:,2:end);
    A=importdata(filename,' ',11+2*Ny+2*3);
    if isstruct(A)
        txkick1=A.data;
        txkick1=txkick1(:,2:end);
    else
        txkick1=0*txkick;
    end
    A=importdata(filename,' ',11+3*Ny+3*3);
    if isstruct(A)
        tykick1=A.data;
        tykick1=tykick1(:,2:end);
    else
        tykick1=0*tykick;
    end

    xkick=factor*txkick;
    ykick=factor*tykick;
    
    xkick1=factor1*txkick1;
    ykick1=factor1*tykick1;
    
    % Sort arrays in ascending order (needed for "IdTablePass.c")
    [y, indy]=sort(y);
    [x, indx]=sort(x);
    x=x';
    xkick=xkick(indy,indx);
    ykick=ykick(indy,indx);
    
end

ElemData.Length= L;
ElemData.NumX = Ny;
ElemData.NumY = Nx;
ElemData.xtable = x;
ElemData.ytable = y;
ElemData.xkick = xkick;
ElemData.ykick = ykick;
ElemData.xkick1 = xkick1;
ElemData.ykick1 = ykick1;
ElemData.PolynomA= [0 0 0 0];	 
ElemData.PolynomB= [0 0 0 0];

% figure;
% h1=subplot(1,2,1);
% pcolor(x,y,1000*xkick);
% xlabel('x (m)');
% ylabel('y (m)');
% shading interp;
% title ('X Kick Map (mrad)');
% caxis([min([min(min(1000*xkick)) min(min(1000*ykick))]) max([max(max(1000*xkick)) max(max(1000*ykick))])]);
% h2=subplot(1,2,2);
% pcolor(x,y,1000*ykick);
% xlabel('x (m)');
% ylabel('y (m)');
% shading interp;
% title ('Y Kick Map (mrad)');
% caxis([min([min(min(1000*xkick)) min(min(1000*ykick))]) max([max(max(1000*xkick)) max(max(1000*ykick))])]);
% hp=colorbar('SouthOutside');
% p1 = get(h1,'position'); p2 = get(h2,'position'); p0 = get(hp,'position');
% set(h1,'position',[p1(1) p2(2) p1(3) p2(4)]);
% set(hp,'position',[p1(1) p0(2)-0.06 p2(1)-p1(1)+p2(3) p0(4)/1.6]);
% set(h2,'position',p2);
% SetPlotSize(20,10);
% print('-djpeg ','KickMap');

global FAMLIST
z = length(FAMLIST)+1; % number of declare families including this one
FAMLIST{z}.FamName = fname;
FAMLIST{z}.NumKids = 0;
FAMLIST{z}.KidsList= [];
FAMLIST{z}.ElemData= ElemData;
