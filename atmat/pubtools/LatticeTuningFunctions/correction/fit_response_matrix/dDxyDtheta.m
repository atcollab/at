function [N,S,Dx,Dy]=dDxyDtheta(r,indbpm,indbend,varargin)
% function [N,S,Dx,Dy]=dDxyDtheta(r,indbpm,indbend,varargin)
%
% Computes the derivative of the dispersion with respect to dipole angles.
% units are [(m/[%])/rad]
% 
% 
% Inputs:
% r         : AT lattice    
% indbpm    : BPM indexes
% indbend   : Bending manget indexes (all the dipoles in r)
% magmodel  : '1slice', 'thick' (default)
%
% Outputs:
% N: horizontal dispersion derivative with respect to bending angles.
% S: vertical dispersion derivative with respect to vertical bending angles.
% 
% Example:
% [N,S,Dx,Dy]=dDxyDtheta(r,....
%                  find(atgetcells(r,'Class','Monitor'))',....
%                  find(atgetcells(r,'Class','Bend')),...
%                  'magmodel','thick');
% 
% https://arxiv.org/pdf/1711.06589.pdf
%see also: atlinopt 

% parse inputs
expectedmodels={'thick','1slice'};

p=inputParser;

addRequired(p,'r');
addRequired(p,'indbpm');
addRequired(p,'indbend');
addParameter(p,'magmodel',expectedmodels{1},...
    @(x) any(validatestring(x,expectedmodels)));
  
parse(p,r,indbpm,indbend,varargin{:});
   
magmodel = p.Results.magmodel; 

switch magmodel
    case '1slice'
        % optics at center
        machx2=cellfun(@(a)atdivelem(a,[0.5,0.5])',r,'un',0);
        machx2=[machx2{:}, atmarker('end')]';
        [lx2,~,~]=atlinopt(machx2,0,1:length(machx2)+1);
        l=lx2(2:2:length(machx2)+1);
    case 'thick'
        % optics at entrance
        [l,~,~]=atlinopt(r,0,1:length(r)+1);
end

% dipoles angle and length
ang_dip = cellfun(@(a)a.BendingAngle,r(indbend),'un',1)';
len_dip = cellfun(@(a)a.Length,r(indbend),'un',1)';

% beta x at bpm and bends
bx_bpm = arrayfun(@(a)a.beta(1),l(indbpm));
bx_dip = arrayfun(@(a)a.beta(1),l(indbend));
by_bpm = arrayfun(@(a)a.beta(2),l(indbpm));
by_dip = arrayfun(@(a)a.beta(2),l(indbend));

% alpha x at bends
ax_dip = arrayfun(@(a)a.alpha(1),l(indbend));
ay_dip = arrayfun(@(a)a.alpha(2),l(indbend));

% phase advance x at bpm and bends
phx_bpm = arrayfun(@(a)a.mu(1),l(indbpm));
phx_dip = arrayfun(@(a)a.mu(1),l(indbend));
phy_bpm = arrayfun(@(a)a.mu(2),l(indbpm));
phy_dip = arrayfun(@(a)a.mu(2),l(indbend));

% tunes including integer part
Q=l(end).mu/2/pi;

% define all quantities as matrices,
[PXdip,PXbpm] = meshgrid(     phx_dip,phx_bpm);
[BXdip,BXbpm] = meshgrid(sqrt(bx_dip),sqrt(bx_bpm)); % square root here for speed
[AXdip,~]     = meshgrid(      ax_dip,bx_bpm);
[Tdip,~]      = meshgrid(     ang_dip,bx_bpm);
[sTdip,~]     = meshgrid(sin(ang_dip),bx_bpm);% sin/cos here for speed
[cTdip,~]     = meshgrid(cos(ang_dip),bx_bpm);
[Ldip,~]      = meshgrid(     len_dip,bx_bpm);

[PYdip,PYbpm] = meshgrid(     phy_dip,phy_bpm);
[BYdip,BYbpm] = meshgrid(sqrt(by_dip),sqrt(by_bpm)); % square root here for speed
[AYdip,~]     = meshgrid(      ay_dip,by_bpm);


% define phase distance
dphX = (PXbpm-PXdip);
dphY = (PYbpm-PYdip);

tau_xmj = abs(dphX) - pi*Q(1);
tau_ymj = abs(dphY) - pi*Q(2);

    function y=signtilde(x)
        y=sign(x)-double(x==0);
    end

% phase term
switch magmodel
    case '1slice'
        
        % hor disp/ dtheta hor
        JCmjx = BXdip.*cos(tau_xmj);
        
        % ver disp/ dtheta ver
        JCmjy = BYdip.*cos(tau_ymj);
        
    case 'thick'
        
        % hor disp/ dtheta hor
        TSm = Ldip .* sTdip ./ (2 * Tdip .* BXdip);
        TCm = BXdip .* cTdip - AXdip .* Ldip .* sTdip ./ (2 * Tdip .* BXdip);
        JCmjx = TCm.*cos(tau_xmj) + TSm.*signtilde(dphX).*sin(tau_xmj);
        
        % ver disp/ dtheta ver
        TSm = Ldip .* sTdip ./ (2 * Tdip .* BYdip);
        TCm = BYdip .* cTdip - AYdip .* Ldip .* sTdip ./ (2 * Tdip .* BYdip);
        JCmjy = TCm.*cos(tau_ymj) + TSm.*signtilde(dphY).*sin(tau_ymj);
end

% dispersion derivative respect to bending angles
N =  BXbpm./(2*sin(pi*Q(1))) .* ( JCmjx );
S = -BYbpm./(2*sin(pi*Q(2))) .* ( JCmjy );

% dispersion
Dx=N*ang_dip';
Dy=S*ang_dip'*0;

end
