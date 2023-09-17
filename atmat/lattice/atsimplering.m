function ring = atsimplering(energy,circumference,hnum,...
    qx,qy,vrf,alpha,varargin)
%ATSIMPLERING   Creates a "simple ring"
%
% A "simple ring" consists in:
% - A RF cavity
% - a 6x6 linear transfer map
% - a detuning and chromaticity element
% - a simplified quantum diffusion element which provides the equilibrium
% emittance, radiation damping and energy loss
%
%RING=ATSIMPLERING(ENERGY,CIRCUMFERENCE,HARMONIC_NUMBER,HTUNE,VTUNE,...
%                  VRF,ALPHAC)
% ENERGY:           nominal energy [eV]
% CIRCUMFERENCE:    ring circumference [m]
% HARMONIC_NUMBER:  harmonic number
% HTUNE:            horizontal tune
% VTUNE:            vertical tune
% VRF:              RF voltage [V]
% ALPHAC:           momentum compaction factor
%
% RING:             generated simple lattice
%
%RING=ATSIMPLERING(...,'name',name,...)
%       specify the name of the generated lattice. Default: ''
%
%RING=ATSIMPLERING(...,'Particle',particle,...)
%       specify the particle. Default: 'relativistic'
%
%RING=ATSIMPLERING(...,'Qpx',xsi_h,'Qpy',xsi_v,...)
%       specify the chromaticites. Default: 0
%
%RING=ATSIMPLERING(...,'A1',A1,'A2',A2,'A3',A3,...)
%       specify the amplitude detuning coefficients. Default: 0
%
%RING=ATSIMPLERING(...,'betax',betax,'betay',betay,...)
%       specify the beta values at origin [m]. Default: 1
%
%RING=ATSIMPLERING(...,'betax',betax,'betay',betay,...)
%       specify the alpha values at origin. Default: 0
%
%RING=ATSIMPLERING(...,'emitx',emitx,'emity',emity,...)
%       specify the equilibrium emittance. Default: 0, meaning no quantum
%       diffusion
%
%RING=ATSIMPLERING(...,'taux',taux,'tauy',tauy,'tauz',tauz,...)
%       specify the damping times. Default: 0, meaning no damping
%
%RING=ATSIMPLERING(...,'U0',U0,...)
%       specify the emergy loss per turn [eV]. Default: 0

clight=PhysConstant.speed_of_light_in_vacuum.value;

[particle,varargs]=getoption(varargin,'Particle',atparticle('relativistic'));
[qpx,varargs]=getoption(varargs,'Qpx',0.0);
[qpy,varargs]=getoption(varargs,'Qpy',0.0);
[A1,varargs]=getoption(varargs,'A1',0.0);
[A2,varargs]=getoption(varargs,'A2',0.0);
[A3,varargs]=getoption(varargs,'A3',0.0);
[taux,varargs]=getoption(varargs,'taux',0.0);
[tauy,varargs]=getoption(varargs,'tauy',0.0);
[tauz,varargs]=getoption(varargs,'tauz',0.0);
[sigma,varargs]=getoption(varargs,'sigma',[]);
[espread,varargs]=getoption(varargs,'Espread',0.0);
if isempty(sigma)
    [alphax,varargs]=getoption(varargs,'alphax',0.0);
    [betax,varargs]=getoption(varargs,'betax',1.0);
    [alphay,varargs]=getoption(varargs,'alphay',0.0);
    [betay,varargs]=getoption(varargs,'betay',1.0);
    [emitx,varargs]=getoption(varargs,'emitx',0.0);
    [emity,varargs]=getoption(varargs,'emity',0.0);
    sigma=atsigma(betax,alphax,emitx,betay,alphay,emity);
end

gammainv=particle.rest_energy/energy;
eta = alpha - gammainv*gammainv;
[vrf,hnum]=broadcast(vrf(:),hnum(:));
frf = hnum * clight / circumference;

[am,lambda]=amat(sigma(1:4,1:4)*jmat(2));
rots=blkdiag(rotmat(qx), rotmat(qy));
[alphax,betax]=ab(am(1,1),am(2,1),am(2,2));
[alphay,betay]=ab(am(3,3),am(4,3),am(4,4));
emit=real(-1i*lambda);

m66=eye(6);
m66(1:4,1:4)=am*rots/am;
m66(6,5)=eta*circumference;

rf_cavity=cellfun(@makerf,num2cell(vrf),num2cell(frf),num2cell(hnum),'UniformOutput',false);
lin_elem=atM66('Linear',m66,'Length',circumference);
nonlin_elem=atbaselem('NonLinear','DeltaQPass',...
    'Betax',betax,'Betay',betay,...
    'Alphax',alphax,'Alphay',alphay,...
    'Qpx',qpx,'Qpy',qpy,...
    'A1',A1,'A2',A2,'A3',A3);
quantdiff=atSimpleQuantDiff('SQD','Betax',betax,'Betay',betay,...
    'Emitx',emit(1),'Emity',emit(2),'Espread',espread,...
    'Taux',taux,'Tauy',tauy,'Tauz',tauz);
ring=atSetRingProperties([rf_cavity;{lin_elem;nonlin_elem;quantdiff}],...
    'Energy',energy,'Periodicity',1,'Particle',particle,varargs{:});

    function rot=rotmat(q)
        cs=cos(2*pi*q);
        sn=sin(2*pi*q);
        rot=[cs, sn; -sn, cs];
    end

    function [alpha,beta]=ab(a11,a21,a22)
        alpha=-a21*a11;
        beta=a11/a22;
    end

    function [v,h]=broadcast(v,h)
        nv=max(numel(v),numel(h));
        h(1:nv,1)=h;
        v(1:nv,1)=v;
    end

    function cavity=makerf(v,f,h)
        cavity=atrfcavity('rf',0.0,v,f,h,energy,'PassMethod','RFCavityPass');
    end

end