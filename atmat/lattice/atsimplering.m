function ring = atsimplering(energy,circumference,hnum,...
    qx,qy,vrf,alpha,varargin)
%ATSIMPLERING   Creates a "simple ring"
%
% A "simple ring" consists of:
% - A RF cavity
% - a 6x6 linear transfer map
% - a simple radiation damping element
% - a detuning and chromaticity element
% - a simplified quantum diffusion element which provides the equilibrium
% emittance
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
%RING=ATSIMPLERING(...,'alphax',alphax,'alphay',alpha,...)
%       specify the alpha values at origin. Default: 0
%
%RING=ATSIMPLERING(...,'dispx',dispx,'dispxp',dispxp,...)
%       specify the horizontal dispersion and derivative at origin. Default: 0
%
%RING=ATSIMPLERING(...,'dispy',dispy,'dispyp',dispyp,...)
%       specify the vertical dispersion and derivative at origin. Default: 0
%
%RING=ATSIMPLERING(...,'emitx',emitx,'emity',emity,...)
%       specify the equilibrium emittance. Default: 0, meaning no quantum
%       diffusion
%
%RING=ATSIMPLERING(...,'taux',taux,'tauy',tauy,'tauz',tauz,...)
%       specify the damping times. Default: 0, meaning no damping
%
%RING=ATSIMPLERING(...,'espread',espread,...)
%       specify the energy spread. Default: 0
%
%RING=ATSIMPLERING(...,'U0',U0,...)
%       specify the emergy loss per turn [eV]. Default: 0

clight=PhysConstant.speed_of_light_in_vacuum.value;

[qpx,varargs]=getoption(varargin,'Qpx',0.0);
[qpy,varargs]=getoption(varargs,'Qpy',0.0);
[A1,varargs]=getoption(varargs,'A1',0.0);
[A2,varargs]=getoption(varargs,'A2',0.0);
[A3,varargs]=getoption(varargs,'A3',0.0);
[taux,varargs]=getoption(varargs,'taux',0.0);
[tauy,varargs]=getoption(varargs,'tauy',0.0);
[tauz,varargs]=getoption(varargs,'tauz',0.0);
[sigma,varargs]=getoption(varargs,'sigma',[]);
[espread,varargs]=getoption(varargs,'espread',0.0);
[dispx, varargs]=getoption(varargs,'dispx',0.0);
[dispxp, varargs]=getoption(varargs,'dispxp',0.0);
[dispy, varargs]=getoption(varargs,'dispy',0.0);
[dispyp, varargs]=getoption(varargs,'dispyp',0.0);
[U0,varargs]=getoption(varargs,'U0',0.0);
if isempty(sigma)
    [alphax,varargs]=getoption(varargs,'alphax',0.0);
    [betax,varargs]=getoption(varargs,'betax',1.0);
    [alphay,varargs]=getoption(varargs,'alphay',0.0);
    [betay,varargs]=getoption(varargs,'betay',1.0);
    [emitx,varargs]=getoption(varargs,'emitx',0.0);
    [emity,varargs]=getoption(varargs,'emity',0.0);
    sigma=atsigma(betax,alphax,10.0,betay,alphay,1.0);
    [am,~]=amat(sigma(1:4,1:4)*jmat(2));
else
    [am,lambda]=amat(sigma(1:4,1:4)*jmat(2));
    [alphax,betax]=ab(am(1,1),am(2,1),am(2,2));
    [alphay,betay]=ab(am(3,3),am(4,3),am(4,4));
    emitx=real(-1i*lambda(1));
    emity=real(-1i*lambda(2));
end

[vrf,hnum]=broadcast(vrf(:),hnum(:));
frf = hnum * clight / circumference;

rots=blkdiag(rotmat(qx), rotmat(qy));

dampx = damp(taux);
dampy = damp(tauy);
dampz = damp(tauz);

m66=eye(6);
m66(1:4,1:4)=am*rots/am;
m66(1,5)=(1.0-m66(1,1))*dispx - m66(1,2)*dispxp;
m66(2,5)=-m66(2,1)*dispx + (1.0-m66(2,2))*dispxp;
m66(3,5)=(1.0-m66(3,3))*dispy - m66(3,4)*dispyp;
m66(4,5)=-m66(4,3)*dispy + (1.0-m66(4,4))*dispyp;
m66(6,5)=alpha*circumference;

rf_cavity=cellfun(@makerf,num2cell(vrf),num2cell(frf),num2cell(hnum),'UniformOutput',false);
lin_elem=atM66('Linear',m66,'Length',circumference);
nonlin_elem=atbaselem('NonLinear','DeltaQPass','Class','NonLinear',...
    'Betax',betax,'Betay',betay,...
    'Alphax',alphax,'Alphay',alphay,...
    'Qpx',qpx,'Qpy',qpy,...
    'A1',A1,'A2',A2,'A3',A3);
simple_rad = atbaselem('SR','SimpleRadiationPass','Class','SimpleRadiation',...
    'dispx',dispx,'dispxp',dispxp,'dispy',dispy,'dispyp',dispyp,...
    'damp_mat_diag',[dampx, dampx, dampy, dampy, dampz, dampz], 'U0', U0);
quantdiff=atSimpleQuantDiff('SQD','betax',betax,'betay',betay,...
    'emitx',emitx,'emity',emity,'espread',espread,...
    'taux',taux,'tauy',tauy,'tauz',tauz);
ring=atSetRingProperties([rf_cavity;{lin_elem; nonlin_elem; simple_rad; quantdiff}],...
    'Energy',energy,'Periodicity',1,varargs{:});

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

    function d = damp(tau)
        if tau == 0.0
            d = 1.0;
        else
            d = exp(-1.0/tau);
        end
    end

    function cavity=makerf(v,f,h)
        cavity=atrfcavity('rf',0.0,v,f,h,energy,'PassMethod','RFCavityPass');
    end

end
