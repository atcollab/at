function rerr=atsetwaveerrors(rerr,magindex,indBPM,W,A,type)
%function rerr=atsetwaveerrors(...
% rerr, lattice
% magindex, indexes of magnts to set errors
% indBPM,  index of bpm to mode reference using offsets
% W, wave length
% A, amplitude
% type)
%
% Set Error wave
%
% is magindex is a cellarray, then also type, W and A must be. This will be
% applied summing.
%
% if W and A are vectors the resulting sinusoids are summed. 
% 
% type may be: x, y, s, psi, theta, phi, bpm'
%              x.y, x.y.psi, x.y.s.psi,
%              x.y.s.psi.theta.phi'
%
%see also: setANYshift setDSerr setTiltAbout seterrorwave

if iscell(magindex)
    
    if length(magindex)~=length(W) || length(magindex)~=length(A) || length(magindex)~=length(type)
        error('magindex, A, W, type ,must be cell array of the same length, or numeric.')
    end
    
    % loop magindex to apply errors
    for im=1:length(magindex)
        
       disp([...
            'magindex ' num2str(im)...
            ', type:' type{im}...
            ', A:' num2str(A{im}(1)) ' ... ' num2str(A{im}(end))...
            ', W: ' num2str(W{im}(1)) ' ... ' num2str(W{im}(end))]...
            );
        
        rerr=atsetwaveerrors(...
            rerr,...
            magindex{im},...
            indBPM,...
            W{im},...
            A{im},...
            type{im});
    end
    
    
else
    
    % set error wave to zero if no error existing
    [X0,Y0,S0,T0,R0,P0,bpmerr]=GetExistingErrors(rerr,magindex);
    
    if (std(X0)==0 && std(Y0)==0)
        % sets all field of T1 to 0
        errfun=@(r,po,er)setANYshift(r,po,1:6,er);
        rerr=seterrorwave(rerr,magindex,errfun,1,0);
    end
    
    % set error torequired amplitude and wavelength
    switch type
        case 'x'
            
            disp('X error');
            
            errfun=@(r,po,er)setANYshift(r,po,1,er); % sets X errors
            rerr=seterrorwave(rerr,magindex,errfun,W,A,X0);
            
            errfun=@(r,po,er)setcellstruct(r,'Offset',po,er,1,1); % sets X bpm errors
            rerr=seterrorwave(rerr,indBPM,errfun,W,-A,bpmerr.offsetx);
        case 'y'
            
            disp('Y error');
            
            errfun=@(r,po,er)setANYshift(r,po,3,er); % sets Y errors
            rerr=seterrorwave(rerr,magindex,errfun,W,A,Y0);
            
            errfun=@(r,po,er)setcellstruct(r,'Offset',po,er,1,2); % sets Y bpm errors
            rerr=seterrorwave(rerr,indBPM,errfun,W,-A,bpmerr.offsety);
            
        case 's' % longitudinal displacement
            
            disp('S error');
            
            errfun=@(r,po,er)setDSerr(r,po,er); % sets Y errors
            rerr=seterrorwave(rerr,magindex,errfun,W,A,S0);
            
        case 'phi' % about x
            
            disp('phi error');
            
            errfun=@(r,po,er)setTiltAbout(r,po,'x',er); % sets Y errors
            rerr=seterrorwave(rerr,magindex,errfun,W,A,R0);
            
        case 'theta' % about y
            
            disp('theta error');
            
            errfun=@(r,po,er)setTiltAbout(r,po,'y',er); % sets Y errors
            rerr=seterrorwave(rerr,magindex,errfun,W,A,P0);
            
        case 'psi' % about s
            
            disp('s-asix rotation error');
            
            errfun=@(r,po,er)setTiltAbout(r,po,'s',er); % sets Y errors
            rerr=seterrorwave(rerr,magindex,errfun,W,A,T0);
            
            %  errfun=@(r,po,er)setcellstruct(r,'Rotation',po,er,1); % sets Y bpm errors
            %  rerr=seterrorwave(rerr,indBPM,errfun,W,-A);
            
        case 'bpm'
            
            % bpm
            % errfun=@(r,po,er)setcellstruct(r,'Rotation',po,er,1); % sets psi bpm errors
            % rerr=seterrorwave(rerr,indBPM,errfun,W,-A);
            errfun=@(r,po,er)setcellstruct(r,'Offset',po,er,1,1); % sets x bpm errors
            rerr=seterrorwave(rerr,indBPM,errfun,W,-A,bpmerr.offsetx);
            errfun=@(r,po,er)setcellstruct(r,'Offset',po,er,1,2); % sets Y bpm errors
            rerr=seterrorwave(rerr,indBPM,errfun,W,-A,bpmerr.offsety);
            
        case 'x.y.s.psi'
            disp('x,y,s misal, bpm and s-asix rotation error (psi only, no phi and theta)');
            
            % rotation
            errfun=@(r,po,er)settilt_THERING_Dipole(r,po,er); % sets rotation about s
            rerr=seterrorwave(rerr,magindex,errfun,W,A,T0);
            
            % alignment
            errfun=@(r,po,er)setANYshift(r,po,1,er); % sets x errors
            rerr=seterrorwave(rerr,magindex,errfun,W,A,X0);
            
            errfun=@(r,po,er)setANYshift(r,po,3,er); % sets Y errors
            rerr=seterrorwave(rerr,magindex,errfun,W,A,Y0);
            
            errfun=@(r,po,er)setDSerr(r,po,er); % sets S errors
            rerr=seterrorwave(rerr,magindex,errfun,W,A,S0);
            
            % bpm
            errfun=@(r,po,er)setcellstruct(r,'Offset',po,er,1,1); % sets x bpm errors
            rerr=seterrorwave(rerr,indBPM,errfun,W,-A,bpmerr.offsetx);
            errfun=@(r,po,er)setcellstruct(r,'Offset',po,er,1,2); % sets Y bpm errors
            rerr=seterrorwave(rerr,indBPM,errfun,W,-A,bpmerr.offsety);
            
        case 'x.y.psi'
            disp('x,y misal, bpm and s-asix rotation error (psi only, no phi and theta)');
            
            % rotation
            errfun=@(r,po,er)settilt_THERING_Dipole(r,po,er); % sets rotation about s
            rerr=seterrorwave(rerr,magindex,errfun,W,A,T0);
            
            % alignment
            errfun=@(r,po,er)setANYshift(r,po,1,er); % sets x errors
            rerr=seterrorwave(rerr,magindex,errfun,W,A,X0);
            
            errfun=@(r,po,er)setANYshift(r,po,3,er); % sets Y errors
            rerr=seterrorwave(rerr,magindex,errfun,W,A,Y0);
            
            % bpm
            errfun=@(r,po,er)setcellstruct(r,'Offset',po,er,1,1); % sets x bpm errors
            rerr=seterrorwave(rerr,indBPM,errfun,W,-A,bpmerr.offsetx);
            errfun=@(r,po,er)setcellstruct(r,'Offset',po,er,1,2); % sets Y bpm errors
            rerr=seterrorwave(rerr,indBPM,errfun,W,-A,bpmerr.offsety);
            
        case 'x.y'
            disp('x,y misal, bpm ');
            
            % alignment
            errfun=@(r,po,er)setANYshift(r,po,1,er); % sets x errors
            rerr=seterrorwave(rerr,magindex,errfun,W,A,X0);
            
            errfun=@(r,po,er)setANYshift(r,po,3,er); % sets Y errors
            rerr=seterrorwave(rerr,magindex,errfun,W,A,Y0);
            
            % bpm
            errfun=@(r,po,er)setcellstruct(r,'Offset',po,er,1,1); % sets x bpm errors
            rerr=seterrorwave(rerr,indBPM,errfun,W,-A,bpmerr.offsetx);
            errfun=@(r,po,er)setcellstruct(r,'Offset',po,er,1,2); % sets Y bpm errors
            rerr=seterrorwave(rerr,indBPM,errfun,W,-A,bpmerr.offsety);
            
        case 'x.y.s.psi.theta.phi'
            disp('x,y,s misal, bpm and rotation errors (psi, phi and theta)');
            
            % rotation
            errfun=@(r,po,er)settilt_THERING_Dipole(r,po,er); % sets rotation about s
            rerr=seterrorwave(rerr,magindex,errfun,W,A,T0);
            
            errfun=@(r,po,er)setTiltAbout(r,po,'y',er); % sets rotation about Y
            rerr=seterrorwave(rerr,magindex,errfun,W,A,P0);
            
            errfun=@(r,po,er)setTiltAbout(r,po,'x',er); % sets rotation about X
            rerr=seterrorwave(rerr,magindex,errfun,W,A,R0);
            
            % alignment
            errfun=@(r,po,er)setANYshift(r,po,1,er); % sets x errors
            rerr=seterrorwave(rerr,magindex,errfun,W,A,X0);
            
            errfun=@(r,po,er)setANYshift(r,po,3,er); % sets Y errors
            rerr=seterrorwave(rerr,magindex,errfun,W,A,Y0);
            
            errfun=@(r,po,er)setDSerr(r,po,er); % sets S errors
            rerr=seterrorwave(rerr,magindex,errfun,W,A,S0);
            
            % bpm
            % errfun=@(r,po,er)setcellstruct(r,'Rotation',po,er,1); % sets psi bpm errors
            % rerr=seterrorwave(rerr,indBPM,errfun,W,-A);
            errfun=@(r,po,er)setcellstruct(r,'Offset',po,er,1,1); % sets x bpm errors
            rerr=seterrorwave(rerr,indBPM,errfun,W,-A,bpmerr.offsetx);
            errfun=@(r,po,er)setcellstruct(r,'Offset',po,er,1,2); % sets Y bpm errors
            rerr=seterrorwave(rerr,indBPM,errfun,W,-A,bpmerr.offsety);
            
        case 'gx'
            
            disp('girder X error');
            
            rerr=atsetwaveerrors(rerr,magindex,indBPM,W,A,'x');
            rerr=UniformGirderErrors(rerr);
            
        case 'gy'
            
            disp('girder Y error');
            
            rerr=atsetwaveerrors(rerr,magindex,indBPM,W,A,'y');
            rerr=UniformGirderErrors(rerr);
            
        case 'gpsi'
            
            disp('girder PSI error');
            
            rerr=atsetwaveerrors(rerr,magindex,indBPM,W,A,'psi');
            rerr=UniformGirderErrors(rerr);
            
        case 'gx.gy'
            
            disp('girder X Y error');
            
            rerr=atsetwaveerrors(rerr,magindex,indBPM,W,A,'x.y');
            rerr=UniformGirderErrors(rerr);
            
        case 'gx.gy.gpsi'
            disp('girder X Y PSI error');
            
            rerr=atsetwaveerrors(rerr,magindex,indBPM,W,A,'x.y.psi');
            rerr=UniformGirderErrors(rerr);
            
        case 'gx.gy.gpsi.x.y.psi'
            
            disp('girder wave and individual rand mag. X Y PSI error');
            
            rerr=atsetwaveerrors(rerr,magindex,indBPM,W,A,'x.y.psi');
            rerr=UniformGirderErrors(rerr);
            
            % get errors
            [X,Y,T]=GetMisalignments(rerr,magindex);
            
            seed=1;
            nsig=2;
            disp(['Girder wave amplitude: ' num2str(A*1e6) ' um'])
            betweengird=max(abs(diff(X)));
            disp(['Max distance between two adjacent girders: ' num2str(betweengird*1e6) ' um'])
            Ar=betweengird/sqrt(2)/2;% get distance between 2 girders
            disp(['Random errors max amplitude: ' num2str(Ar*1e6) ' um'])
            Ar=Ar/nsig; 
            % rotation
            errfun=@(r,po,er)settilt_THERING_Dipole(r,po,er); % sets rotation about s
            rerr=seterrorrand(rerr,magindex,errfun,seed,Ar,nsig,T);
            
            % alignment
            errfun=@(r,po,er)setANYshift(r,po,1,er); % sets x errors
            rerr=seterrorrand(rerr,magindex,errfun,seed*0,Ar,nsig,X);
            
            errfun=@(r,po,er)setANYshift(r,po,3,er); % sets Y errors
            rerr=seterrorrand(rerr,magindex,errfun,seed*0,Ar,nsig,Y);
            %uniform errors in sliced magnets
            rerr=UniformMagGroupsErrors(rerr);
            
        case 'gx.gy.x.y'
            
            disp('girder wave and individual rand mag. X Y error');
            
            rerr=atsetwaveerrors(rerr,magindex,indBPM,W,A,'x.y');
            rerr=UniformGirderErrors(rerr);
            
            % get errors
            [X,Y,T]=GetMisalignments(rerr,magindex);
            
            seed=1;
            nsig=2;
            disp(['Girder wave amplitude: ' num2str(A*1e6) ' um'])
            betweengird=max(abs(diff(X)));
            disp(['Max distance between two adjacent girders: ' num2str(betweengird*1e6) ' um'])
            Ar=betweengird/sqrt(2)/2;% get distance between 2 girders
            disp(['Random errors max amplitude: ' num2str(Ar*1e6) ' um'])
            Ar=Ar/nsig; 
           
            % alignment
            errfun=@(r,po,er)setANYshift(r,po,1,er); % sets x errors
            rerr=seterrorrand(rerr,magindex,errfun,seed,Ar,nsig,X);
            
            errfun=@(r,po,er)setANYshift(r,po,3,er); % sets Y errors
            rerr=seterrorrand(rerr,magindex,errfun,seed*0,Ar,nsig,Y);
            %uniform errors in sliced magnets
            rerr=UniformMagGroupsErrors(rerr);
              
        otherwise
            disp('type may be: x, y, s, psi, theta, phi, bpm');
            disp('             x.y, x.y.psi, x.y.s.psi, x.y.s.psi.theta.phi');
            disp('             gx, gy, gx.gy, gx.gy.gpsi,gx.gy.x.y, gx.gy.gpsi.x.y.psi');
    end
    
end

return