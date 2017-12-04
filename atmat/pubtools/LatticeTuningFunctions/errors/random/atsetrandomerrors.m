function rerr=atsetrandomerrors(...
    rerr,...
    magindex,...
    indBPM,...
    seed,...
    sigma,...
    nsig,...
    type)
%function rerr=atsetrandomerrors(...
% rerr, %lattice
% magindex, %indexes of magnts to set errors
% indBPM,  %index of bpm to mode reference using offsets
%  seed    % error seed
%  sigma   % number of sigma for truncation
% type)    % type of error (see list)
%
% Set random errors for a given seed.
% Considers also magnets grouped by MagNum field and
% girders specified by GS and GE markers.
% BPMs are moved with girders.
%
% If magindx is a cellarray of index vectors, 
% then also sigma and type must be and they will be applyed summing. 
% 
% type may be: x, y, s, psi (roll), theta (yaw), phi (pitch), (individual magnets)
%               bpm (bpm.offset)    
%               bpm.scale, bpm.read   
%              x.y, x.y.psi, x.y.s.psi,
%              x.y.s.psi.theta.phi
%              gx, gy, gpsi, gtheta, gphi,  gx.gy,             (girder)
%              gx.gy.gpsi, gx.gy.gpsi.x.y.psi
%              dpb1, dpb2, dpb3, dpb4           (PolynomB random error)
%
%
%see also: atsetshift atset_s_shift setTiltAbout seterrorrand
%UniformMagGroupsErrors UniformGirderErrors ApplyErrorWave

r0=rerr; % store nominal lattice at this stage

if iscell(magindex)% if list of errors is given
    if length(magindex)~=length(sigma) || length(magindex)~=length(type)
        error('magindex, sigma, type, must be cell array of the same length, or numeric.')
    end
    
    % loop magindex to apply errors
    for im=1:length(magindex)
        
        if im>0, seed=0; end; % do not change or reset seed anymore (seterrorrand)
        
        disp([...
            'magindex ' num2str(im)...
            ', type:' type{im}...
            ', sig:' num2str(sigma{im})...
            ', seed: ' num2str(seed)]...
            );
        
        rerr=atsetrandomerrors(...
            rerr,...
            magindex{im},...
            indBPM,...
            seed,...
            sigma{im},...
            nsig,...
            type{im});
        
    end
    
else
    
    % set error wave to zero if no error existing
    [X0,Y0,S0,T0,R0,P0,bpm0]=GetExistingErrors(rerr,magindex);
    
    if isnan(Y0)
    disp('NAN Y0! There are NaN in the previous set of errors!')
    end
    
    % initialize errors if not inizialized
    %if (std(X0)==0 && std(Y0)==0) 
    if any(X0==0 & Y0==0) 
        % sets all field of T1 to 0
        errfun=@(r,po,er)setANYshift(r,po,1:6,er);
        rerr=seterrorrand(rerr,magindex(X0==0 & Y0==0),errfun,0,0,nsig);
    end
    
    % set error to required amplitude and wavelength
    switch type
        case 'zero'
            disp('setting all errors to 0');
            
            % sets all field of T1 to 0
            errfun=@(r,po,er)setANYshift(r,po,1:6,er);
            rerr=seterrorrand(rerr,magindex,errfun,1,0,nsig);
            errfun=@(r,po,er)atset_s_shift(r,po,er); % sets DS errors to zero
            rerr=seterrorrand(rerr,magindex,errfun,0,0,nsig);
            
            % rotation
            errfun=@(r,po,er)atsettiltdipole(r,po,er); % sets rotation about s
            rerr=seterrorrand(rerr,magindex,errfun,0,0,nsig);
            errfun=@(r,po,er)setTiltAbout(r,po,'y',er); % sets rotation about Y
            rerr=seterrorrand(rerr,magindex,errfun,0,0,nsig);
            errfun=@(r,po,er)setTiltAbout(r,po,'x',er); % sets rotation about X
            rerr=seterrorrand(rerr,magindex,errfun,0,0,nsig);
            
            % bpm
            errfun=@(r,po,er)setcellstruct(r,'Offset',po,er,1,1); % sets x bpm Offset errors
            rerr=seterrorrand(rerr,indBPM,errfun,0,0,nsig);
            errfun=@(r,po,er)setcellstruct(r,'Offset',po,er,1,2); % sets Y bpm errors
            rerr=seterrorrand(rerr,indBPM,errfun,0,0,nsig);
            errfun=@(r,po,er)setcellstruct(r,'Scale',po,er,1,1); % sets x bpm Scale errors
            rerr=seterrorrand(rerr,indBPM,errfun,0,0,nsig);
            errfun=@(r,po,er)setcellstruct(r,'Scale',po,er,1,2); % sets Y bpm errors
            rerr=seterrorrand(rerr,indBPM,errfun,0,0,nsig);
            errfun=@(r,po,er)setcellstruct(r,'Reading',po,er,1,1); % sets x bpm Reading errors
            rerr=seterrorrand(rerr,indBPM,errfun,0,0,nsig);
            errfun=@(r,po,er)setcellstruct(r,'Reading',po,er,1,2); % sets Y bpm errors
            rerr=seterrorrand(rerr,indBPM,errfun,0,0,nsig);
            errfun=@(r,po,er)setcellstruct(r,'Rotation',po,er,1,1); % sets x bpm Rotation errors
            rerr=seterrorrand(rerr,indBPM,errfun,0,0,nsig);
            
        case 'dpb1'
            
            disp('PolynomB(1) error');
            rerr=setFieldIntegralError(r0,rerr,magindex,1,nsig,sigma);
            
        case 'dpb2'
            disp('PolynomB(2) error');
            
            rerr=setFieldIntegralError(r0,rerr,magindex,2,nsig,sigma);
            
        case 'dpb3'
            disp('PolynomB(3) error');
            
            rerr=setFieldIntegralError(r0,rerr,magindex,3,nsig,sigma);
            
        case 'dpb4'
            disp('PolynomB(4) error');
            
            rerr=setFieldIntegralError(r0,rerr,magindex,4,nsig,sigma);
            
        case 'x'
            
            disp('X error');
           
            errfun=@(r,po,er)setANYshift(r,po,1,er); % sets X errors
            rerr=seterrorrand(rerr,magindex,errfun,seed,sigma,nsig,X0);
            
            %uniform errors in sliced magnets
            rerr=UniformMagGroupsErrors(rerr);
            
        case 'y'
            
            disp('Y error');
            
            errfun=@(r,po,er)setANYshift(r,po,3,er); % sets Y errors
            rerr=seterrorrand(rerr,magindex,errfun,seed,sigma,nsig,Y0);
            %uniform errors in sliced magnets
            rerr=UniformMagGroupsErrors(rerr);
            
        case 's' % longitudinal displacement
            
            disp('S error');
            
            errfun=@(r,po,er)atset_s_shift(r,po,er); % sets Y errors
            rerr=seterrorrand(rerr,magindex,errfun,seed,sigma,nsig,S0);
            %uniform errors in sliced magnets
            rerr=UniformMagGroupsErrors(rerr);
            
        case 'phi' % about x
            
            disp('phi error');
            
            errfun=@(r,po,er)setTiltAbout(r,po,'x',er); % sets Y errors
            rerr=seterrorrand(rerr,magindex,errfun,seed,sigma,nsig,R0);
            %uniform errors in sliced magnets
            rerr=UniformMagGroupsErrors(rerr);
            
        case 'theta' % about y
            
            disp('theta error');
            
            errfun=@(r,po,er)setTiltAbout(r,po,'y',er); % sets Y errors
            rerr=seterrorrand(rerr,magindex,errfun,seed,sigma,nsig,P0);
            %uniform errors in sliced magnets
            rerr=UniformMagGroupsErrors(rerr);
            
        case 'psi' % about s
            
            disp('s-asix rotation error');
            
            errfun=@(r,po,er)setTiltAbout(r,po,'s',er); % sets Y errors
            rerr=seterrorrand(rerr,magindex,errfun,seed,sigma,nsig,T0);
            %uniform errors in sliced magnets
            rerr=UniformMagGroupsErrors(rerr);
            
        case 'bpm'
            
            disp('bpm offset error');
            % bpm
            % errfun=@(r,po,er)setcellstruct(r,'Rotation',po,er,1); % sets psi bpm errors
            % rerr=seterrorrand(rerr,indBPM,errfun,W,-A);
          
            errfun=@(r,po,er)setcellstruct(r,'Offset',po,er,1,1); % sets x bpm errors
            rerr=seterrorrand(rerr,indBPM,errfun,seed,-sigma,nsig,bpm0.offsetx);
            errfun=@(r,po,er)setcellstruct(r,'Offset',po,er,1,2); % sets Y bpm errors
            rerr=seterrorrand(rerr,indBPM,errfun,seed*0,-sigma,nsig,bpm0.offsety);
            
         case 'bpm.offset'
            
            disp('bpm offset error');
            % bpm
            % errfun=@(r,po,er)setcellstruct(r,'Rotation',po,er,1); % sets psi bpm errors
            % rerr=seterrorrand(rerr,indBPM,errfun,W,-A);
          
            errfun=@(r,po,er)setcellstruct(r,'Offset',po,er,1,1); % sets x bpm errors
            rerr=seterrorrand(rerr,indBPM,errfun,seed,-sigma,nsig,bpm0.offsetx);
            errfun=@(r,po,er)setcellstruct(r,'Offset',po,er,1,2); % sets Y bpm errors
            rerr=seterrorrand(rerr,indBPM,errfun,seed*0,-sigma,nsig,bpm0.offsety);
            
        case 'bpm.scale'
            
            disp('bpm scale error (1+...)');
            % bpm
            % errfun=@(r,po,er)setcellstruct(r,'Rotation',po,er,1); % sets psi bpm errors
            % rerr=seterrorrand(rerr,indBPM,errfun,W,-A);
          
            errfun=@(r,po,er)setcellstruct(r,'Scale',po,1+er,1,1); % sets x bpm errors
            rerr=seterrorrand(rerr,indBPM,errfun,seed,-sigma,nsig,zeros(size(bpm0.offsetx)));
            errfun=@(r,po,er)setcellstruct(r,'Scale',po,1+er,1,2); % sets Y bpm errors
            rerr=seterrorrand(rerr,indBPM,errfun,seed*0,-sigma,nsig,zeros(size(bpm0.offsety)));
            
        case 'bpm.read'
            
            disp('bpm reading error');
            % bpm
            % errfun=@(r,po,er)setcellstruct(r,'Rotation',po,er,1); % sets psi bpm errors
            % rerr=seterrorrand(rerr,indBPM,errfun,W,-A);
          
            errfun=@(r,po,er)setcellstruct(r,'Reading',po,er,1,1); % sets x bpm errors
            rerr=seterrorrand(rerr,indBPM,errfun,seed,-sigma,nsig,zeros(size(bpm0.offsetx)));
             
        case 'x.y.s.psi'
            disp('x,y,s misal, bpm and s-asix rotation error (psi only, no phi and theta)');
            
            % rotation
            errfun=@(r,po,er)atsettiltdipole(r,po,er); % sets rotation about s
            rerr=seterrorrand(rerr,magindex,errfun,seed,sigma,nsig,T0);
            
            % alignment
            errfun=@(r,po,er)setANYshift(r,po,1,er); % sets x errors
            rerr=seterrorrand(rerr,magindex,errfun,seed*0,sigma,nsig,X0);
            
            errfun=@(r,po,er)setANYshift(r,po,3,er); % sets Y errors
            rerr=seterrorrand(rerr,magindex,errfun,seed*0,sigma,nsig,Y0);
            
            %uniform errors in sliced magnets
            rerr=UniformMagGroupsErrors(rerr);
            
            errfun=@(r,po,er)atset_s_shift(r,po,er); % sets S errors
            rerr=seterrorrand(rerr,magindex,errfun,seed*0,sigma,nsig,S0);
            
            
        case 'x.y.psi'
            disp('x,y misal, bpm and s-asix rotation error (psi only, no phi and theta)');
            
            % rotation
            errfun=@(r,po,er)atsettiltdipole(r,po,er); % sets rotation about s
            rerr=seterrorrand(rerr,magindex,errfun,seed,sigma,nsig,T0);
            
            % alignment
            errfun=@(r,po,er)setANYshift(r,po,1,er); % sets x errors
            rerr=seterrorrand(rerr,magindex,errfun,seed*0,sigma,nsig,X0);
            
            errfun=@(r,po,er)setANYshift(r,po,3,er); % sets Y errors
            rerr=seterrorrand(rerr,magindex,errfun,seed*0,sigma,nsig,Y0);
            %uniform errors in sliced magnets
            rerr=UniformMagGroupsErrors(rerr);
            
            
        case 'x.y'
            disp('x,y misal');
            
            % alignment
            errfun=@(r,po,er)setANYshift(r,po,1,er); % sets x errors
            rerr=seterrorrand(rerr,magindex,errfun,seed,sigma,nsig,X0);
            
            errfun=@(r,po,er)setANYshift(r,po,3,er); % sets Y errors
            rerr=seterrorrand(rerr,magindex,errfun,0,sigma,nsig,Y0);
            %uniform errors in sliced magnets
            rerr=UniformMagGroupsErrors(rerr);
            
            
        case 'x.y.s.psi.theta.phi'
            disp('x,y,s misal, bpm and rotation errors (psi, phi and theta)');
            
            % rotation
            errfun=@(r,po,er)atsettiltdipole(r,po,er); % sets rotation about s
            rerr=seterrorrand(rerr,magindex,errfun,seed,sigma,nsig,T0);
            
            errfun=@(r,po,er)setTiltAbout(r,po,'y',er); % sets rotation about Y
            rerr=seterrorrand(rerr,magindex,errfun,seed*0,sigma,nsig,R0);
            
            errfun=@(r,po,er)setTiltAbout(r,po,'x',er); % sets rotation about X
            rerr=seterrorrand(rerr,magindex,errfun,seed*0,sigma,nsig,P0);
            
            % alignment
            errfun=@(r,po,er)setANYshift(r,po,1,er); % sets x errors
            rerr=seterrorrand(rerr,magindex,errfun,seed*0,sigma,nsig,X0);
            
            errfun=@(r,po,er)setANYshift(r,po,3,er); % sets Y errors
            rerr=seterrorrand(rerr,magindex,errfun,seed*0,sigma,nsig,Y0);
            
            %uniform errors in sliced magnets before placing s errors
            rerr=UniformMagGroupsErrors(rerr);
            
            errfun=@(r,po,er)atset_s_shift(r,po,er); % sets S errors
            rerr=seterrorrand(rerr,magindex,errfun,seed*0,sigma,nsig,S0);
            
        case 'gx'
            
            disp('girder X error');
            
            % set all errors to 0
            rZ=atsetrandomerrors(r0,magindex,indBPM,0,0,nsig,'zero'); 

            % assign girder errors on zero err lattice
            rgerr=atsetrandomerrors(rZ,magindex,indBPM,seed,sigma,nsig,'x');
            [rgerr,mag_gr]=UniformGirderErrors(rgerr); % unifroms all errors!
            
            % sum errors 
            magindex=[mag_gr{:}];
            rerr=SumErrors(rerr,rgerr,magindex,indBPM);
            
        case 'gy'
           
            disp('girder Y error');
            
            % set all errors to 0
            rZ=atsetrandomerrors(r0,magindex,indBPM,0,0,nsig,'zero'); 

            % assign girder errors on zero err lattice
            rgerr=atsetrandomerrors(rZ,magindex,indBPM,seed,sigma,nsig,'y');
            [rgerr,mag_gr]=UniformGirderErrors(rgerr); % unifroms all errors!
            
            % sum errors 
            magindex=[mag_gr{:}];
            rerr=SumErrors(rerr,rgerr,magindex,indBPM);
            
        case 'gpsi'
            
            disp('girder PSI error');
            
            % set all errors to 0
            rZ=atsetrandomerrors(r0,magindex,indBPM,0,0,nsig,'zero'); 

            % assign girder errors on zero err lattice
            rgerr=atsetrandomerrors(rZ,magindex,indBPM,seed,sigma,nsig,'psi');
            [rgerr,mag_gr]=UniformGirderErrors(rgerr); % unifroms all errors!
            
            % get errors on all gider magnets
            magindex=[mag_gr{:}];
            rerr=SumErrors(rerr,rgerr,magindex,indBPM);
        
        case 'gtheta'
            
            disp('girder THETA error');
            
            % set all errors to 0
            rZ=atsetrandomerrors(r0,magindex,indBPM,0,0,nsig,'zero'); 

            % assign girder errors on zero err lattice
            rgerr=atsetrandomerrors(rZ,magindex,indBPM,seed,sigma,nsig,'theta');
            [rgerr,mag_gr]=UniformGirderErrors(rgerr); % unifroms all errors!
            
            rgerr=ThetaPhiGirder(rgerr,mag_gr);
            
            % get errors on all gider magnets
            magindex=[mag_gr{:}];
            rerr=SumErrors(rerr,rgerr,magindex,indBPM);
              
        case 'gphi'
            disp('girder phi error');
            
            % set all errors to 0
            rZ=atsetrandomerrors(r0,magindex,indBPM,0,0,nsig,'zero'); 

            % assign girder errors on zero err lattice
            rgerr=atsetrandomerrors(rZ,magindex,indBPM,seed,sigma,nsig,'phi');
            [rgerr,mag_gr]=UniformGirderErrors(rgerr); % unifroms all errors!
            
            rgerr=ThetaPhiGirder(rgerr,mag_gr);
            
            % sum errors
            magindex=[mag_gr{:}];
            rerr=SumErrors(rerr,rgerr,magindex,indBPM);
              
            
        case 'gx.gy'
            
            disp('girder X Y error');
            
            rerr=atsetrandomerrors(rerr,magindex,indBPM,seed,sigma,nsig,'gy');
            rerr=atsetrandomerrors(rerr,magindex,indBPM,seed+1000,sigma,nsig,'gx');

        case 'gx.gy.gpsi'
            disp('girder X Y PSI error');
            rerr=atsetrandomerrors(rerr,magindex,indBPM,seed,sigma,nsig,'gx.gy');
            rerr=atsetrandomerrors(rerr,magindex,indBPM,seed+2000,sigma,nsig,'gpsi');

        otherwise
            disp('type may be: x, y, s, psi, theta, phi, bpm');
            disp('             x.y, x.y.psi, x.y.s.psi, x.y.s.psi.theta.phi');
            disp('             gx, gy, gpsi, gtheta,gphi, gx.gy, gx.gy.gpsi,');
            disp('             dpb1, dpb2, dpb3, dpb4');
    end
    
%     [Xe,Ye,Se,Te,Re,Pe]=GetExistingErrors(rerr,magindex);
%     disp(['X: ' num2str(std(X0),'%2.2e') ' -> ' num2str(std(Xe),'%2.2e')]);
%     disp(['Y: ' num2str(std(Y0),'%2.2e') ' -> ' num2str(std(Ye),'%2.2e')]);
%     disp(['S: ' num2str(std(S0),'%2.2e') ' -> ' num2str(std(Se),'%2.2e')]);
%     disp(['T: ' num2str(std(T0),'%2.2e') ' -> ' num2str(std(Te),'%2.2e')]);
%     disp(['R: ' num2str(std(R0),'%2.2e') ' -> ' num2str(std(Re),'%2.2e')]);
%     disp(['P: ' num2str(std(P0),'%2.2e') ' -> ' num2str(std(Pe),'%2.2e')]);

end% if cell

%rerr=setBpmOffsetOnDipoleRef(rerr);

return


