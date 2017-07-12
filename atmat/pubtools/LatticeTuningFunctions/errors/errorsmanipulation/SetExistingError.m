function rerr=SetExistingError(rerr,magindex,indBPM,X0,Y0,S0,T0,R0,P0,bpm0)
% function SetExistingError(rerr,magindex,X0,Y0,S0,T0,R0,P0,bpm0)
% function SetExistingError(rerr,magindex,rerrset)
% rerrset is then a lattice with errors and the errors are taken from it
% sets a given set of errors.
%
% errors are overwriting the existing ones.
%
%see also: GetExistingErrors

nsig=2;

if iscell(X0)
    [X,Y,S,T,R,P,bpm]=GetExistingErrors(X0,magindex);
    rerr=SetExistingError(rerr,magindex,indBPM,X,Y,S,T,R,P,bpm);
else
    if ~isempty(magindex)
        % alignment
        if nargin>3
            errfun=@(r,po,er)setANYshift(r,po,1,er); % sets x errors
            rerr=seterrorrand(rerr,magindex,errfun,0,0,nsig,X0);
        end
        
        if nargin>4
            errfun=@(r,po,er)setANYshift(r,po,3,er); % sets Y errors
            rerr=seterrorrand(rerr,magindex,errfun,0,0,nsig,Y0);
        end
        
        if nargin>5
            errfun=@(r,po,er)atset_s_shift(r,po,er); % sets S errors
            rerr=seterrorrand(rerr,magindex,errfun,0,0,nsig,S0);
        end
        
        % rotation
        if nargin>6
            errfun=@(r,po,er)atsettiltdipole(r,po,er); % sets rotation about s
            rerr=seterrorrand(rerr,magindex,errfun,0,0,nsig,T0);
        end
        
        % represented by x and y displacements, already implemented
        if nargin>7
            %    errfun=@(r,po,er)setTiltAbout(r,po,'y',er); % sets rotation about Y
            %    rerr=seterrorrand(rerr,magindex,errfun,0,0,nsig,R0);
        end
        
        if nargin>8
            %    errfun=@(r,po,er)setTiltAbout(r,po,'x',er); % sets rotation about X
            %    rerr=seterrorrand(rerr,magindex,errfun,0,0,nsig,P0);
        end
    end
    
    % bpm
    if nargin>9
        errfun=@(r,po,er)setcellstruct(r,'Offset',po,er,1,1); % sets x bpm errors
        rerr=seterrorrand(rerr,indBPM,errfun,0,0,nsig,bpm0.offsetx);
        errfun=@(r,po,er)setcellstruct(r,'Offset',po,er,1,2); % sets Y bpm errors
        rerr=seterrorrand(rerr,indBPM,errfun,0,0,nsig,bpm0.offsety);
        errfun=@(r,po,er)setcellstruct(r,'Rotation',po,er,1,1); % sets rot bpm errors
        rerr=seterrorrand(rerr,indBPM,errfun,0,0,nsig,bpm0.rotation);
    end
end

return