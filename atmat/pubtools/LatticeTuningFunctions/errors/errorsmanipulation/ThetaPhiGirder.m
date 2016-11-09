function rtp=ThetaPhiGirder(r,mag_gr)
%rtp=ThetaPhiGirder(r,mag_gr)
%
% sets misalignment to model theta, phi errors for magnets on girder
%
% mag_gr is the output of getMagGroupsFromGirderIndex(ring)
% 
%see also: GetExistingErrors setANYshift setTiltAbout seterrorrand

%    get girder centers
gm=cellfun(@(mg)gmisal(r,mg),mag_gr,'un',0);

magindex=[mag_gr{:}];
allmisg=[gm{:}];

% get existing alignment errors
%[X0,Y0]=GetExistingErrors(r,magindex);

rtp=r;
nsig=2;

errfun=@(r,po,er)setANYshift(r,po,1,er); % sets x errors
%rtp=seterrorrand(rtp,magindex,errfun,0,0,nsig,X0+allmisg(1,:));
rtp=seterrorrand(rtp,magindex,errfun,0,0,nsig,allmisg(1,:));

errfun=@(r,po,er)setANYshift(r,po,3,er); % sets Y errors
%rtp=seterrorrand(rtp,magindex,errfun,0,0,nsig,Y0+allmisg(2,:));
rtp=seterrorrand(rtp,magindex,errfun,0,0,nsig,allmisg(2,:));


return

function gm=gmisal(r,mg)

sg=findspos(r,mg);
gc=(max(sg)+min(sg))/2;

dg=sg-gc; % distance from girder center

tg=atgetfieldvalues(r,mg,'RotAboutX',{1,1}); % sets a vertical misalignment
pg=atgetfieldvalues(r,mg,'RotAboutY',{1,1}); % sets a horizontal misalignment

gm=zeros(length(mg),2);

if ~isempty(tg)
    gm(:,2)=dg'.*sin(tg);
end
if ~isempty(pg)
    gm(:,1)=dg'.*sin(pg);
end

gm=gm';
return


