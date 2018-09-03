function cortab = atcorrectiontable(rerr,ring)
%ATCORRECTIONTABLE retrives a correction table from an AT lattice 
%
% input
% rcor: lattice with errors corrected
% rerr: lattice with same errors of rerr and without correction
%
%see also:

if length(ring)~=length(rerr)
error('lattices provided must be identical. first input with errors corrected, second input without correction');
end


% create correctors table
FamNames = atgetfieldvalues(rerr,1:length(ring),'FamName');
%DeviceNames = atgetfieldvalues(rerr,1:length(ring),'Device');
Lengths = atgetfieldvalues(rerr,1:length(ring),'Length');
KL0n = zeros(size(ring)); % hor steerer
KL0s = KL0n; % ver steerer
KL1n = KL0n; % normal quadruple
KL1s = KL0n; % skew quadrupole
KL2n = KL0n; % sextupole
KL3n = KL0n; % octupole

% values for lattice with errors no correction
ch0  = atgetfieldvalues(ring,b.indQCor,'PolynomB',{1,1});
cv0  = atgetfieldvalues(ring,b.indQCor,'PolynomB',{1,1});
cq0  = atgetfieldvalues(ring,b.indQCor,'PolynomB',{1,2});
csk0 = atgetfieldvalues(ring,b.indQCor,'PolynomA',{1,2});
cs0  = atgetfieldvalues(ring,b.indSext,'PolynomB',{1,3});
co0  = atgetfieldvalues(ring,b.indOctu,'PolynomB',{1,4});

% values for lattice with errors and correction
che  = atgetfieldvalues(rerr,b.indQCor,'PolynomB',{1,1});
cve  = atgetfieldvalues(rerr,b.indQCor,'PolynomB',{1,1});
cqe  = atgetfieldvalues(rerr,b.indQCor,'PolynomB',{1,2});
cske = atgetfieldvalues(rerr,b.indQCor,'PolynomA',{1,2});
cse  = atgetfieldvalues(rerr,b.indSext,'PolynomB',{1,3});
coe  = atgetfieldvalues(rerr,b.indOctu,'PolynomB',{1,4});

KL0n(b.indHCor) = (che-ch0).*Lengths(b.indHCor);
KL0s(b.indVCor) = (cve-cv0).*Lengths(b.indVCor);
KL1n(b.indQCor) = (cqe-cq0).*Lengths(b.indQCor);
KL1s(b.indSCor) = (cske-csk0).*Lengths(b.indSCor);
KL2n(b.indSext) = (cse-cs0).*Lengths(b.indSext);
KL3n(b.indOctu) = (coe-co0).*Lengths(b.indOctu);

% for ii = 1:length(DeviceNames)
%     if isempty(DeviceNames{ii})
%         DeviceNames{ii} = '-';
%     end
% end

FamNames = categorical(FamNames);
%DevNames = categorical(DeviceNames);
cortab = table(...%    DevNames,...
    FamNames,...
    Lengths,...
    KL0n,...
    KL0s,...
    KL1n,...
    KL1s,...
    KL2n,...
    KL3n);

end