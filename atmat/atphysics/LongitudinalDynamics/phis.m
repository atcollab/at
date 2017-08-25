function phase = phis(U0MeV,VrfMV)
% phase = phis(U0MeV,VrfMV)
%
% this function returns the synchronous phase in radians
% input:
% U0MeV is energy loss per turn in MeV
% VrfMV is the RF voltage in MV

phase=pi - asin(U0MeV./VrfMV);

end