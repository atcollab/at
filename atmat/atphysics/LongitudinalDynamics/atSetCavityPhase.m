function newring=atSetCavityPhase(ring,varargin)
%SETCAVITYPHASE     Set the TimeLag attribute of RF cavities
%
%NEWRING=SETCAVITYPHASE(RING)
%   Adjust the TimeLag attribute of RF cavities based on frequency,
%   voltage and energy loss per turn, so that the synchronous phase is zero.
%   An error occurs if all cavities do not have the same frequency.
%
%NEWRING=SETCAVITYPHASE(...,'refpts',CAVPTS)
%   CAVPTS is the location of RF cavities. This allows to ignore harmonic
%   cavities.
%
%NEWRING=SETCAVITYPHASE(...,'method',METHOD)
%   Choose the method for computing the energy loss per turn
%
% METHOD:   'integral': (default) The losses are obtained from
%                       Losses = Cgamma / 2pi * EGeV^4 * I2
%                       Takes into account bending magnets and wigglers.
%           'tracking': The losses are obtained by tracking without cavities.
%                       Needs radiation ON, takes into account all radiating elements.

clight=PhysConstant.speed_of_light_in_vacuum.value;

[cavities,varargs]=getoption(varargin,'refpts',atgetcells(ring,'Frequency'));
rfv=sum(atgetfieldvalues(ring(cavities),'Voltage'));
freq=unique(atgetfieldvalues(ring(cavities),'Frequency'));
if length(freq) > 1
    error('AT:NoFrequency','RF frequency not equal for all cavities');
end
warning('AT:CavityTimeLag',...
    ['\nThis function modifies the time reference\n',...
       'This should be avoided, you have been warned\n']);
U0=atgetU0(ring,'periods',1,varargs);
timelag=clight / (2*pi*freq) * asin(U0/rfv);
newring=atsetfieldvalues(ring,cavities,'TimeLag',timelag);
end
