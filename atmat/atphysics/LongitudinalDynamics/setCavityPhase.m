function newring=setCavityPhase(ring)
%SETCAVITYPHASE
%   Detailed explanation goes here

% Speed of light
CLIGHT=PhysConstant.speed_of_light_in_vacuum.value;

cavities=atgetcells(ring,'Frequency');
rfv=sum(atgetfieldvalues(ring(cavities),'Voltage'));
freq=atgetfieldvalues(ring(cavities),'Frequency');
U0=atgetU0(ring);
timelag=CLIGHT / (2*pi*freq) * asin(U0/rfv);
newring=atsetfieldvalues(ring,cavities,'TimeLag',timelag);
end
