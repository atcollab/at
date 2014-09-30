function elem=atdampMatElem(fname,ring,varargin)
%   atdampMatElem creates an element that applies the global damping matrix
%ELEM=ATDAMPMATELEM(FAMNAME,RING,CAVIPASS,BENDPASS,QUADPASS)
%
%FAMNAME:   family name
%RING:		initial AT structure, without radiation passmethods
%CAVIPASS:	pass method for cavities (default ThinCavityPass)
%BENDPASS:	pass method for bending magnets. Special values:
%           '' makes no change,
%           'auto' wille substitute 'Pass' with 'RadPass' in any method
%           (default: 'auto')
%QUADPASS:	pass method for quadrupoles
%           '' makes no change,
%           'auto' wille substitute 'Pass' with 'RadPass' in any method
%           (default: '')
%


ringrad=atradon(ring,varargin{:});

m66_norad=findm66(ring);
m66_rad=findm66(ringrad);

elem=atM66(fname,m66_norad\m66_rad);
end
