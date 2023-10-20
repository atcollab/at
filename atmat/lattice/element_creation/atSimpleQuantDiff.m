function elem = atSimpleQuantDiff(fname,varargin)
%SimpleQuantDiff creates a simple quantum difusion element
%
%ELEM=ATSIMPLEQUANTDIFF(FAMNAME,...)
%   FAMNAME:   family name
%
%ELEM=ATSIMPLEQUANTDIFF(FAMNAME,...,'Betax',BETAX,...)
%   BETAX:   Horizontal beta function. Default: 1.0
%
%ELEM=ATSIMPLEQUANTDIFF(FAMNAME,...,'Betay',BETAY,...)
%   BETAY:   Vertical beta function. Default: 1.0
%
%ELEM=ATSIMPLEQUANTDIFF(FAMNAME,...,'Emitx',EMITX,...)
%   EMITX:   Horizontal equilibrium emittance. Default: 0.0
%
%ELEM=ATSIMPLEQUANTDIFF(FAMNAME,...,'Emity',EMITY,...)
%   EMITY:   Vertical equilibrium emittance. Default: 0.0
%
%ELEM=ATSIMPLEQUANTDIFF(FAMNAME,...,'Espread',ESPREAD,...)
%   ESPREAD: Equilibrium momentum spread. Default: 0.0
%
%ELEM=ATSIMPLEQUANTDIFF(FAMNAME,...,'Taux',TAU_X,...)
%   TAU_X: Horizontal damping time. Default: 0.0
%
%ELEM=ATSIMPLEQUANTDIFF(FAMNAME,...,'Tauy',TAU_Y,...)
%   TAU_Y: Vertical damping time. Default: 0.0
%
%ELEM=ATSIMPLEQUANTDIFF(FAMNAME,...,'Tauz',TAU_Z,...)
%   TAU_Z: Longitudinal damping time. Default: 0.0
%
%ELEM=ATSIMPLEQUANTDIFF(FAMNAME,...,'U0',U0,...)
%   U0:     Energy loss [eV]. Default: 0.0

[rsrc,method]=decodeatargs({'SimpleQuantDiffPass'},varargin);
[method,rsrc]=getoption(rsrc,'PassMethod',method);
[cl,rsrc]=getoption(rsrc,'Class','SimpleQuantDiff');
[betax,rsrc]      = getoption(rsrc,'betax',1.0);
[betay,rsrc]      = getoption(rsrc,'betay',1.0);
[emitx,rsrc]      = getoption(rsrc,'emitx',0.0);
[emity,rsrc]      = getoption(rsrc,'emity',0.0);
[espread,rsrc]    = getoption(rsrc,'espread',0.0);
[taux,rsrc]       = getoption(rsrc,'taux',0.0);
[tauy,rsrc]       = getoption(rsrc,'tauy',0.0);
[tauz,rsrc]       = getoption(rsrc,'tauz',0.0);
[U0,rsrc]         = getoption(rsrc,'U0',0.0);
elem=atbaselem(fname,method,'Class',cl,'betax',betax,'betay',betay,...
    'emitx',emitx,'emity',emity,'espread',espread,...
    'taux',taux,'tauy',tauy,'tauz',tauz,...
    'U0',U0,rsrc{:});
end