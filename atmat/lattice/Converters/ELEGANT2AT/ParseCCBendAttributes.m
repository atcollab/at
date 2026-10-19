function R=ParseCCBendAttributes(R,name,key,value)
%PARSECCBENDATTRIBUTES Map supported elegant CCBEND fields without dropping physics.
% The legacy converter visits overlapping token pairs, including value tokens.
if ~contains(key,'='), return; end
key=upper(regexprep(key,'[=,\s]',''));
value=regexprep(value,'[=,\s]','');
if isempty(value), return; end
keys={'L','ANGLE','N_SLICES','INTEGRATION_ORDER','FRINGEMODEL','EDGE_ORDER', ...
    'EDGE_FLIP','YAW','OPTIMIZE_FSE','OPTIMIZE_DX','COMPENSATE_KN', ...
    'REFERENCE_CORRECTION','FSE','FSE_DIPOLE','FSE_QUADRUPOLE','FINT','FINT1','FINT2', ...
    'SYNCH_RAD','ISR','ISR1PART'};
fields={'Length','BendingAngle','NumIntSteps','IntegrationOrder','FringeModel','EdgeOrder', ...
    'EdgeFlip','Yaw','OptimizeFSE','OptimizeDX','CompensateKN', ...
    'ReferenceCorrection','FSE','FSEDipole','FSEQuadrupole','FringeInt1','FringeInt1','FringeInt2', ...
    'CCBSynchRad','CCBISR','QuantumSingleParticle'};
pos=find(strcmp(key,keys),1);
if ~isempty(pos)
    R=[R name '.' fields{pos} '=' value ';\n'];
    if strcmp(key,'FINT'), R=[R name '.FringeInt2=' value ';\n']; end
elseif ~isempty(regexp(key,'^K[1-8]$','once'))
    order=str2double(key(2));
    R=[R sprintf('%s.PolynomB(%d)=(%s)/%d;\n',name,order+1,value,factorial(order))];
    if order==1, R=[R name '.K=' value ';\n']; end
elseif strcmp(key,'HGAP')
    R=[R name '.FullGap=2*(' value ');\n'];
elseif ~isempty(regexp(key,'^FRINGE[12](K0|I0|K2|I1|K4|K5|K6|K7)$','once'))
    index=find(strcmp(key(8:end),{'K0','K2','K4','K5','K6','K7','I0','I1'}));
    % Initialize the complete array before assigning individual entries.
    R=[R sprintf('if ~isfield(%s,''fringeInt%s''), %s.fringeInt%s=zeros(1,8); end\n', ...
        name,key(7),name,key(7))];
    R=[R sprintf('%s.fringeInt%s(%d)=%s;\n',name,key(7),index,value)];
elseif strcmp(key,'TILT')
    R=[R name '=attiltelem(' name ',' value ');\n'];
elseif any(strcmp(key,{'GROUP','VERBOSE'}))
    return;
elseif ~isequal(str2double(value),0)
    error('AT:CCBendUnsupported','CCBEND %s: unsupported %s=%s',name,key,value);
end
end
