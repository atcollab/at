function P = mkparamgroup(LATTICE,varargin)
%MKPARAMGROUP simplifies creation of AT parameter groups
% It group one or more elements in the
% same family and simultaneously vary
% 
% MKPARAMGROUP(LATTICE,ELEMINDEX,PARAMSTR)
% MKPARAMGROUP(LATTICE,FAMNAMESTR,PARAMSTR)
% MKPARAMGROUP(LATTICE,FAMNAMESTR,KIDNUM,PARAMSTR)
% 
% LATTICE 
% FAMNAMESTR
% 
%
% PARAMSTR: 'TILT','K1','K2','K3'
% wjc 2-09-04 changed index 'i' to 'k'

if isnumeric(varargin{1})
    if ~((nargin==3)& ischar(varargin{2}))
        error('The third argument must be a string')
    else
        INDEX = varargin{1};
        KIDNUM = 1:length(INDEX);
        PARAMSTR = varargin{2};
    end
else
    FAMNAMESTR = varargin{1};
    INDEX = findcells(LATTICE,'FamName',FAMNAMESTR);
    if(isempty(INDEX))
        error(['No elements that belong to the family ''',FAMNAMESTR,...
                ''' found in the lattice ',inputname(1)]);  
    end
    if isnumeric(varargin{2})
        KIDNUM = varargin{2};
        PARAMSTR = varargin{3};
    else
        KIDNUM = 1:length(INDEX);
        PARAMSTR = varargin{2};
    end 
end

switch lower(PARAMSTR)
case {'k1','k'}
    if ~isfield(LATTICE{INDEX(KIDNUM(1))},'K')
        error('Element ',int2str(KIDNUM(i)),' does not have field ''K''');
    end
    P1 = struct('ElemIndex',num2cell(INDEX(KIDNUM)),'FieldName','K','Function',inline('x'));
    [P1.FieldIndex]=deal({1,1});
    [P1.Args]=deal({});
    
    if ~isfield(LATTICE{INDEX(KIDNUM(1))},'PolynomB')
        error('Element ',int2str(KIDNUM(1)),' does not have field ''PolynomB''');
    end
    P2 = struct('ElemIndex',num2cell(INDEX(KIDNUM)),'FieldName','PolynomB','Function',inline('x'));
    [P2.FieldIndex]=deal({1,2});
    [P2.Args]=deal({});
    
    for k = 1:length(KIDNUM)
        P1(k).SavedValue = LATTICE{INDEX(KIDNUM(k))}.K;
        P2(k).SavedValue = LATTICE{INDEX(KIDNUM(k))}.PolynomB(2);
    end
    P = reshape([P1;P2],1,2*length(P1));
case 'k2'
    
    if ~isfield(LATTICE{INDEX(KIDNUM(1))},'PolynomB')
        error('Element ',int2str(KIDNUM(1)),' does not have field ''PolynomB''');
    end
    P = struct('ElemIndex',num2cell(INDEX(KIDNUM)),'FieldName','PolynomB','Function',inline('x'));
    [P.FieldIndex]=deal({1,3});
    [P.Args]=deal({}); 
    for k = 1:length(KIDNUM)
        P(k).SavedValue = LATTICE{INDEX(KIDNUM(k))}.PolynomB(3);
    end    
case 'k3'
    if ~isfield(LATTICE{INDEX(KIDNUM(1))},'PolynomB')
        error('Element ',int2str(KIDNUM(1)),' does not have field ''PolynomB''');
    end
    P = struct('ElemIndex',num2cell(INDEX(KIDNUM)),'FieldName','PolynomB','Function',inline('x'));
    [P.FieldIndex]=deal({1,4});
    [P.Args]=deal({}); 
    for k = 1:length(KIDNUM)
        P(k).SavedValue = LATTICE{INDEX(KIDNUM(k))}.PolynomB(4);
    end  

case 'tilt'
    if ~isfield(LATTICE{INDEX(KIDNUM(1))},'R1')
        error('Element ',int2str(KIDNUM(1)),' does not have field ''R1''');
    end
    if ~isfield(LATTICE{INDEX(KIDNUM(1))},'R2')
        error('Element ',int2str(KIDNUM(1)),' does not have field ''R2''');
    end
    
    P1 = struct('ElemIndex',num2cell(INDEX(KIDNUM)),'FieldName','R1','Function',inline('mksrollmat(x)'));
    [P1.FieldIndex]=deal({1:6,1:6});
    [P1.Args]=deal({});
    
    P2 = struct('ElemIndex',num2cell(INDEX(KIDNUM)),'FieldName','R2','Function',inline('mksrollmat(-x)'));
    [P2.FieldIndex]=deal({1:6,1:6});
    [P2.Args]=deal({});
    
    for k = 1:length(KIDNUM)
        P1(k).SavedValue = LATTICE{INDEX(KIDNUM(k))}.R1;
        P2(k).SavedValue = LATTICE{INDEX(KIDNUM(k))}.R2;
    end
    P = reshape([P1;P2],1,2*length(P1)); 
    
    case {'k1','k'}
    if ~isfield(LATTICE{INDEX(KIDNUM(1))},'K')
        error('Element ',int2str(KIDNUM(k)),' does not have field ''K''');
    end
    
    P1 = struct('ElemIndex',num2cell(INDEX(KIDNUM)),'FieldName','K','Function',inline('x'));
    [P1.FieldIndex]=deal({1,1});
    [P1.Args]=deal({});
    
    if ~isfield(LATTICE{INDEX(KIDNUM(1))},'PolynomB')
        error('Element ',int2str(KIDNUM(1)),' does not have field ''PolynomB''');
    end
    P2 = struct('ElemIndex',num2cell(INDEX(KIDNUM)),'FieldName','PolynomB','Function',inline('x'));
    [P2.FieldIndex]=deal({1,2});
    [P2.Args]=deal({});
    
    for k = 1:length(KIDNUM)
        P1(k).SavedValue = LATTICE{INDEX(KIDNUM(k))}.K;
        P2(k).SavedValue = LATTICE{INDEX(KIDNUM(k))}.PolynomB(2);
    end
    P = reshape([P1;P2],1,2*length(P1));
    
    case {'s','s1'}
    if ~isfield(LATTICE{INDEX(KIDNUM(1))},'PolynomA')
        error('Element ',int2str(KIDNUM(1)),' does not have field ''PolynomA''');
    end
    P = struct('ElemIndex',num2cell(INDEX(KIDNUM)),'FieldName','PolynomA','Function',inline('x'));
    [P.FieldIndex]=deal({1,2});
    [P.Args]=deal({}); 
    for k = 1:length(KIDNUM)
        P(k).SavedValue = LATTICE{INDEX(KIDNUM(k))}.PolynomA(2);
    end  

end

