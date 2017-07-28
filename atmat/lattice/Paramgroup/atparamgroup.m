% PARAMETER GROUP in AT is a general way
% to define a functional relationship between a single
% scalar (control) paramer and mulltiple
% (dependent) paramters in the accelerator lattice.
%
% Examples:
% 1. Simultaneously varying the strength of all 
%    quadrupole magnets in the same family
% 2. Moving a magnet longitudinally. This requires 
%    changing the lengh of the drift spaces before and 
%    after the element while constrining the total length.
% 3. Setting the R1, R2, T1 and  T2 field for misaligned magnets.
%    All elements of these matrixes are parametrized with  
%    three angles and two displacemets.
%
% To use Parameter Group - first define a 
% MATLAB structure array that describes it.
% Each element in this structure arraay for each of the
% dependent paramanters
% 
% This structure must have the following fields:
%   ElemIndex  - Index of an element in the lattice
%   FieldName  - Name of the field to be modified
%   Function   - MATLAB inline objct the defines the
%                functional relationship with
%                the control parameter
%   FieldIndex - If the data in the field given by FieldName 
%                is a matrix or a vector, FieldIndex specifies
%                inexes of elements to be varied
%   SavedValue - saved values that can be restored with RESTOREPARAMGROUP
%   Args       - Cell array of aditional arguments to be passed 
%                to the function
%
% After the parameter group structure is defined
% it is used as an argument to AT function SETPARAMGROUP:
% NEWRING = setparamgroup(RING,PARAMGROUPSTRUCT,PARAMVALUE)

% Example 1
clear all
spear2;
QFI = findcells(THERING,'FamName','QF');
InitialValue = THERING{QFI(1)}.K;

% Create parameter group structure P
P = struct('ElemIndex',num2cell(QFI),'FieldName','K',...
    'Function',inline('x'));
% In this case the field 'K' in quadrupoles is a scalar
[P.FieldIndex]=deal({1,1});
[P.Args]=deal({});


% Inspect P:
disp(P(1))
% Change the K field by 0.1%
NewValue = InitialValue*1.001;
THERING = setparamgroup(THERING,P,NewValue);


% Example 2
clear all
spear2;
QFI = findcells(THERING,'FamName','QF');
FirstQF = QFI(1);

L1 = THERING{FirstQF-1}.Length;
L2 = THERING{FirstQF+1}.Length;

P(1).ElemIndex=FirstQF-1;
P(2).ElemIndex=FirstQF+1;

P(1).Function=inline('P1+x',1);
P(1).Args={L1};

P(2).Function=inline('P1-x',1);
P(2).Args={L2};

[P.FieldName]=deal('Length');
[P.FieldIndex]=deal({1,1});

THERING = setparamgroup(THERING,P,0.1);