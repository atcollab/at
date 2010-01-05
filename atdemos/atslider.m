function fig = atslider(varargin)
%ATSLIDER is an example of a GUI control of multiple parameters in THERING
% by mapping from one-dimensional control parameter space.
% The value of the control parameter is set
% with a slider or typed in a text window. When it is modified
% each of the controlled variables is modified accordingly
%
% HANDLE = ATSLIDER(KnobData,KnobName)
% 	creates a new knob figure identified with KnobName
%	with  a slider and editable text. The initial value is set to 0
%  
%  KnobData is MATLAB structure arrray where each element corresponds to
%       one controlled paramater in THERING and controls what gets modified
%       in the ring when the slider value changes
%	     
%       Position 			Index of an element to be modified in THERING
%	     FieldName       Name of the field in the element data structure 
%                       whos value will bi modified
%	     M,N             When a field is an array M,N index an element (PolynomA)
%                       When a field is a scalar set M=N=1 
%	     Weight          Mapping coefficient between slider position and controlled value
%
%	HANDLE = ATSLIDER(KnobData,KnobName,COMMAND)
%			evaluates COMMAND in MATLAB base workspace
%        for example try COMMAND = 'plotbeta'  for dynamically updating betafunctions
%
% ATSLIDER(action) with proper 'action' stringis is recursively called 
% from inside the ATSLIDER function to evaluate callbacks
%		'set'
%		'reset' 

global THERING

if nargin >1 % initial call
	
 

KnobData = varargin{1};
KnobName = varargin{2};

% make a copy of original values of all controlled parametes
NumVar = length(KnobData);
I = [KnobData.Position];
F = [KnobData.FieldName];
M = [KnobData.M];
N = [KnobData.N];
W = [KnobData.Weight];
for i = 1:NumVar
   OriginalValues(i) = getfield(THERING{I(i)},F(i),{M(i),N(i)});
end


h0 = figure('Color', [0.8 0.8 0.8], ...
	'FileName','D:\MATLABR11\work\atslider.m', ...
   'HandleVisibility', 'Callback' , ...
	'PaperPosition',[18 180 576 432], ...
	'PaperUnits','points', ...
	'Position',[520 400 300 140], ...
	'Tag','Fig1', ...
	'ToolBar','none');
if nargout > 0 
   fig = h0; 
end

s1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.752941176470588 0.752941176470588 0.752941176470588], ...
	'ListboxTop',0, ...
	'Position',[20 20 100 15], ...
   'Style','slider', ...
   'Max',1, 'Min', -1, ...
   'Callback','atslider set', ...
   'Tag','Slider1');

e1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'ListboxTop',0, ...
	'Position',[140 20 65 15], ...
   'Style','edit', ...
   'Callback','atslider set', ...
   'Tag','EditText1');
 

if nargin == 2 
	set(h0,'UserData',{e1 s1 KnobData OriginalValues});
elseif nargin == 3
	set(h0,'UserData',{e1 s1 KnobData OriginalValues varargin{3}});
end


t1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.8 0.8 0.8], ...
	'FontSize',12, ...
	'ListboxTop',0, ...
	'Position',[20 60 165 20], ...
	'String',KnobName, ...
	'Style','text', ...
	'Tag',KnobName);

elseif nargin == 1
   action = varargin{1};
   

   UD = get(gcbf,'UserData');
   eh = UD{1};
	sh = UD{2};
	OV = UD{4};  
	if(length(UD)==5) 
		COMMAND = UD{5};
		% quick fix
		COMMAND = strcat('figure(3);',COMMAND);
	else 
		COMMAND = '';
   end
	

   I = [UD{3}.Position];
	F = [UD{3}.FieldName];
	M = [UD{3}.M];
	N = [UD{3}.N];
	W = [UD{3}.Weight];

   switch lower(action)
   case 'set'
      whoiscalling = gcbo;
      %synchronize slider and text
      switch(whoiscalling)
      case sh 
         newvalue = get(sh,'Value');
         set(eh,'String',newvalue)
      case eh
         newvalue = sscanf(get(eh,'String'),'%f');
         set(sh,'Value',newvalue);
      end
      
      
      %update controlled variables in the worksapce
      for i = 1:length(OV)
         THERING{I(i)}=setfield(THERING{I(i)},F(i),{M(i),N(i)},OV(i)+newvalue*W(i));
      end
	   evalin('base',COMMAND);  
   case 'reset'
      %Do nothing for now
   end
   
end

   


