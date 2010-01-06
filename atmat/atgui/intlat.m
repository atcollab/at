function intlat(varargin)
%INTLAT - interactive AT lattice editor
% INTLAT(DIRECTION)
%  Direction is the initial angle[rad] of the orbit with respect
%  to the plot axis

global THERING
if nargin < 1 | isnumeric(varargin{1})

    if nargin == 1
        STARTANGLE = varargin{1};
    else
        STARTANGLE = 0;
    end
    
    % LAUNCH GUI

	fig = openfig(mfilename,'reuse');
    
    set(fig,'ToolBar','none','HandleVisibility','callback');

    handles = guihandles(fig);
    
    guidata(fig,handles);
    NumElements = length(THERING);
    AllFamNames  = getcellstruct(THERING,'FamName',1:NumElements);
    NumFamilies = 0;
    
    FamNames = {};
    for i=1:NumElements
        if ~any(strcmp(FamNames,THERING{i}.FamName))
            NumFamilies = NumFamilies + 1;
            FamNames{NumFamilies} = THERING{i}.FamName;
        end
    end
    FamNames = sort(FamNames);
    [Families(1:NumFamilies).FamName] = deal(FamNames{:});
    
    set(handles.FamilyPMenu,'String',FamNames);
    set(handles.IconTypePMenu,'String',{'line','rectangle','o','x'});
    
    [x2d, y2d, a2d] = Survey2D(THERING,STARTANGLE);
   	XScale=max(x2d)-min(x2d);
	YScale=max(y2d)-min(y2d);

    set(handles.Axes,'DataAspectRatioMode','manual', ... 
        'DataAspectRatio',[1 1 1],...
        'PlotBoxAspectRatioMode','manual', ...
        'PlotBoxAspectRatio',[XScale YScale 1]);
    
    
    FamNumbers = zeros(1,NumElements);    
    
    for i =1:NumFamilies
        Families(i).KidsList = find(strcmp(Families(i).FamName, AllFamNames));
        Families(i).Display = 1;
        Families(i).Color = [0 0 0];
        Families(i).IconType = 'line';
        Families(i).FieldsList = fieldnames(THERING{Families(i).KidsList(1)});
                
        Families(i).SelectedFields = [find(strcmp(Families(i).FieldsList,'FamName')),...
                find(strcmp(Families(i).FieldsList,'Length')),...
                find(strcmp(Families(i).FieldsList,'PassMethod'))];
                
        Families(i).IconWidth = 0;
        FamNumbers(Families(i).KidsList)=i;
        
    end
    
    Elements = struct('FamIndex',num2cell(FamNumbers),'IconHandle',0);
    Families = SetDefaultIcons(THERING,Families,Elements);
    
    setappdata(fig,'Families',Families);
    setappdata(fig,'Elements',Elements);
    setappdata(fig,'X2D',x2d);
    setappdata(fig,'Y2D',y2d);
    setappdata(fig,'A2D',a2d);

 
        
    ShowFamilyDisplayMode(handles.FamilyPMenu, [], handles);
    
    % Plot all elements
    PlotElements(handles,1:NumElements);
    
        
elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		if (nargout)
			[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
		else
			feval(varargin{:}); % FEVAL switchyard
		end
	catch
		disp(lasterr);
	end

end


function f = SetDefaultIcons(Lattice,Families,Elements)
    L = findspos(Lattice,length(Lattice)+1);
    f = Families;
    % Make default icons for elements of different physical types
    for i=1:length(Families)
        Elem = Lattice{Families(i).KidsList(1)};
        % make icons for bending magnets
        if isfield(Elem,'BendingAngle') & Elem.BendingAngle
            f(i).Display = 1;
            f(i).Color = [1 1 0];
            f(i).IconType = 'rectangle';
            f(i).IconWidth = L/300;
        % Quadrupoles    
        elseif isfield(Elem,'K') & Elem.K
            if Elem.K > 0 % focusing
                f(i).Display = 1;
                f(i).Color = [1 0 0];
                f(i).IconType = 'rectangle';
                f(i).IconWidth = L/400;
            else
                f(i).Display = 1;
                f(i).Color = [0 0 1];
                f(i).IconType = 'rectangle';
                f(i).IconWidth = L/400;
            end
        elseif isfield(Elem,'PolynomB') & length(Elem.PolynomB)>2 & Elem.PolynomB(3)
            if Elem.PolynomB(3)>0 % focusing sextupole
                f(i).Display = 1;
                f(i).Color = [1 0 1];
                f(i).IconType = 'rectangle';
                f(i).IconWidth = L/500;
            else
                f(i).Display = 1;
                f(i).Color = [0 1 0];
                f(i).IconType = 'rectangle';
                f(i).IconWidth = L/500;
            end
        elseif isfield(Elem,'Frequency') & isfield(Elem,'Voltage') % RF cavity
            f(i).Display = 1;
            f(i).Color = [1 0.5 0];
            f(i).IconType = 'o';
            f(i).IconWidth = 0;
            
        end
        
    end
            
        


% --------------------------------------------------------------------
function PlotElements(GUIhandles,INDEX)
    % Retrieve application data
    
    % Temporarily set handle visibility to 'on' for drawing elements
    set(GUIhandles.IntlatMainFigure,'HandleVisibility','on');
    
    figure(GUIhandles.IntlatMainFigure);
    Families = getappdata(GUIhandles.IntlatMainFigure,'Families');
    Elements = getappdata(GUIhandles.IntlatMainFigure,'Elements');
    x2d = getappdata(GUIhandles.IntlatMainFigure,'X2D');
    y2d = getappdata(GUIhandles.IntlatMainFigure,'Y2D');
    a2d = getappdata(GUIhandles.IntlatMainFigure,'A2D');

    
    xcorners = [-1 -1  1  1];
	ycorners = [ 1  1 -1 -1];
    
    
    for i=[INDEX(:)]'
        
        FamIndex = Elements(i).FamIndex;
        
        % If Icon already exists
        if Elements(i).IconHandle
            delete(Elements(i).IconHandle);
            Elements(i).IconHandle = 0;
        end
        
        if  Families(FamIndex).Display
            switch Families(FamIndex).IconType
            case 'rectangle'
                % compute vertex coordinates
                IconWidth = Families(FamIndex).IconWidth;
                vx = [ x2d(i), x2d(i+1), x2d(i+1), x2d(i)] + IconWidth*xcorners*sin((a2d(i)+a2d(i+1))/2);
      	        vy = [ y2d(i), y2d(i+1), y2d(i+1), y2d(i)] + IconWidth*ycorners*cos((a2d(i)+a2d(i+1))/2);
      	        Elements(i).IconHandle = patch(vx,vy,Families(FamIndex).Color);
            case 'line'
                Elements(i).IconHandle = line([x2d(i) x2d(i+1)],[y2d(i) y2d(i+1)]);
                set(Elements(i).IconHandle,'Color',Families(FamIndex).Color);   
            case 'o'
                Elements(i).IconHandle = line([x2d(i) x2d(i+1)],[y2d(i) y2d(i+1)]);
                set(Elements(i).IconHandle,'Color',Families(FamIndex).Color,...
                    'Marker','o','MarkerFaceColor',Families(FamIndex).Color);
            case 'x'
                Elements(i).IconHandle = line([x2d(i) x2d(i+1)],[y2d(i) y2d(i+1)]);
                set(Elements(i).IconHandle,'Color',Families(FamIndex).Color,...
                    'Marker','x');
            end
            % Assign Callback
            set(Elements(i).IconHandle,'UserData',i,...
                'ButtonDownFcn','intlat(''ElementCallback'',gcbo)');
        end
    end
    
    setappdata(GUIhandles.IntlatMainFigure,'Elements',Elements);
    set(GUIhandles.IntlatMainFigure,'HandleVisibility','callback');
    
            
            
            
    
        

function ShowFamilyDisplayMode(h, eventdata, handles, varargin)
    Families = getappdata(handles.IntlatMainFigure,'Families');
    FNum = get(h,'Value');    
    
    set(handles.DisplayCheckBox,'Value',Families(FNum).Display);
    
    
    
    set(handles.ColorSelectionBox,'BackgroundColor',Families(FNum).Color);
    
    set(handles.FieldListBox,'String',Families(FNum).FieldsList);
    set(handles.FieldListBox,'Value',Families(FNum).SelectedFields);
    
    set(handles.WidthEditBox,'String',num2str(Families(FNum).IconWidth));
    
    PossibleIconTypes = get(handles.IconTypePMenu,'String');
    
    set(handles.IconTypePMenu,'value',...
        find(strcmp(Families(FNum).IconType,PossibleIconTypes)));
    
    
    if Families(FNum).Display 
        Visible = 'on';
    else
        Visible = 'off';
    end
    set(handles.IconTypePMenu,'Visible',Visible);
    set(handles.ColorSelectionBox,'Visible',Visible);
    set(handles.WidthEditBox,'Visible',Visible);
    set(handles.FieldListBox,'Visible',Visible);
    set(handles.DisplayIconLabel,'Visible',Visible);
    set(handles.IconColorLabel,'Visible',Visible);
    set(handles.IconWidthLabel,'Visible',Visible);
    set(handles.FieldsLabel,'Visible',Visible);

        
function SelectColor(h, eventdata, handles, varargin)
    Families = getappdata(handles.IntlatMainFigure,'Families');
    FNum = get(handles.FamilyPMenu,'Value');
    NewColor = uisetcolor(['Select icon color for ',Families(FNum).FamName,' family']);
    
    set(handles.ColorSelectionBox,'BackgroundColor',NewColor);
    Families(FNum).Color = NewColor;
    setappdata(handles.IntlatMainFigure,'Families',Families);
    
    PlotElements(handles,Families(FNum).KidsList);

function SetIconType(h, eventdata, handles, varargin)
    Families = getappdata(handles.IntlatMainFigure,'Families');
    FNum = get(handles.FamilyPMenu,'Value');

    PossibleIconTypes = get(handles.IconTypePMenu,'String');
    
    NewIconType = PossibleIconTypes{get(h,'Value')};
    
    Families(FNum).IconType = NewIconType;
    setappdata(handles.IntlatMainFigure,'Families',Families);
    
    
    PlotElements(handles,Families(FNum).KidsList);
    
function SetDisplay(h, eventdata, handles, varargin)
    Families = getappdata(handles.IntlatMainFigure,'Families');
    FNum = get(handles.FamilyPMenu,'Value');

    Families(FNum).Display = get(h,'Value');
    
    
    if Families(FNum).Display 
        Visible = 'on';
    else
        Visible = 'off';
    end
    set(handles.IconTypePMenu,'Visible',Visible);
    set(handles.ColorSelectionBox,'Visible',Visible);
    set(handles.WidthEditBox,'Visible',Visible);
    set(handles.FieldListBox,'Visible',Visible);
    set(handles.DisplayIconLabel,'Visible',Visible);
    set(handles.IconColorLabel,'Visible',Visible);
    set(handles.IconWidthLabel,'Visible',Visible);
    set(handles.FieldsLabel,'Visible',Visible);
    
    setappdata(handles.IntlatMainFigure,'Families',Families);
    PlotElements(handles,Families(FNum).KidsList);
% --------------------------------------------------------------------
function SelectFieldsFromList(h, eventdata, handles, varargin)

    Families = getappdata(handles.IntlatMainFigure,'Families');
    FNum = get(handles.FamilyPMenu,'Value');
    Families(FNum).SelectedFields =  get(h,'Value');
    setappdata(handles.IntlatMainFigure,'Families',Families);
    
    
function SetIconWidth(h, eventdata, handles, varargin)

    Families = getappdata(handles.IntlatMainFigure,'Families');
    FNum = get(handles.FamilyPMenu,'Value');
    NewIconWidth = str2double(get(h,'String'));
    Families(FNum).IconWidth =  NewIconWidth;
    setappdata(handles.IntlatMainFigure,'Families',Families);
    PlotElements(handles,Families(FNum).KidsList);

    
    


% --------------------------------------------------------------------
function varargout = ElementCallback(h)

index = get(h,'UserData');
handles = guidata(gcbo);
Families = getappdata(handles.IntlatMainFigure,'Families');
Elements = getappdata(handles.IntlatMainFigure,'Elements');

FamIndex= Elements(index).FamIndex;
FieldsList =  Families(FamIndex).FieldsList;
Fields2Edit = FieldsList(Families(FamIndex).SelectedFields);

intelem(index,Fields2Edit);

% --------------------------------------------------------------------
function varargout = ColorSelection_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = listbox1_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = Untitled_1_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = ZoomButtonCallback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function [x2d, y2d, a2d] = Survey2D(LATTICE,STARTANGLE)
% Determine 2-d geometry of the LATTICE 
    NumElements = length(LATTICE);
    x2d = zeros(1,NumElements+1);
	y2d = zeros(1,NumElements+1);
	a2d = zeros(1,NumElements+1); % angle of orbit in radians 
	a2d(1) = STARTANGLE;
    for en = 1:NumElements-1
	   if isfield(LATTICE{en},'BendingAngle') 
	      ba = LATTICE{en}.BendingAngle(1); % bending angle in radians
   	   else
          ba = 0;
   	   end
   	 
       if ba == 0
           Lt = LATTICE{en}.Length;
      	   Lp = 0;
       else
      	   Lt = LATTICE{en}.Length*sin(ba)/ba;
      	   Lp = -LATTICE{en}.Length*(1-cos(ba))/ba;
       end
    
       x2d(en+1) = x2d(en) + Lt*cos(a2d(en)) - Lp*sin(a2d(en));
   	   y2d(en+1) = y2d(en) + Lt*sin(a2d(en)) + Lp*cos(a2d(en));
   	   a2d(en+1)=a2d(en) - ba;
        
    end
	x2d(NumElements+1) = x2d(1);
	y2d(NumElements+1) = y2d(1);
	a2d(NumElements+1) = a2d(1);

	X0 = (max(x2d)+min(x2d))/2;
	Y0 = (max(y2d)+min(y2d))/2;	
	x2d = x2d - X0;
	y2d = y2d - Y0;
