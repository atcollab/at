function h0 = intelem(varargin)

%INTELEM - interactive element editor. 
%	
%	INTELEM(INDEX) retrieves THERING{INDEX} from the
%		main workspace and displays the values of all fields for that element
%		Fields that are 1xN vectors or MxN matrixies 
%		such as multipole field data stored in 'PolynomA' are displayed
%		in M raws and N columns, each element in a separate text box.
%
%	INTELEM(INDEX, Fields2Display)
%		Some element models/definitions contain large number of
%		parameters. It may be desired to interactively control only few of them 
%		A cell array of strings Fields2Display allows to select which 
%		element parameters are included in the GUI. 
%		When Fields2Display contains a field name that does not exist for 
% 		an elemet no error is generated ,that field is ignored.
%     For example 
%			Fields2Display = {'FamName' 'Length' 'K' 'BendingAngle'}  
% 			INELEM displays 'FamName' and 'Length' when called for a drift
%			'FamName','Length','K' when called for a quadrupole
%			'FamName','Length','BendingAngle' for a bending magnet etc.
%		
%	INTELEM('action') when the first argument is a string
%	recursively called from inside the INTELEM GUI to evaluate callbacks
%	Possible values for action are
%	'set'
%	'reset'
%	'synch'

global THERING  

if isnumeric(varargin{1})	%initial call
    index = varargin{1};
    UD.LatticeIndex = index;

    ElementRecordCopy = THERING{index};
    if nargin > 1
        NumFields = 0;
        Names = {};
        for i = 1:length(varargin{2})
            if isfield(THERING{index},varargin{2}{i})
                NumFields = NumFields+1;
                Names{NumFields} = varargin{2}{i};
            end
        end
    else
        Names = fieldnames(THERING{index});
        NumFields = length(Names);
    end



    NameBoxWidth = 70;
    NameBoxHeight = 14;

    EditBoxWidth = 60;
    EditBoxWidth2 = 40;
    EditBoxHeight = 14;

    SpaceX =20;
    SpaceY = 15;

    FamilyIndexes = findcells(THERING,'FamName',THERING{index}.FamName);
    KidNum = find(FamilyIndexes == index);
    h0 = figure('Color', [0.8 0.8 0.8], ...
	    'PaperPosition',[18 180 576 432], 'Units','points', 'Position',[30 30 600 200], ...
	    'ToolBar','none','MenuBar','none','NumberTitle','off','Visible','off',...
        'Name',['Lattice Position: ',int2str(index),'      Elemenet # ',int2str(KidNum),...
            '  Element Family: ',THERING{index}.FamName]);

    Handles = cell(1,NumFields);
    TextHandles = zeros(1,NumFields);

    % Create editable text controls for each field
    % If a field is an MxN  matrix (Multipole coefficients) 
    % create MxN text controls for each element of the matrix

    LastPos = 0;
    
    for i = 1:NumFields
        
        FieldData = getfield(THERING{index},Names{NumFields-i+1});
        if ~isempty(FieldData)
            [M,N] = size(FieldData);
            Name = Names{NumFields-i+1};
            UD.FieldName = Name;
            
            LastPos = LastPos + SpaceY  + (M-1)*EditBoxHeight;
            
            % One Static Text control per field 
            TextHandles(i) = uicontrol('Parent',h0, 'Units','points', ...
                'BackgroundColor',[0.8 0.8 0.8], ...
                'FontSize',8, ...
                'FontSize',8 , ...
                'ListboxTop',0, ...
                'Position',[SpaceX  LastPos  NameBoxWidth  NameBoxHeight], ...
                'String',Name, ...
                'HorizontalAlignment','right', ...
                'Style','text', ...
                'Tag','StaticText1');
            
            
            if isnumeric(FieldData)
                for m = 1:M
                    UD.M = m;
                    for n = 1:N
                        UD.N = n;
                        EditHandles{i}(m,n)=uicontrol('Parent',h0, 'Units','points', ...
                            'BackgroundColor',[1 1 1], 'FontSize',8 , ...
                            'Position',[2*SpaceX+NameBoxWidth+(n-1)*EditBoxWidth2 ,  LastPos-(m-1)*EditBoxHeight,  EditBoxWidth2, EditBoxHeight], ...
                            'Style','edit', ...
                            'String',sprintf('%.6f',FieldData(m,n)),'HorizontalAlignment','right', ...      
                            'UserData',UD,...
                            'Callback','intelem sync', ...
                            'Tag','EditText1');
                    end
                end  
            elseif ischar(FieldData)
                UD.M = 1;
                UD.N = 1;
                EditHandles{i}=uicontrol('Parent',h0,'Units','points', ...
                    'BackgroundColor',[1 1 1],'FontSize',8 , ...
                    'Position',[2*SpaceX+NameBoxWidth LastPos  100 EditBoxHeight],'Style','edit', ...
                    'String',FieldData, 'HorizontalAlignment','left', ...
                    'UserData',UD, ...
                    'Callback','intelem sync', ...
                    'Tag','EditText1');
            end
        end
    end

    H = get(h0,'Position');
    H(4) = LastPos+40;
    set(h0,'Position',H);
    set(h0,'HandleVisibility','off','Visible','on');

elseif ischar(varargin{1})

    switch varargin{1} 
    case 'sync'
        UD = get(gcbo,'UserData');
        OldValue = getfield(THERING{UD.LatticeIndex},UD.FieldName);
        if ischar(OldValue)
            THERING{UD.LatticeIndex}=setfield(THERING{UD.LatticeIndex},UD.FieldName,get(gcbo,'String'));
        elseif isnumeric(OldValue)
            st = get(gcbo,'String');
            NewValue = sscanf(st,'%f');
            THERING{UD.LatticeIndex}=setfield(THERING{UD.LatticeIndex},UD.FieldName,{UD.M,UD.N},NewValue);
        end

    end
end