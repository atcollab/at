function varargout = fma_calc(fma,handles)

disp('In fma_calc and updating a gui handle/object');
set(handles.comments,'String','This has been set from fma_calc');

varargout{1} = 1;