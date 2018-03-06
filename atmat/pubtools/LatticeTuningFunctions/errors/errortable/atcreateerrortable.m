function ErrTab = atcreateerrortable(r)
% ATCREATEERRORTABLE create empty error table for lattice r
%
% USAGE:
% ErrTab = atcreateerrortable(r)
%
% INPUTS:
% r : AT lattice structure
% 
% EXAMPLE:
% r = esrf;
% % create table of errors
% ErrTab = atcreateerrortable(r);
% % enter some errors in the table
% % hor misal
% ErrTab.X = randn(size(ErrTab.X)).*1e-6 ;
% % roll
% ErrTab.Roll = randn(size(ErrTab.X)).*1e-6 ;
%
%see also: atseterrortable

% create table of errors
FamNamesCell = atgetfieldvalues(r,1:length(r),'FamName');
Index = (1:length(r))';
X=zeros(size(r));
Y=zeros(size(r));
S=zeros(size(r));
Roll=zeros(size(r));
Pitch=zeros(size(r));
Yaw = zeros(size(r));
DAngle_Angle = zeros(length(r),1); 
Bend_Roll = zeros(length(r),1); 
DK_K = zeros(length(r),20); 
b_n_systematic = zeros(length(r),20);
b_n_random = zeros(length(r),20);
a_n_systematic=zeros(length(r),20);
a_n_random=zeros(length(r),20);
BPM_Offset = zeros(length(r),2); 
BPM_Gain = zeros(length(r),2); 
BPM_Reading = zeros(length(r),2); 
BPM_Rotation = zeros(size(r)); 

FamNames = categorical(FamNamesCell);
ErrTab = table(...
    FamNames,...
    Index,...
    X,...
    Y,...
    S,...
    Roll,...
    Pitch,...
    Yaw,...
    DAngle_Angle,...
    Bend_Roll,...
    DK_K,...
    b_n_systematic,...
    a_n_systematic,...
    b_n_random,...
    a_n_random,...
    BPM_Offset,...
    BPM_Gain,...
    BPM_Reading,...
    BPM_Rotation);

%ErrTab(1:10,:)

% apply them in the lattice
%rerr = atseterrortable(r,ErrTab);



