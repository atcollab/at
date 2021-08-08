%% Lattice manipulation
% 
%%
% <matlab:doc('isatlattice') isatlattice> - Tests if an input argument is a valid AT lattice
%
%% Access elements
% 
%%
% <matlab:doc('atindex') atindex> - Extracts the information about element families and
%
% <matlab:doc('atgetcells') atgetcells> - Performs a search on MATLAB cell arrays of structures
%
% <matlab:doc('atgetfieldvalues') atgetfieldvalues> - Retrieves the field values AT cell array of elements
%
% <matlab:doc('atsetfieldvalues') atsetfieldvalues> - Sets the field values of MATLAB cell array of structures
%
%% Insert elements
% 
%%
% <matlab:doc('atinsertelems') atinsertelems> - Insert elements at given locations in a line
%
% <matlab:doc('atdivelem') atdivelem> - LINE=ATDIVELEM(ELEM,FRAC) divide an element into pieces
%
% <matlab:doc('atsplitelem') atsplitelem> - Creates a line by inserting one or more elements into a base element
%
% <matlab:doc('insertindrift') insertindrift> - Inserts one or more elements into a drift element
%
% <matlab:doc('atsbreak') atsbreak> - Insert markers at given s positions in a lattice
%
%% Join elements
% 
%%
% <matlab:doc('atreduce') atreduce> - Remove useless elements from an AT structure
%
% <matlab:doc('mergedrift') mergedrift> - Removes a lattice element and merges the two adjacent drift spaces
%
% <matlab:doc('combinebypassmethod') combinebypassmethod> - Combines adjacent elements that have the same specified pass method
%
% <matlab:doc('combinelinear45') combinelinear45> - Combines adjacent  elements that use 4-by-5 PassMethods
%
%% Other
% 
%%
% <matlab:doc('atloadfielderrs') atloadfielderrs> - Will load a field error structure into a ring
%
% <matlab:doc('atsetRFCavity') atsetRFCavity> - - sets the RF Cavity with the passmethod RFCavityPass
%
% <matlab:doc('atsetshift') atsetshift> - Sets the misalignment vectors
%
% <matlab:doc('atsettilt') atsettilt> - Sets the entrance and exit rotation matrices
%
% <matlab:doc('settags') settags> - Sets the 'Tag' field in AT lattice elements
%
% <matlab:doc('findtags') findtags> - Looks for string matches in 'Tag' field of AT lattice elements
%
% <matlab:doc('mvelem') mvelem> - (ELEMPOS, DIST) moves an element  located at ELEMPOS in THERING
%
% <matlab:doc('mvfield') mvfield> - Move fields from one structure to another
%
% <matlab:doc('reverse') reverse> - Reverses the order of elements in a one-dimensional MATLAB ARRAY
%
% <matlab:doc('splitdrift') splitdrift> - Inserts an element into a drift space
%
