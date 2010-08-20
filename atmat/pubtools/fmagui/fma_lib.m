function varargout = fma_lib(varargin)

% This library contains functions needed by FMA_GUI which sets up the
% initial conditions/parameters and provides an interface for the analysis.
% The functions in this library are listed below.
%
%   readfmadata: Reading data file containing FMA related data.
%   writefmadata: Writing a file containing FMA related data.
%   analysis: The routine to calculate the frequency map.
%   freqsearch: Routine to apply modified Fourier Transform to extract
%               the dominant frequency term.
%   plotfm: Plot the results stored in the FMA structure.
%
% COMMANDS:
%    'read'
%    'write'
%    'analysis'
%    'plotfm'
%    'freqsearch'
%    'lineanalysis'
%    'interactive'
%
% USAGE:
%    [fmadatastructure status] = fma_lib('read',filename);
%    success = fma_lib('write',filename,fma);
%    [fmadatastructure status] = fma_lib('analysis',fma,fmagui);
%    freq = fma_lib('freqsearch',data,[RANGE, DT, TOLERANCE, SEARCHWINDOW]);
%    fighandles = fma_lib('plotfm',fma,['file']);
%
% 16/11/2004 Eugene Tan: Created
% 26/11/2009 Eugene Tan: updated for new version of Matlab and MiddleLayer.

if nargin < 1
    error('See help on fma_lib');
end

% Switch yard
switch(varargin{1})
    case 'read'
        if nargin >=2 && ischar(varargin{2})
            [fma status] = readfmadata(varargin{2});
            varargout{1} = fma;
            varargout{2} = status;
        else
            error('Error using ''read'' option. See help on fma_lib');
        end
    case 'write'
        if nargin >= 3 && ischar(varargin{2}) && isstruct(varargin{3})
            varargout{1} = writefmadata(varargin{2},varargin{3});
        else
            error('Error using ''write'' option. See help on fma_lib');
        end 
    case 'analysis'
        try
            [fma status] = analysis(varargin{2:end});
            varargout{1} = fma;
            varargout{2} = status;
            if ~status
                disp(lasterr);
            end
        catch
            disp(lasterr)
            varargout{1} = [];
            varargout{2} = 0;
        end
    case 'freqsearch'
        varargout{1} = local_freqsearch(varargin{2:end});
    case 'plotfm'
        handles = plotfm(varargin{2:end});
        varargout{1} = handles;
    case 'lineanalysis'
        fma_points_around_resonance(varargin{2:end});
    case 'interactive'
        interactive_selection(varargin{2:end});
    otherwise
        disp(['Unknown option: ' varargin{1}]);
end



% ==============
%  Read FMA struct saved in a text file
% ==============
function [fma, status] = readfmadata(filename)

[fin msg] = fopen(filename,'r');
status = 1;

if fin < 0
    fma = [];
    status = 0;
    lasterr(msg);
    return
end

try
    % This is for first declaration to simplify process. Will be removed later.
    fma.temp = NaN;
    while ~feof(fin)
        % Read in variable name
        varname = fgetl(fin);
        % Read in type
        vartype = fgetl(fin);
        % Read number of items to read in
        numvar = str2double(fgetl(fin));

        if numvar == 0
            fgetl(fin);
        else
            if strcmpi(vartype,'char')
                celltemp = cell(1,numvar);
                % Read in the data
                for i=1:numvar
                    celltemp{i} = fgetl(fin);
                end
                temp = char(celltemp);
            else
                temp = zeros(1,numvar);
                % Read in the data
                for i=1:numvar
                    temp(i) = str2num(fgetl(fin));
                end
            end

            eval([varname '=' 'temp;']);
        end

    end
    fma = rmfield(fma,'temp');
catch
    lasterr(['Problem reading from file: ' filename])
    fma = [];
    status = 0;
end

fclose(fin);
% ==============
% End of readfmadata
% ==============


% ==============
%  Write FMA struct to a text file
% ==============
function varargout = writefmadata(filename, fma)

% Check if file already exists. If so, rename appropriately. 
if exist(filename,'file')
    [pathstr namestr ext] = fileparts(filename);
    list = cellstr(ls(sprintf('%s(*',namestr)));
	if isempty(list{1})
        namestr = [namestr '(1)'];
    else
        for i=1:length(list)
            token = regexpi(list{i},'\w\((.)\).dat$','tokens');
            if ~isempty(token)
                num(i) = str2num(token{1}{1});
            end
        end
        namestr = sprintf('%s(%d)',namestr,max(num)+1);
    end
    filename = [pathstr filesep namestr ext];
end
[fout msg] = fopen(filename,'w');
varargout{1} = 1;

if fout < 0
    varargout{1} = 0;
    lasterr(msg);
    return
end

% Print out basic parameters
namelist = {};
namelist = {'version'...
    'latticefile'...
    'resultsfile'...
    'datadir'...
    'starttime'...
    'endtime'...
    'comments'...
    'nturn1'...
    'nturn2'...
    'totaltime'...
    'aperture'};
for i=1:length(namelist)
    fprintf(fout,'fma.%s\n',namelist{i});
    vartype = class(fma.(namelist{i}));
    if strcmpi(vartype,'char')
        % Convert char array into cells
        temp = cellstr(fma.(namelist{i}));
    else
        % Save all other types, typically doubles here.
        temp = fma.(namelist{i});
    end
    % Print type
    fprintf(fout,'%s\n',vartype);
    % Print number/lnegth of variable
    fprintf(fout,'%d\n',numel(temp));
    % Print the variable
    if strcmpi(vartype,'char')
        % Print each cell/string on a different line.
        for j=1:numel(temp)
            fprintf(fout,'%s\n',temp{j});
        end
    else
        % Will automatically print on each line if type is numeric.
        fprintf(fout,'%17.15f\n',temp);
    end
    clear temp;
end

% Mesh
namelist = {};
namelist = {'type'...
    'hor_pos'...
    'maxx_p'...
    'minx_p'...
    'maxdx_p'...
    'mindx_p'...
    'hor_neg'...
    'maxx_n'...
    'minx_n'...
    'maxdx_n'...
    'mindx_n'...
    'ver_pos'...
    'maxy_p'...
    'miny_p'...
    'maxdy_p'...
    'mindy_p'...
    'ver_neg'...
    'maxy_n'...
    'miny_n'...
    'maxdy_n'...
    'mindy_n'...
    'x'...
    'y'};
for i=1:length(namelist)
    fprintf(fout,'fma.mesh.%s\n',namelist{i});
    vartype = class(fma.mesh.(namelist{i}));
    if strcmpi(vartype,'char')
        % Convert char array into cells
        temp = cellstr(fma.mesh.(namelist{i}));
    else
        % Save all other types, typically doubles here.
        temp = fma.mesh.(namelist{i});
    end
    % Print type
    fprintf(fout,'%s\n',vartype);
    % Print number/lnegth of variable
    fprintf(fout,'%d\n',numel(temp));
    % Print the variable
    if strcmpi(vartype,'char')
        % Print each cell/string on a different line.
        for j=1:numel(temp)
            fprintf(fout,'%s\n',temp{j});
        end
    else
        % Will automatically print on each line if type is numeric.
        fprintf(fout,'%17.15f\n',temp);
    end
    clear temp;
end

% Data
namelist = {};
namelist = {'x'...
    'y'...
    'nux1'...
    'nuy1'...
    'nux2'...
    'nuy2'...
    'dindex'...
    'curr'...
    'max'};
for i=1:length(namelist)
    fprintf(fout,'fma.data.%s\n',namelist{i});
    vartype = class(fma.data.(namelist{i}));
    if strcmpi(vartype,'char')
        % Convert char array into cells
        temp = cellstr(fma.data.(namelist{i}));
    else
        % Save all other types, typically doubles here.
        temp = fma.data.(namelist{i});
    end
    % Print type
    fprintf(fout,'%s\n',vartype);
    % Print number/lnegth of variable
    fprintf(fout,'%d\n',numel(temp));
    % Print the variable
    if strcmpi(vartype,'char')
        % Print each cell/string on a different line.
        for j=1:numel(temp)
            fprintf(fout,'%s\n',temp{j});
        end
    else
        % Will automatically print on each line if type is numeric.
        fprintf(fout,'%17.15f\n',temp);
    end
    clear temp;
end

fclose(fout);
% ==============
%  End of writefmadata
% ==============


% ==============
%  Perform FM analysis
% ==============
function [fma, status] = analysis(varargin)

status = 1;

% Parse the input. The first has to be either a struct with the FMA data or
% the filename to the file that contains it.
fma = varargin{1};

boundary = 0;
useTHERING = 0;
if nargin > 1
    for i=nargin:-1:2
        if ischar(varargin{i})
            switch varargin{i}
                case 'boundary'
                    boundary = 1;
                case 'useTHERING'
                    global THERING
                    useTHERING = 1;
            end
        end
    end
end

if ~isstruct(fma) & ischar(fma)
    fprintf('\nOpening FMA input file for reading: %s\n', fma);
    [fma readstatus] = fma_lib('read',fma);
    if ~readstatus
        error(lasterr);
        status = 0;
        return
    end
end

% Check if datadir exists. If so then go to that directory else stay in the
% current directory and generate all files here. It should be in the same
% directory where the input file resides.
if ~isempty(fma.datadir) && exist(fma.datadir,'dir')
    cd(fma.datadir);
    pwd
end

if useTHERING
    if isempty(whos('global','THERING'))
        lasterr('Lattice file not loaded yet. Can''t find THERING.');
        status = 0;
        return
    else
        global THERING;
    end
else
    % Expects the data file to be either in fma.datadir or have the full
    % path specified. Can lead to portability problems if latter option
    % used.
    eval(strtok(fma.latticefile,'.'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEGIN CALCULATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[temp1, resfilename, temp2] = fileparts(fma.resultsfile);

t0 = clock;
fma.starttime = datestr(t0,31);
startindex = fma.data.curr;

% Create an array containing the index points that we want to calculate and
% at which point we want to update progress printouts and output data into
% files.
if boundary
    tocalculate = find(fma.data.x(startindex:end) < fma.mesh.x(3)/1000 |...
                       fma.data.y(startindex:end) < fma.mesh.y(3)/1000);
else
    tocalculate = startindex:fma.data.max;
end

% Save every [updateinc] hours and the next autosave is in [time2update]
% hour(s) time.
updateinc = 0.5; %2.777777777777778e-004; % Update every hour
time2update = 0.5; %2.777777777777778e-004; % Initial update
% Create autosave directory
mkdir autosave;

fprintf('\n********************\n');
fprintf('Beginning calculations... using the following files and parameters:\n');
fprintf('Wrorking DIR: %s\n',fma.datadir);
fprintf('Total particles to track: %d\n',fma.data.max);
fprintf('Number of turns to evaluate: %d %d\n',fma.nturn1, fma.nturn2);
fprintf('Comments: %s\n',fma.comments);
fprintf('Flags -> boundary   = %d\n',boundary);
fprintf('Flags -> useTHERING = %d\n',useTHERING);

switch fma.aperture
    case 'Transverse'
        horiz_ind = 1;
        vert_ind = 3;
    case 'Momentum'
        horiz_ind = 5;
        vert_ind = 1;
end

for i=tocalculate
    fma.data.curr = i;
    
    % Calculate first set of orbits and tunes
    partpos = zeros(6,1);
    partpos(horiz_ind) = fma.data.x(i);
    partpos(vert_ind) = fma.data.y(i);
    if exist('atversion.m','file')
        [T lost] = ringpass(THERING,partpos,fma.nturn1);
    else
        [T lost] = ringpass2(THERING,partpos,fma.nturn1);
    end
    if ~lost
        xorb = complex(T(1,1:end),-T(2,1:end));
        yorb = complex(T(3,1:end),-T(4,1:end));
        fma.data.nux1(i) = local_freqsearch(xorb);
        fma.data.nuy1(i) = local_freqsearch(yorb);

        % DEBUG
        %         fprintf('%15.13f %15.13f', xtune(1,ind), ytune(1,ind));
        %         return

        % Calculate remaining turns for diffusion calculations
        if exist('atversion.m','file')
            [T lost] = ringpass(THERING,T(:,end),fma.nturn2);
        else
            [T lost] = ringpass2(THERING,T(:,end),fma.nturn2); % turn by turn tracking
        end
        if ~lost
            xorb = complex(T(1,1:end),-T(2,1:end));
            yorb = complex(T(3,1:end),-T(4,1:end));
            fma.data.nux2(i) = local_freqsearch(xorb);
            fma.data.nuy2(i) = local_freqsearch(yorb);

            % Calculate the diffusion index
            fma.data.dindex(i) = sqrt((fma.data.nux2(i) - fma.data.nux1(i))^2 +...
                (fma.data.nuy2(i) - fma.data.nuy1(i))^2);
        end
    end

    if (etime(clock,t0)/3600) > time2update
        fma.totaltime = etime(clock,t0)/3600; % In hours
        % Print status in default output window.
        fprintf('\nSimulation started at %s\n',fma.starttime);
        fprintf('Progress printed at %s\n',datestr(clock,31));
        fprintf('=== %5.2f %% completed ===\n', i/fma.data.max*100);
        fprintf('%d out of %d particles computed in %f hours.\n\n', i, fma.data.max,fma.totaltime);
        % Write data to file in case something goes wrong and for progress
        % plots.
        fma.endtime = datestr(clock,31);
        cd autosave;
        writefmadata([resfilename '_' strrep(datestr(now,13),':','_') '.dat'],fma);
        cd ..;
        
        % Next update of progress
        time2update = time2update + updateinc;
    end
    
end

%% Finished!!!
fma.totaltime = etime(clock,t0)/3600;
fma.endtime = datestr(clock,31);
writefmadata(fma.resultsfile,fma);

fprintf('\n\n ************* FINISHED CALCULATIONS *************\n');



function varargout = local_freqsearch(data, varargin)
% =========================================================================
% Find the dominant frequency in the data set.
%   FREQ = LOCAL_FREQSEARCH(DATA, [RANGE, DT, TOLERANCE, SEARCHWINDOW])
%
% RANGE : is the range of the TUNE ie. normalised to 2pi.

MINIT = 6;
MAXIT = 23;
initSigmaLen = 1000;

% Tune range
if nargin >= 2
    range = varargin{1};
else
    range = [0 1];
end

if nargin >= 3
    dt = varargin{2};
else
    dt = 1;
end

if nargin >= 4
    tolerance = varargin{3};
else
    tolerance = 1e-9;
end

% Determines how much to zoom in when narrowing the search.
if nargin >= 5
    searchwindow = varargin{4};
else
    searchwindow = 0.05;
end


% Define the time or running parameter against which to calculate the
% omega.
neval = length(data);
T2 = dt*(neval-1)/2;  %-T/2 --> T/2
t = -T2:dt:T2;

% Will scan this range of frequency.
sigma = [range(1):(range(2)-range(1))/initSigmaLen:range(2)]*2*pi;

% Window function that increases hight of the fundamental peak to make it
% easier to pickout.
p = 1; % cosine window order
kai = complex( 2^p*(factorial(p))^2*(1+cos(pi*(t/T2))).^p/factorial(2*p) );

% This is the spectrum/frequency scan.
psi = zeros(1,length(sigma));
freq = zeros(1,MAXIT);

% Some initial variables.
freq(1) = median(sigma);
difference = 1;
j = 2;

% determine the frequency spectrum
while difference > tolerance || j <= MINIT
    if j >= MAXIT
        break
    end
    % Do the integral that calculates the average <f(t), e^i*sigma*t>. Not
    % including multiplication of some factors like dt, since we only
    % need to find where psi is a maximum and extract the corresponding
    % sigma. Vectorising this loop does not help, evaluated already.
    psi = midpointsum1(data,kai,t,sigma);

    % Plot the frequency spectrum
%     if j >= 2 & j <=2
%         figure; plot(sigma,abs(psi));
%         xlabel('Sigma / Frequency'); ylabel('Arb. Units');
%     end
    
    % Calculate the value of sigma for the maximum psi. Horizontal
    [maxpsi maxind] = max(psi(:));
    % Keep the sign of the maximum sigma found since by using MOD you loose
    % all sign information.
    freq(j) = sigma(maxind);

    difference = abs(freq(j) - freq(j-1));
    
    % If the 'zero' component dominates, remove the component from the
    % function.
%     freq(j) = 0;
    if mod(freq(j),2*pi) ~= 0
        % Find new range to seach, zoom in.
        halfwidth = (sigma(end) - sigma(1))*searchwindow; % fraction of previous range
        if freq(j) == sigma(1)
            lowerbound = freq(j) - halfwidth*2/searchwindow;
            upperbound = freq(j);
        elseif freq(j) == sigma(end)
            lowerbound = freq(j);
            upperbound = freq(j) + halfwidth*2/searchwindow;
        else
            lowerbound = freq(j) - halfwidth;
            upperbound = freq(j) + halfwidth;
        end
%         sigmalen = length(sigma);
%         clear sigma;
        sigma = [lowerbound:(upperbound-lowerbound)/40:upperbound];
        if length(sigma) ~= length(psi)
            % Resize psi since initial sigma is 1000 while subsequent
            % devisions are smaller.
            clear psi
            psi = zeros(1,length(sigma));
        end
    else
        % Orthogonal projection to determine the coeffients. Since
        % e^i*sigma*t where sigma = 1, means e^i*0*t = 1.

%         a = (0.5/T2)*midpointsum(data.*exp(-complex(0,1)*freq(j)*2*pi.*t));
%         e = exp(complex(0,1)*0.3006877345*2*pi*t);
        a = (0.5/T2)*midpointsum(data);
        e = 1;
        
        % Subtract the component from 'f' function.
        data = data - a*e;
    end

    j = j + 1;
end
% fprintf('(%g > %g | %f <= %f) & %f <= %f\n',difference,tolerance,j,MINIT,j,MAXIT);
temp = freq(find(freq ~= 0));
varargout{1} = temp(end)/(2*pi);
% DEBUG
% fprintf('%i    %17.15g\n',j, difference);

% Critical function (computation timewise)
function fctnsum = midpointsum(fctn)

% [a b] = size(fctn);
% if a > 1 && b > 1
%     disp('SIMSUM only takes a one dimensional function');
%     fctnsum = 0;
%     return
% end

% Vectorise the midpoint method of integrating a numerical function.
% f_n      = [fctn 0];
% f_nplus1 = [0 fctn];  % shift all the numbers one "space" to the right.
% midpoints = 0.5*(f_n + f_nplus1);

midpoints = (fctn(1:end-1) + fctn(2:end));
fctnsum = 0.5*sum(midpoints);

% Critical function (computation timewise)
function psi = midpointsum1(data,kai,t,sigma)

len_sigma = 0;
len_sigma = length(sigma);
fctn = zeros(1,length(t));
midpoints = zeros(1,length(t)-1);
temp = data.*kai;

for k=1:len_sigma
    fctn = temp.*complex(cos(sigma(k).*t),-sin(sigma(k).*t));
    midpoints = (fctn(1:end-1) + fctn(2:end));
    psi(k) = 0.5*sum(midpoints);
end




function varargout = plotfm(varargin)

% Defaults/Init
output = 0;
inputstruct = {}; ninput = 0;
tunediagbounds = [13.1  13.4   5.0   5.3];
workingpoint = [13.29 5.226];

for i=1:nargin
    if ischar(varargin{i})
        if strcmpi(varargin{i},'file')
            output = 1;
        else
            % Assumes that the string input is a filename.
            ninput = ninput + 1;
            inputstruct{ninput} = fma_lib('read',varargin{i});
        end
    elseif isstruct(varargin{i})
        ninput = ninput + 1;
        inputstruct{ninput} = varargin{i};
    elseif isnumeric(varargin{i}) && length(varargin{i}) == 4
        tunediagbounds = varargin{i};
    elseif isnumeric(varargin{i}) && length(varargin{i}) == 2
        workingpoint = varargin{i};
    else
        fprintf('\nInput parameter number %d not recognised',i);
    end
end

h_coor = figure;
h_freq = figure;

nrows = floor(sqrt(ninput));
ncols = ceil(ninput/nrows);

for i=1:ninput
    fma = inputstruct{i};
    
    % Log of zero is undefined, so make it a really small number 1e-19.
    fma.data.dindex(find(fma.data.dindex == 0)) = 1e-19;

    % Plot dynamic aperture in coordinate space where the colour scale has
    % been indexed by the diffusion index.
    figure(h_coor); subplot(nrows,ncols,i);
    m = length(fma.mesh.y);
    n = length(fma.mesh.x);
    
    % Check if the definition of the mesh intervals are sorted. Older
    % versions of the code that generated the FMA input files did not sort
    % this and it led to problems when plotting results using SURFACE.
    [val minindx] = min(fma.mesh.x);
    [val maxindx] = max(fma.mesh.x);
    [val minindy] = min(fma.mesh.y);
    [val maxindy] = max(fma.mesh.y);
    if minindx ~= 1 || maxindx ~= n || minindy ~= 1 || maxindy ~= m
        % Sort the mesh intervals
        temp = zeros(m,n);
        temp = log10(reshape(fma.data.dindex,m,n));
        [sortedx xind] = sort(fma.mesh.x);
        [sortedy yind] = sort(fma.mesh.y);
        surface(sortedx,sortedy,temp(yind,xind));
    else
        surface(fma.mesh.x, fma.mesh.y, log10(reshape(fma.data.dindex,m,n)));
    end
    title('Diffusion');
    if isfield(fma,'aperture')
        switch fma.aperture
            case 'Transverse'
                xlabel('Horizontal Displacement (mm)')
                ylabel('Vertical Displacement (mm)')
                axis equal;
            case 'Momentum'
                xlabel('Momentum (%)')
                ylabel('Horizontal Displacement (mm)')
                axis normal;
        end
    else
        xlabel('Horizontal Displacement (mm)')
        ylabel('Vertical Displacement (mm)')
        axis equal;
    end
    set(gca,'XLim',[0.0 1.1*max(fma.mesh.x)]);
    set(gca,'YLim',[0.0 1.1*max(fma.mesh.y)]);
    shading flat; 
    caxis([-12 -0.17]);

    nux = fma.data.nux2+13;
    nux(find(nux > 13.5)) = nux(find(nux > 13.5))-1;
    nuy = fma.data.nuy2+5;
    nuy(find(nuy > 5.5)) = nuy(find(nuy > 5.5))-1;

    figure(h_freq); subplot(nrows,ncols,i);
    plottunediag(7,1, [12.95 13.55 4.95 6.05],workingpoint,'-nospawn');%9,14
    hold on;
    title('');
    scatter(nux, nuy, 5, log10(fma.data.dindex),'filled');
    axis square; 
    caxis([-12 -0.17]);
    axis(tunediagbounds);
end

if output
    print(h_coor,'-depsc2',sprintf('dynamic_aperture_%s.eps',fma.aperture));
end
if output
    print(h_freq,'-depsc2',sprintf('frequency_map_%s.eps',fma.aperture));
end
handles.difhandle = h_coor;
handles.tunehandle= h_freq;

varargout{1} = handles;




function varargout = fma_points_around_resonance(a,b,c,periodicity,maxd,varargin)
% Function will plot only the orbits that have tunes along the resonance
% line defined by a, b and c up to a maximum perpendicular deviation of
% maxd from the line. Can either specify the filename containing the FM
% data or the data struct.

% Resonance line a*nu_x + b*nu_y = c*14 where 14 is the periodicity of the
% ASP lattice.

ninput = 0;
inputstruct = {};
for i=1:nargin-5
    if ischar(varargin{i})
        if strcmpi(varargin{i},'file')
            output = 1;
        else
            % Assumes that the string input is a filename.
            ninput = ninput + 1;
            inputstruct{ninput} = fma_lib('read',varargin{i});
        end
    elseif isstruct(varargin{i})
        ninput = ninput + 1;
        inputstruct{ninput} = varargin{i};
    else
        fprintf('\nInput parameter number %d not recognised',i);
    end
end


for i=1:ninput
    fma = inputstruct{i};
    
    % Log of zero is undefined, so make it a really small number 1e-19.
    fma.data.dindex(fma.data.dindex == 0) = 1e-19;

    ind = []; nind = 0;
    for i=1:length(fma.data.x)
        % Only plot region around resonance line specified above.
        if abs(a*(13 + fma.data.nux2(i)) + b*(5 + fma.data.nuy2(i)) - c*periodicity)/sqrt(a^2 + b^2) < maxd
            nind = nind + 1;
            ind(nind) = i;
        else
            fma.data.dindex(i) = NaN;
        end
    end

    figure; subplot(2,1,1);
    m = length(fma.mesh.y);
    n = length(fma.mesh.x);
    
    % Check if the definition of the mesh intervals are sorted. Older
    % versions of the code that generated the FMA input files did not sort
    % this and it led to problems when plotting results using SURFACE.
    [val minindx] = min(fma.mesh.x);
    [val maxindx] = max(fma.mesh.x);
    [val minindy] = min(fma.mesh.y);
    [val maxindy] = max(fma.mesh.y);
    if minindx ~= 1 || maxindx ~= n || minindy ~= 1 || maxindy ~= m
        % Sort the mesh intervals
        temp = zeros(m,n);
        temp = log10(reshape(fma.data.dindex,m,n));
        [sortedx xind] = sort(fma.mesh.x);
        [sortedy yind] = sort(fma.mesh.y);
        surface(sortedx,sortedy,temp(yind,xind));
    else
        surface(fma.mesh.x, fma.mesh.y, log10(reshape(fma.data.dindex,m,n)));
    end
    title('Diffusion');
    xlabel('Horizontal Displacement (mm)')
    ylabel('Vertical Displacement (mm)')
    axis equal;
    set(gca,'XLim',[0.0 1.1*max(fma.mesh.x)]);
    set(gca,'YLim',[0.0 1.1*max(fma.mesh.y)]);
    shading flat; caxis([-12 -0.17]);
    % print(h,'-depsc2',[pre 'dynamic_aperture.eps'])

    subplot(2,1,2);
    plottunediag(7,14, [12.95 13.55 4.95 6.05],[13.29 5.216],'-nospawn');
    hold on;
    title('');
    scatter(fma.data.nux2+13, fma.data.nuy2+5, 5, log10(fma.data.dindex),'filled');
    axis square; caxis([-12 -0.17]);
    axis([13.1  13.4   5.0   5.3]);
end


function varargout = interactive_selection(varargin)
% Interactively select the points in tunespace and only plot the
% corresponding points in coordinate space.

ninput = 0;
inputstruct = {};
for i=1:nargin
    if ischar(varargin{i})
        if strcmpi(varargin{i},'file')
            output = 1;
        else
            % Assumes that the string input is a filename.
            ninput = ninput + 1;
            inputstruct{ninput} = fma_lib('read',varargin{i});
        end
    elseif isstruct(varargin{i})
        ninput = ninput + 1;
        inputstruct{ninput} = varargin{i};
    else
        fprintf('\nInput parameter number %d not recognised',i);
    end
end

for i=1:ninput
    fma = inputstruct{i};
    
    % Log of zero is undefined, so make it a really small number 1e-19.
    fma.data.dindex(find(fma.data.dindex == 0)) = 1e-19;

    % Pick the points that you want to plot. Lasso is a user written utility
    % currently found at the Matlab file exchange area.
    temph = figure;
    [xtunesel,ytunesel,indsel] = lasso(fma.data.nux2+13, fma.data.nuy2+5);
    delete(temph);
    % Discard; Outside polygon
    tempdindex = zeros(size(fma.data.dindex));
    tempdindex(:) = NaN;
    tempdindex(indsel) = fma.data.dindex(indsel);

    figure; subplot(2,1,1);
    m = length(fma.mesh.y);
    n = length(fma.mesh.x);
        
    % Check if the definition of the mesh intervals are sorted. Older
    % versions of the code that generated the FMA input files did not sort
    % this and it led to problems when plotting results using SURFACE.
    [val minindx] = min(fma.mesh.x);
    [val maxindx] = max(fma.mesh.x);
    [val minindy] = min(fma.mesh.y);
    [val maxindy] = max(fma.mesh.y);
    if minindx ~= 1 || maxindx ~= n || minindy ~= 1 || maxindy ~= m
        % Sort the mesh intervals
        temp = zeros(m,n);
        temp = log10(reshape(tempdindex,m,n));
        [sortedx xind] = sort(fma.mesh.x);
        [sortedy yind] = sort(fma.mesh.y);
        surface(sortedx,sortedy,temp(yind,xind));
    else
        surface(fma.mesh.x, fma.mesh.y, log10(reshape(tempdindex,m,n)));
    end
    title('Diffusion');
    xlabel('Horizontal Displacement (m)')
    ylabel('Vertical Displacement (m)')
    axis equal;
    set(gca,'XLim',[1.1*min(fma.mesh.x) 1.1*max(fma.mesh.x)]);
    set(gca,'YLim',[1.1*min(fma.mesh.y) 1.1*max(fma.mesh.y)]);
    shading flat; caxis([-12 -0.17]);
    % print(h,'-depsc2',[pre 'dynamic_aperture.eps'])

    subplot(2,1,2);
    plottunediag(5,1, [12.95 13.55 4.95 6.05],[13.3 5.2],'nospawn');
    hold on;
    title('');
    scatter(fma.data.nux2+13, fma.data.nuy2+5, 5, log10(tempdindex),'filled');
    axis square; caxis([-12 -0.17]);
    axis([13.1  13.4   5.0   5.3]);
end