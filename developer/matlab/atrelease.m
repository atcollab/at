function atrelease()
%RELEASE    build, test and package a the AT package
% The release version is assumed to be up to date in the 'Contents.m' file. 
%
%  Copyright 2016-2020 The MathWorks, Inc.

%% Set toolbox name
v = ver(atroot);
shortname = 'atmat';
longname = v.Name;

%% Get AT root directory
rootdir = fullfile(atroot,'..');

%% Check MATLAB and related tools, e.g.:
assert( ~verLessThan( 'MATLAB', '9.6' ), 'MATLAB R2019a or higher is required to use Toolbox Tools.' )

%% Check installation
fprintf( 1, 'Checking installation...' );
switch numel( v )
    case 0
        fprintf( 1, ' failed.\n' );
        error( '%s not found.', shortname );
    case 1
        % OK so far
        version = v.Version;
        fprintf( 1, ' Done.\n' );
    otherwise
        fprintf( 1, ' failed.\n' );
        error( 'There are multiple copies of ''%s'' on the MATLAB path.', shortname );
end

%% Build documentation & examples
fprintf( 1, 'Generating documentation & examples...' );
try
    % Generate custom documentation
    gen_help()
    gen_toc();
    % Remove all mex-files
    atclearmex(atroot);
    fprintf( 1, ' Done.\n' );
catch e
    fprintf( 1, ' failed.\n' );
    e.rethrow()
end

%Build doc search database
try
    builddocsearchdb(fullfile(rootdir,'docs','atdocs','matlab'));
catch me
    warning( me.message )
end

%% Run tests
% fprintf( 1, 'Running tests...' );
% [log, results] = evalc( 'runtests( fullfile( cfdir, "tests" ) )' );
% if ~any( [results.Failed] )
%     fprintf( 1, ' Done.\n' );
% else
%     fprintf( 1, ' failed.\n' );
%     error( '%s', log )
% end
%%
return
%%  Package and rename.
fprintf( 1, 'Packaging...' );
try
    package_name = fullfile(rootdir,'matlab_releases',sprintf('%s.v%s.mltbx',tbxname,version));
    prj = fullfile( rootdir, 'ToolboxPackagingConfiguration.prj');
    matlab.addons.toolbox.packageToolbox( prj);
    fprintf( 1, ' Done.\n' );
catch e
    fprintf( 1, ' failed.\n' );
    e.rethrow()
end

%% Check package
fprintf( 1, 'Checking package...' );
tver = matlab.addons.toolbox.toolboxVersion( package_name );

if strcmp( tver, version )
    fprintf( 1, ' Done.\n' );
else
    fprintf( 1, ' failed.\n' );
    error( 'Package version ''%s'' does not match code version ''%s''.', tver, version )
end

%% Show message
fprintf( 1, 'Created package ''%s''.\n', newMltbx );

end %release
