## MATLAB development

To enable the Matlab developer functions, use

```Matlab
>> addpath(fullfile(atroot,'..','developer','matlab'));
```

### Versioning

#### development version

A new release should be set on each pull request, using the pull request ID:
```Matlab
>> setversion()
```
This will keep the version number unchanged and set the release to the number
of the Pull Request.

#### Release version
```Matlab
>> setversion('major.minor')
```
This sets the version number in _atroot_/Contents.m and in _atroot_/at.m.
The release is set to 'atcollab'

Tagged version should use `'major.minor'` or `'major.minor.patch'`. Extensions are
not allowed by the toolbox packager.

### Packaging
Tagged versions must be packaged and uploaded on Matlab Central.
```Matlab
>> atclearmex;
```
Deletes all the mex-files (necessary for publishing).

The command `matlab.addons.toolbox.packageToolbox(prj)` generates a wrong path for the toobox.
Instead, use the "package" button in the Matlab toolbox packaging interface.

Matlab versions <= R2021a include useless files in the toolbox. Use R2021b for packaging.

### Documentation
All Matlab help files are located in _atroot_/../docs/matlab

Define the chapters by editing `atchapters.m`. Then

```Matlab
gen_help();
```
Generates _atroot_/at.m ("help at")

```Matlab
gen_toc();
```
- Generates _atroot_/developer/matlab/m/*.m files*
- Publishes these files in _atroot_/../docs/matlab/*.html
- Publishes _atroot_/developer/matlab/mlx/*.mlx in _atroot_/../docs/matlab/*.html
- Generates the User's guide files for the help browser

**Warning:**

.mlx files must be run with Matlab 2021a to get plots with the right size. Otherwise,
at the end of the data of each plot, replace `style="height: auto;"` with `style="width: 560;"`

#### Search database
```Matlab
>> builddocsearchdb(fullfile(atroot,'..','docs','matlab'))
```