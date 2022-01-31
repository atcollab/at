## MATLAB development

To enable the Matlab developer functions, use

```Matlab
>> addpath(fullfile(atroot,'..','Developer','matlab'));
```

### Versioning
A new version should be set on each pull request, using the pull request ID:
```Matlab
>> setversion('major.minor-dev.xxx');
```
where 'xxx' is the pull request ID.
Tagged version should use `'major.minor'` or `'major.minor.patch'`.

The command updates the version number in _atroot_/Contents.m and in _atroot_/at.m.

### Packaging
Tagged versions must be packaged and uploaded on Matlab Central.
```Matlab
>> atclearmex;
```
Deletes all the mex-files (necessary for publishing).

The command `matlab.addons.toolbox.packageToolbox(prj)` generates a wrong path for the toobox.
Instead, use the "package" button in the Matlab toolbox packaging interface.

Matlab versions <= R2021a include useless files in the toolbox. Use R2021b.

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
- Generates the User's guide files for the help browser

#### Custom files
All _atroot_/developer/matlab/mlx/*.mlx must be manually exported as html into _atroot_/../docs/matlab

#### Search database
```Matlab
>> builddocsearchdb(fullfile(atroot,'..','docs','matlab'))
```