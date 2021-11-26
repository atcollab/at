## MATLAB help files

All Matlab help files are located in atroot/../docs/matlab

### Set a new version and release number
in atroot/Contents.m

### Automatically generated files
Define the chapters by editing atchapters.m

gen_help:
Generate _atroot_/at.m ("help at")

gen_toc:
Generate _atroot_/developer/mdoc/m/*.m files*
Publish these files in _atroot_/../docs/matlab/*.html
Generate the User's guide files for the help browser

### Custom files
All _atroot_/developer/mdoc/mlx/*.mlx must be manually exported as html into _atroot_/../docs/matlab

### Search database
builddocsearchdb(_atroot_/../docs/matlab)