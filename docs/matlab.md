---
title: Matlab
layout: page
permalink: /matlab.html
sidebar: matlab_sidebar
---
Automatic HTML documentation is generated from mfile using the [m2html package](https://www.artefact.tk/software/matlab/m2html/) by Guillaume Flandin.

**[Online documentation of all functions](https://cdn.rawgit.com/atcollab/atdoc/aa9b9f58/doc_html/index.html)**

## Requirements:
AT is compatible with:
- **Matlab** release >= R2016b (Matlab 9.1)
- **Octave** >= 6.0

## Installation:
### From Matlab Central
Available soon…
### From GitHub
1. Install git on your computer.

2. Download the latest version of AT:
    ```shell
    $ git clone https://github.com/atcollab/at.git
    ```

3. Insert recursively the directories `<at_installation>/atmat` and
`<at_installation>/atintegrators` in the Matlab path. This can be done by:
    - Using the GUI:
        Open the "Set Path" window, press "Add with subfolders", select
        `<at_installation>/atmat`; repeat the operation for
        `<at_installation>/atintegrators`.
    - Using the `startup.m` file:
        Insert a line `addpath(genpath(<at_installation>/atmat))` and a similar
        one for atintegrators.
    - Temporarily modifying the path by running:
      ```matlab
      >> cd at_installation/atmat
      >> atpath
      ```

4. Compile all mexfunctions:
   ```matlab
   >> atmexall
   ```

### Notes on Octave ###

For the most part AT code works with [Octave](https://www.gnu.org/software/octave/).
At least version 6 is required.

In order to prepare octave environment you'll need to run `atoctave/bootstrap.m` function.
It will install required Octave modules (Internet connection is needed), compile C sources and add
paths.

When in `at` folder, you may start an Octave CLI session with command:
```shell
$ octave --eval 'cd atoctave;bootstrap;cd ..;' --persist
```
For GUI session, run command:
```shell
$ octave --gui --eval 'cd atoctave;bootstrap;cd ..;' --persist
```
These commands will run `atoctave/bootstrap.m` for you.
