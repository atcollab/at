AT Installation
---------------

The AT code in this repository evolved from version 1.3 of the original code
and is now maintained by the 'atcollab' collaboration, involving people from
different research institutes.

The latest release can be found [on Github](https://github.com/atcollab/at/releases).

Installation process:
---------------------

1. Install git on your computer.

2. Download the latest version of AT:
    ```
    $ git clone https://github.com/atcollab/at.git
    ```

3. Insert recursively the directories `at_installation/atmat` and
`at_installation/atintegrators` in the Matlab path. This can be done by:
    - Using the GUI:
        Open the "Set Path" window, press "Add with subfolders", select
        `at_installation/atmat`; repeat the operation for
        `at_installation/atintegrators`.
    - Using the `startup.m` file:
        Insert a line `addpath(genpath(at_installation/atmat)` and a similar
        one for atintegrators.
    - Temporarily modifying the path by running:
        ```
        >> cd  at_installation/atmat
        >> atpath
        ```

4. Compile all mexfunctions:
    ```
    >> atmexall
    ```

5. Update helpfiles:
    ```
    >> atupdateContents
    ```
    You can now use `athelp` to list all main functions.

6. Update html doc - not yet documented.

### Notes on Octave ###

For the most part AT code works with [Octave](https://www.gnu.org/software/octave/).
At least version 6 is required.

In order to prepare octave environment you'll need to run `atoctave/bootstrap.m` function.
It will install required Octave modules (Internet connection is needed), compile C sources and add
paths.

When in `at` folder, you may start an Octave CLI session with command:
```bash
> octave --eval 'cd atoctave;bootstrap;cd ..;' --persist
```
For GUI session, run command:
```bash
> octave --gui --eval 'cd atoctave;bootstrap;cd ..;' --persist
```
These commands will run `atoctave/bootstrap.m` for you.
