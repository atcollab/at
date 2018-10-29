AT Installation
---------------

The AT code in this repository evolved from version 1.3 of the original code
and is now maintained by the 'atcollab' collaboration, involving people from
different research institutes.

The latest release can be found [on Github](https://github.com/atcollab/at/releases).

Installation process:

1. Install git on your computer.

2. Download the latest version of AT:
    `$ git clone https://github.com/atcollab/at.git`

3. Insert recursively the directories at_installation/atmat and at_installation/atintegrators
in the Matlab path. This can be done by:
    - Using the GUI:
        Open the "Set Path" window, press "Add with subfolders", select
        at_installation/atmat; repeat the operation for
        at_installation/atintegrators.
    - Using the startup.m file:
        Insert a line addpath(genpath(at_installation/atmat) and a similar one
        for atintegrators.
    - Temporarily modifying the path by running:
        `>> cd  at_installation/atmat`
        `>> atpath`

4. Go in atmat directory:
    `>> cd at/atmat`

5. Compile all mexfunctions:
    `>> atmexall`

6. Update helpfiles:
    `>> atupdateContents`
    You can now use `athelp` to list all main functions.

7. Update html doc - not yet documented.
