.. _atmat_module:

atmat
=====

.. rubric:: Functions


.. list-table::

   * - :func:`atcavityoff`
     - switches RF cavities off
   * - :func:`atcavityon`
     - ATRADON switches RF cavities on
   * - :func:`atdiag`
     - Tests AT intallation
   * - :func:`atdisplay`
     - checks the verbosity level in the global variable GLOBVAL
   * - :func:`athelp`
     - generate the list of Accelerator Toolbox functions
   * - :func:`atm2html`
     - MAKEDOC_HTML - Generate new MML, SOLEIL and AT HTML help files
   * - :func:`atmexall`
     - build all AT platform dependent mex-files from C-sources
   * - :func:`atpath`
     - Adds the AT directories to the MATLAB search path
   * - :func:`atroot`
     - returns Accelerator Toolbox root directory
   * - :func:`getContents`
     - Get the contents of a specified directory
   * - :func:`isOctave`
     - isOctave check if running Octave
   * - :func:`updateContents`
     - Create a Contents.m file including subdirectories

.. py:function:: atcavityoff(ring)

   | switches RF cavities off
   | 
   |   **[ring2, cavitiesindex] = atcavityoff(ring)**
   |     Changes cavity passmethods to turn off acceleration
   | 
   |   INPUTS:
   |   1. RING      initial AT structure
   | 
   |   OUPUTS
   |   1. RING2          output ring with cavities switched off
   |   2. CAVITIESINDEX  indices cavities
   | 
   |   See also ATCAVITYON, ATRADON, ATRADOFF

.. py:function:: atcavityon(ring,cavitypass)

   | ATRADON switches RF cavities on
   | 
   |   **[ring2,cavindex]=atcavityon(ring,cavitypass)**
   |     Changes cavity passmethods to get RF acceleration
   | 
   |   INPUTS
   |   1. RING	     	initial AT structure
   |   2. CAVITYPASS	customed passmethod for cavities (default RFCavityPass)
   | 
   |   OUPUTS
   |   1. RING2          output ring with cavities off
   |   2. CAVITIESINDEX  indices of cavities
   | 
   |   See also ATCAVITYOFF, ATRADON, ATRADOFF

.. py:function:: atdiag

   | Tests AT intallation

.. py:function:: atdisplay

   | checks the verbosity level in the global variable GLOBVAL
   |           and displays message if this is greater than the verbosity
   |           for this message.
   | 
   |   See also NUMDIFPARAMS

.. py:function:: athelp

   | generate the list of Accelerator Toolbox functions
   | 
   |   INPUTS
   |     No argument - open the help file in Matlab browser
   |     'new'       - force the update of the documentation, which requires a few more seconds
   | 
   |   EXAMPLES
   |   1. **athelp**: full help.
   |   2. for selected help, use help directory where directory is
   |       help atintegrators
   |       help atmat
   |       help atdemos
   |       help atgui
   |       help atmatch
   |       help atphysics
   |       help linearoptics
   |       help longitudinaldynamics
   |       help nonlineardynamics
   |       help atplot
   |       help plotfunctions
   |       help lattice
   |       help element_creation
   |       help pubtools
   |       help survey
   |       help lattice_tools
   |       help LatticeTuningFunctions
   |       help machine_date
   |       help tuneandchromaticity
   |       help touschekpiwinski
   |       help radiation
   |       help parametersummaryfunctions
   | 
   | 
   |   See also help

.. py:function:: atm2html

   | MAKEDOC_HTML - Generate new MML, SOLEIL and AT HTML help files
   |   makedoc_html
   | 
   |   HOWTO
   |   1. Make sure to update and run toolboxUpdateHeader.m
   |   2. Update history.txt appropriately, including w current version
   |   3. Update overview.html file with the version/date/link to zip:
   |      edit external/m2html/templates/at/about.html
   |   4. Need to install graphviz fro graph dependency
   |      see: https://graphviz.org/

.. py:function:: atmexall

   | build all AT platform dependent mex-files from C-sources
   | 
   | **atmexall** option1 ... optionN
   | 
   |  AT Options:
   | 
   | 	-missing    Build only the outdated components
   |    -fail       Throw an exception if compiling any passmethod fails
   |                (By defaults compilation goes on)
   | 	-openmp     Build the integrators for OpenMP parallelisation
   | 	-cuda CUDA_PATH Build the GPU tracking support using Cuda
   | 	-opencl OCL_PATH Build the GPU tracking support using OpenCL
   |                Use "-opencl default" for using standard OpenCL install
   |    -c_only     Do no compile C++ passmethods
   |    -DOMP_PARTICLE_THRESHOLD=n
   |                Set the parallelisation threshold to n particles
   |                (Default 10)
   | 
   |  Options forwarded to the mex command:
   | 
   |    -v          Verbose output
   |    -g          Compile with debug options
   |    -O          Optimize the object code (Default)
   |    -n          Display the generated command without executing
   |    ...
   | 

.. py:function:: atpath

   | Adds the AT directories to the MATLAB search path

.. py:function:: atroot

   | returns Accelerator Toolbox root directory

.. py:function:: getContents(directory)

   | Get the contents of a specified directory
   | 
   |    This function returns the contents of a specified directory.
   | 
   |    CONT = IOSR.GENERAL.**getcontents(directory)** returns the files and
   |    folders in a directory and returns them to the cell array cont. It
   |    ignores hidden files and folders (those starting '.'). DIRECTORY must
   |    be a character array (string).
   | 
   |    CONT = IOSR.GENERAL.**getcontents(directory,'parameter',value)** allows
   |    search options to be specified. The options include:
   |        'rec'       {false} | true
   |                    Search recursively within the subfolders of the
   |                    specified directory.
   |        'path'      {'relative'} | 'full'
   |                    Specifies whether returned paths are full or relative
   |                    to the specified directory.
   |        'sort'      {false} | true
   |                    Specify whether the output is sorted alphabetically.
   |        'filter'    {'all'} | 'files' | 'folders' | '*.ext' | str
   |                    This option allows a filter to be specified. 'files'
   |                    returns names of all files in the directory. 'folders'
   |                    returns names of all folders in the directory. '*.ext',
   |                    where 'ext' is a user-specified file extension, returns
   |                    all files with the extension '.ext'. str may be any
   |                    string; only elements that contain str will be returned
   |                    (files or folders). str is case-sensitive.
   | 
   |    [CONT,DIRFLAG] = IOSR.GENERAL.**getcontents(...)** returns a logical array
   |    DIRFLAG, the same size as CONT, indicating whether each element is a
   |    directory.
   | 
   |    Examples
   | 
   |        Ex. 1
   | 
   |        % Return all m-files in the current directory
   | 
   |        cont = iosr.general.**getcontents(cd,'filter','*.m')**
   | 
   |        Ex. 2
   | 
   |        % Return all files in the current directory and its
   |        % sub-directories
   | 
   |        cont = iosr.general.**getcontents(cd,'rec',true)**
   | 
   |        Ex. 3
   | 
   |        % Return all files in current directory with names
   |        % containing 'foo'
   | 
   |        % may return files and folders:
   |        [cont,dirflag] = iosr.general.**getcontents(cd,'filter','foo')**
   | 
   |        % use dirflag to limit:
   |        cont = cont(~dirflag);

.. py:function:: isOctave()

   | isOctave check if running Octave
   | 
   |   **[retval]=isoctave()**
   |     Check if running Octave
   | 
   |   OUTPUTS
   |   1. RETVAL          boolean is running Octave
   | 

.. py:function:: updateContents(folder)

   | Create a Contents.m file including subdirectories
   | 
   |    **updatecontents** scans through the current directory, and
   |    its subdirectories, and builds a Contents file similar to Matlab's
   |    report-generated Contents.m files. Any existing Contents.m file will be
   |    overwritten.
   | 
   |    **updatecontents(folder)** scans through the directory FOLDER.
   | 
   |    Typing
   |        help(FOLDER)
   |    or
   |        help path/to/folder
   | 
   |    will display Contents.m in the Command Window, and display links to the
   |    help for any functions that are in Matlab's search path.
   | 
   |    NB: Do not use Matlab's Contents Report generator to edit the
   |    Contents.m file. Execute this function to update it.

