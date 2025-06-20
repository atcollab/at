.. _atutils_module:

atutils
=======

.. rubric:: Classes


.. list-table::

   * - :class:`setoption`
     - SETOPTION	Set AT preference values
   * - :class:`getflag`
     - GETFLAG Check the presence of a flag in an argument list
   * - :class:`getoption`
     - GETOPTION Extract a keyword argument from an argument list
   * - :class:`opticsoptions`
     - OPTICSOPTIONS  (private) extract arguments for atlinopt
   * - :class:`parseargs`
     - PARSEARGS Check and expands optional argument lists
   * - :class:`getargs`
     - GETARGS Process positional arguments from the input arguments
   * - :class:`frequency_control`
     - FREQUENCY_CONTROL  Private. Handle off-momentum for 6D lattice
   * - :class:`getdparg`
     - GETDPARG Handle positional dp arguments
   * - :class:`getenvopt`
     -    GETENVOPT(NAME, DEFAULTVALUE)
   * - :class:`wrapper6d`
     - WRAPPER6D  Private. Handle off-momentum for 6D lattice

.. py:class:: atoptions


.. rubric:: Functions


.. list-table::

   * - :func:`setoption`
     - SETOPTION	Set AT preference values
   * - :func:`getflag`
     - GETFLAG Check the presence of a flag in an argument list
   * - :func:`getoption`
     - GETOPTION Extract a keyword argument from an argument list
   * - :func:`opticsoptions`
     - OPTICSOPTIONS  (private) extract arguments for atlinopt
   * - :func:`parseargs`
     - PARSEARGS Check and expands optional argument lists
   * - :func:`getargs`
     - GETARGS Process positional arguments from the input arguments
   * - :func:`frequency_control`
     - FREQUENCY_CONTROL  Private. Handle off-momentum for 6D lattice
   * - :func:`getdparg`
     - GETDPARG Handle positional dp arguments
   * - :func:`getenvopt`
     -    GETENVOPT(NAME, DEFAULTVALUE)
   * - :func:`wrapper6d`
     - WRAPPER6D  Private. Handle off-momentum for 6D lattice

.. py:function:: setoption

   | SETOPTION	Set AT preference values
   | 
   | SETOPTION('KEY',DEFAULT)
   |    Set the default value for the given KEY to DEFAULT. It is an error to set
   |    a default for a non-existing KEY. Use GETOPTION() for the list of
   |    predefined keys.
   | 
   |  KEY:      Key name
   |  DEFAULT:  New default value for the key
   | 
   | SETOPTION('KEY') Resets the default value for KEY to its inital setting
   | 
   | SETOPTION()      Resets all values to their initial setting
   | 
   | See also GETOPTION, ATOPTIONS

.. py:function:: getflag

   | GETFLAG Check the presence of a flag in an argument list
   | 
   | OPTION=GETFLAG(ARGS,OPTNAME)
   |    Return a logical value indicating the presence of the flag name in the
   |    argument list. Flag names are case insensitive.
   | 
   | ARGS:      Argument list (cell array)
   | OPTNAME:	Name of the desired option (string)
   | 
   | [OPTION,NEWARGS]=GETFLAG(ARGS,OPTNAME)
   |            Returns the argument list after removing the processed flag
   | 
   | Example:
   | 
   | function testfunc(varargin)
   | 
   | [optflag,args]=getflag(varargin,'option');     % Extract an optional flag
   | [range,args]=getoption(args,'Range', 1:10);	% Extract a keyword argument
   | [width, height]=getargs(args, 210, 297);       % Extract positional arguments
   | 
   | Dee also GETOPTION, GETARGS

.. py:function:: getoption

   | GETOPTION Extract a keyword argument from an argument list
   | 
   | VALUE=GETOPTION(ARGS,'KEY',DEFAULT)
   | VALUE=GETOPTION(ARGS,KEY=DEFAULT)  in Matlab >= R2021a
   |    Extract a keyword argument, in the form of a pair "key,value" from
   |    input arguments ARGS (typically from VARARGIN).
   |    Return DEFAULT value if the keyword is absent
   | 
   |  ARGS:     Argument list: cell array (usually VARARGIN) or structure
   |  KEY:      Key name
   |  DEFAULT:  Value used if "key,value" is absent from the argument list
   | 
   | VALUE=GETOPTION(ARGS,'KEY')
   |    The default value is taken from a list of predefined keys. Use
   |    GETOPTION() for the list of predefined keys
   | 
   | VALUE=GETOPTION(ARGS,{'KEY1','KEY2',...)
   |    Value is the list of key/value pairs matching KEY1 or KEY2 or...
   | 
   | VALUE=GETOPTION('KEY')
   |    Return the default value of a predefined key. Use GETOPTION() for
   |    the list of predefined keys
   | 
   | VALUE=GETOPTION()
   |    Return all the default values
   | 
   | [VALUE,REMARGS]=GETOPTION(ARGS,...)
   |   Return the remaining arguments after removing the processed ones
   | 
   | Example:
   | 
   | function testfunc(varargin)
   | 
   | [flag,args] = getflag(varargin, 'Flag');       % Extract an optional flag
   | [range,args] = getoption(args, 'Range', 1:10); % Extract a keyword argument
   | [width, height] = getargs(args, 210, 297});    % Extract positional arguments
   | 
   | See also GETFLAG, GETARGS, SETOPTION, ATOPTIONS

.. py:function:: opticsoptions

   | OPTICSOPTIONS  (private) extract arguments for atlinopt
   | 
   |  Separate the options for atlinopt

.. py:function:: parseargs

   | PARSEARGS Check and expands optional argument lists
   | ARGOUT=PARSEARGS(DEFAULT_VALUES,ARGIN)
   | [ARG1,ARG2,...]=PARSEARGS(DEFAULT_VALUES,ARGIN)
   | 
   |  obsolete: see GETARGS

.. py:function:: getargs

   | GETARGS Process positional arguments from the input arguments
   | 
   | Process input arguments (typically from VARARGIN), replacing missing
   | arguments by default values. The default values is also used when an
   | argument is "[]" (empty numeric array).
   | 
   | [V1,V2,...,REMARGS]=GETARGS(ARGIN,DEF1,DEF2,...)
   |    Return as many variables as default values, plus a cell array of
   |    remaining arguments.
   | 
   | [...]=GETARGS(ARGIN,...,'check',checkfun)
   |    Check each input arguments using checkfun(arg) and stops processing
   |    at the first failing argument. Remaining arguments are available in
   |    REMARGS. "checkfun" may either:
   |       - return "false": the function continues using the default value,
   |       - throw an exception: the function aborts, control is returned to
   |         the keyboard or an enclosing try/catch block. The default value
   |         is never used but must be specified.
   | 
   |         Matlab R2020b introduces a series of function "mustBe*" that can
   |         be used for checkfun.
   | 
   | Example 1:
   | 
   | [optflag,args]=getflag(varargin,'option');   % Extract an optional flag
   | [range,args]=getoption(args,'Range', 1:10);  % Extract a keyword argument
   | [dname,dvalue]=getargs(args,'abcd',297);     % Extract positional arguments
   | 
   | Example 2:
   | 
   | global THERING
   | [ring,args]=getargs(varargin,THERING,'check',@iscell)
   |    If the 1st argument is a cell array, it will be used as "ring",
   |    otherwise, THERING will be used. In both cases, the remaining
   |    arguments are available in "args".
   | 
   | Example 3:
   | 
   | function checkcell(arg)
   | if ~iscell(A)
   |     throwAsCaller(MException('AT:WrongType','Argument must be a cell array'));
   | end
   | 
   | [ring,args]=getargs(varargin,{},'check',@checkcell)
   |    If the 1st argument is a cell array, it will be used as "ring" and the
   |    remaining arguments are available in "args". Otherwise, the function
   |    aborts with an error message.
   | 
   | See also GETFLAG, GETOPTION

.. py:function:: frequency_control

   | FREQUENCY_CONTROL  Private. Handle off-momentum for 6D lattice
   | 
   |    VARARGOUT=FREQUENCY_CONTROL(FUNC,RING,VARARGIN)
   | 
   | FUNC   Wrapped function, called as FUNC(RING,VARARGIN{:},'is6d',IS6D)

.. py:function:: getdparg

   | GETDPARG Handle positional dp arguments
   | 
   | [DP,VARARGS]=GETDPARG(VARARGS)
   |    If the 1st argument in VARARGS is a scalar numeric less than 1, it is
   |    considered as DP and removed from VARARGS.
   | 
   | VARARGS=GETDPARG(VARARGS)
   |    DP is extracted, and if it is finite and non-zero,
   |    {'DP', DP} is added to VARARGS

.. py:function:: getenvopt

   |    GETENVOPT(NAME, DEFAULTVALUE)
   |    Looks for an environment variable and return a default value if absent

.. py:function:: wrapper6d

   | WRAPPER6D  Private. Handle off-momentum for 6D lattice
   | 
   |    VARARGOUT=WRAPPER6D(RING,FUNC,VARARGIN)
   | 
   | FUNC   Wrapped function, called as FUNC(RING,IS6D,VARARGIN{:})

