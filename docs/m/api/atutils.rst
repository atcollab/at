.. _atutils_module:

atutils
=======

.. rubric:: Classes


.. list-table::

   * - :class:`atoptions`
     - 

.. rubric:: Functions


.. list-table::

   * - :func:`getargs`
     - Process positional arguments from the input arguments
   * - :func:`getdparg`
     - Handle positional dp arguments
   * - :func:`getenvopt`
     - (NAME, DEFAULTVALUE)
   * - :func:`getflag`
     - Check the presence of a flag in an argument list
   * - :func:`getoption`
     - Extract a keyword argument from an argument list
   * - :func:`parseargs`
     - Check and expands optional argument lists
   * - :func:`setoption`
     - Set AT preference values

.. py:function:: atoptions


.. py:function:: getargs(argin,def1,def2,...)

   | Process positional arguments from the input arguments
   
   | Process input arguments (typically from VARARGIN), replacing missing
   | arguments by default values. The default values is also used when an
   | argument is "[]" (empty numeric array).
   
   | **[v1,v2,...,remargs]=getargs(argin,def1,def2,...)**
   |    Return as many variables as default values, plus a cell array of
   |    remaining arguments.
   
   | **[...]=getargs(argin,...,'check',checkfun)**
   |    Check each input arguments using checkfun(arg) and stops processing
   |    at the first failing argument. Remaining arguments are available in
   |    REMARGS. "checkfun" may either:
   |       - return "false": the function continues using the default value,
   |       - throw an exception: the function aborts, control is returned to
   |         the keyboard or an enclosing try/catch block. The default value
   |         is never used but must be specified.
   
   |         Matlab R2020b introduces a series of function "mustBe*" that can
   |         be used for checkfun.
   
   | Example 1:
   
   | [optflag,args]=getflag(varargin,'option');   % Extract an optional flag
   | [range,args]=getoption(args,'Range', 1:10);  % Extract a keyword argument
   | **[dname,dvalue]=getargs(args,'abcd',297)**;     % Extract positional arguments
   
   | Example 2:
   
   | global THERING
   | **[ring,args]=getargs(varargin,thering,'check',@iscell)**
   |    If the 1st argument is a cell array, it will be used as "ring",
   |    otherwise, THERING will be used. In both cases, the remaining
   |    arguments are available in "args".
   
   | Example 3:
   
   | function checkcell(arg)
   | if ~iscell(A)
   |     throwAsCaller(MException('AT:WrongType','Argument must be a cell array'));
   | end
   
   | **[ring,args]=getargs(varargin,{},'check',@checkcell)**
   |    If the 1st argument is a cell array, it will be used as "ring" and the
   |    remaining arguments are available in "args". Otherwise, the function
   |    aborts with an error message.
   
   | See also :func:`getflag`, :func:`getoption`

.. py:function:: getdparg(varargs)

   | Handle positional dp arguments
   
   | **[dp,varargs]=getdparg(varargs)**
   |    If the 1st argument in VARARGS is a scalar numeric less than 1, it is
   |    considered as DP and removed from VARARGS.
   
   | **varargs=getdparg(varargs)**
   |    DP is extracted, and if it is finite and non-zero,
   |    {'DP', DP} is added to VARARGS

.. py:function:: getenvopt

   | (NAME, DEFAULTVALUE)
   |    Looks for an environment variable and return a default value if absent

.. py:function:: getflag(args,optname)

   | Check the presence of a flag in an argument list
   
   | **option=getflag(args,optname)**
   |    Return a logical value indicating the presence of the flag name in the
   |    argument list. Flag names are case insensitive.
   
   | ARGS:      Argument list (cell array)
   | OPTNAME:	Name of the desired option (string)
   
   | **[option,newargs]=getflag(args,optname)**
   |            Returns the argument list after removing the processed flag
   
   | Example:
   
   | function testfunc(varargin)
   
   | **[optflag,args]=getflag(varargin,'option')**;     % Extract an optional flag
   | [range,args]=getoption(args,'Range', 1:10);	% Extract a keyword argument
   | [width, height]=getargs(args, 210, 297);       % Extract positional arguments
   
   | Dee also GETOPTION, GETARGS

.. py:function:: getoption(args,'key',default)

   | Extract a keyword argument from an argument list
   
   | **value=getoption(args,'key',default)**
   | **value=getoption(args,key=default)**  in Matlab >= R2021a
   |    Extract a keyword argument, in the form of a pair "key,value" from
   |    input arguments ARGS (typically from VARARGIN).
   |    Return DEFAULT value if the keyword is absent
   
   |  ARGS:     Argument list: cell array (usually VARARGIN) or structure
   |  KEY:      Key name
   |  DEFAULT:  Value used if "key,value" is absent from the argument list
   
   | **value=getoption(args,'key')**
   |    The default value is taken from a list of predefined keys. Use
   |    **getoption()** for the list of predefined keys
   
   | **value=getoption(args,{'key1','key2',...)**
   |    Value is the list of key/value pairs matching KEY1 or KEY2 or...
   
   | **value=getoption('key')**
   |    Return the default value of a predefined key. Use **getoption()** for
   |    the list of predefined keys
   
   | **value=getoption()**
   |    Return all the default values
   
   | **[value,remargs]=getoption(args,...)**
   |   Return the remaining arguments after removing the processed ones
   
   | Example:
   
   | function testfunc(varargin)
   
   | [flag,args] = getflag(varargin, 'Flag');       % Extract an optional flag
   | **[range,args] = getoption(args, 'range', 1:10)**; % Extract a keyword argument
   | [width, height] = getargs(args, 210, 297});    % Extract positional arguments
   
   | See also :func:`getflag`, :func:`getargs`, :func:`setoption`, :func:`atoptions`

.. py:function:: parseargs(default_values,argin)

   | Check and expands optional argument lists
   | **argout=parseargs(default_values,argin)**
   | **[arg1,arg2,...]=parseargs(default_values,argin)**
   
   |  obsolete: see GETARGS

.. py:function:: setoption('key',default)

   | Set AT preference values
   
   | **setoption('key',default)**
   |    Set the default value for the given KEY to DEFAULT. It is an error to set
   |    a default for a non-existing KEY. Use GETOPTION() for the list of
   |    predefined keys.
   
   |  KEY:      Key name
   |  DEFAULT:  New default value for the key
   
   | **setoption('key')** Resets the default value for KEY to its inital setting
   
   | **setoption()**      Resets all values to their initial setting
   
   | See also :func:`getoption`, :func:`atoptions`

