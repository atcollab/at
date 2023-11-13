# Variables and Parameters
Variables and parameters provide a unified way of varying any scalar quantity affecting
a lattice. They may be used in parameter scans, matching, response matrices…

````{grid} 1 1 1 1
:gutter: 2
```{grid-item-card} Variables
:shadow: md

{py:doc}`notebooks/variables` are **references** to any scalar quantity. AT includes two predefined
variable classes referring to scalar attributes of lattice elements:
- an {py:class}`~.element_variables.ElementVariable` is associated to an element object, and acts on
  all occurences of this object. But it will not affect any copy, neither shallow
  nor deep, of the original object,
- a {py:class}`.RefptsVariable` is not associated to an element object, but to an
  element location in a {py:class}`.Lattice`. It acts on any copy of the initial
  lattice. A *ring* argument must be provided to the *set* and *get* methods to
  identify the lattice.

Variable referring to other quantities may be created by:
- deriving the {py:class}`~.variables.Variable` base class. Usually this consist in
  overloading the abstract methods *_setfun* and *_getfun*
- Using the {py:class}`.CustomVariable` class.
```
```{grid-item-card} Parameters

{py:doc}`notebooks/parameters` are objects of class {py:class}`.Param` which can be used instead of numeric
values as {py:class}`.Element` attributes.

Arithmetic combinations of parameters create new read-only parameters of class
{py:class}`.ParamBase`, whose value is permanently kept up-to-date. This is useful to
introduce correlation between attributes of different elements.
```
````

Variables and parameters share a common interface inherited from the
{py:class}`~.variables.Variable` abstract base class. This includes the
{py:meth}`~.variables.Variable.get` and {py:meth}`~.variables.Variable.set` methods
and the {py:attr}`~.variables.Variable.value` property.

Variables and parameters can be grouped in {py:class}`~.variables.VariableList`
containers providing the vectorised equivalent {py:meth}`~.variables.VariableList.get`
and {py:meth}`~.variables.VariableList.set` methods.

```{toctree}
:maxdepth: 2
:hidden:

notebooks/variables.ipynb
notebooks/parameters.ipynb
```