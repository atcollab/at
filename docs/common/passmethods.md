# AT Passmethods

The tracking through each AT element is performed by an integrator, called its
"passmethod". There is no correlation between an element class and its associated
passmethod. The passmethod must be explicitly specified by the `PassMethod` attribute 
of the element.

Passmethods are usually written in C for the best performance, but thys may also be
written in python or Matlab for prototyping.

```{toctree}
:maxdepth: 1

   "Small angle" passmethods <sa_passmethods>
   "Exact" passmethods <ex_passmethods>
```