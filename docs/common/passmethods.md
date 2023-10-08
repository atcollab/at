# AT Passmethods

Passmethods are usually written in C for the best performance, but they may also be
written in python or Matlab for prototyping.

The same passmethods are used by both the Matlab and python interfaces, ensuring a
strict equality of tracking results in both version.

AT comes with a set of passmethods covering most needs, but it is very simple to
add a custom passmethod: 

```{toctree}
:maxdepth: 1

   "Small angle" passmethods <sa_passmethods>
   "Exact" passmethods <ex_passmethods>
```