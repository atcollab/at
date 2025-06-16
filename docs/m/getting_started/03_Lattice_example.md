
<a id="T_BA7B5DC2"></a>

# LATTICEDEMO self\-running tutorial
<!-- Begin Toc -->

## Table of Contents
&emsp;[Create elements](#H_4CEDD7DE)
 
&emsp;[Build lattices](#H_246C8121)
 
&emsp;[Summary](#H_51B373F1)
 
<!-- End Toc -->
```matlab
clear
```
<a id="H_4CEDD7DE"></a>

# Create elements

An **element** in Accelerator Toolbox is a 1\-by\-1 MATLAB STRUCTURE.


Functions are provided to easily create such structures with adequate fields.


The following code creates a structure `D1` for a drift space and a structure `QF` for a quadrupole.

```matlab
D1 = atdrift('DR01', 3.0);
QF = atquadrupole('QF', 1.0, 0.2);
```

Use **`whos`**, **`disp`** or just type variable's name without closing semicolon to print the element's info:

```matlab
whos D1 QF
```

```matlabTextOutput
  Name      Size            Bytes  Class     Attributes

  D1        1x1               600  struct              
  QF        1x1              1366  struct              
```

```matlab
disp(D1);
```

```matlabTextOutput
       FamName: 'DR01'
    PassMethod: 'DriftPass'
        Length: 3
         Class: 'Drift'
```

```matlab
QF
```

```matlabTextOutput
QF = struct with fields:
        FamName: 'QF'
     PassMethod: 'StrMPoleSymplectic4Pass'
         Length: 1
          Class: 'Quadrupole'
              K: 0.200000000000000
       PolynomB: [0 0.200000000000000]
       PolynomA: [0 0]
       MaxOrder: 1
    NumIntSteps: 10
```


The next few lines will create another drift structure `D2` from the exiting `D1` and modify the values of fields '`FamName`' and '`Length`'

```matlab
D2 = D1;
D2.FamName = 'DR02';
D2.Length = 2;
```

Create another quadrupole element structure `QD` from `QF` and modify the values of fields '`K`' and '`PolynomB`' to make it defocusing 

```matlab
QD = QF;
QD.FamName = 'QD';
QD.K = -0.4;
```

The field '`PolynomB'` is a vector of coefficients of the polynomial field expansion.


The second element (quadrupole coefficient) must be consistent with field '`K`' 

```matlab
QD.PolynomB(2) = QD.K;
disp(QD);
```

```matlabTextOutput
        FamName: 'QD'
     PassMethod: 'StrMPoleSymplectic4Pass'
         Length: 1
          Class: 'Quadrupole'
              K: -0.400000000000000
       PolynomB: [0 -0.400000000000000]
       PolynomA: [0 0]
       MaxOrder: 1
    NumIntSteps: 10
```


We have created four elements:

```matlab
whos
```

```matlabTextOutput
  Name      Size            Bytes  Class     Attributes

  D1        1x1               600  struct              
  D2        1x1               600  struct              
  QD        1x1              1366  struct              
  QF        1x1              1366  struct              
```

<a id="H_246C8121"></a>

# Build lattices

We are ultimately interested in sequences of elements to model storage ring lattices or single\-pass beam transport lines.


The next section illustrates building such sequences. Accelerator Toolbox represents sequences of elements as MATLAB **cell arrays** where individual cells are 1\-by\-1 structures describing elements.


The following command creates a simple FODO cell by copying previously created element structures into a cell array `fodocell`:

```matlab
fodocell = {QF D1 QD D2 QF};
whos fodocell
```

```matlabTextOutput
  Name          Size            Bytes  Class    Attributes

  fodocell      1x5              5898  cell               
```


**length** is useful to find the number of elements in a sequence:

```matlab
L = length(fodocell) 
```

```matlabTextOutput
L =      5
```


Use the  {:} **cell array** syntax to print some or all elements:

```matlab
fodocell{1}
```

```matlabTextOutput
ans = struct with fields:
        FamName: 'QF'
     PassMethod: 'StrMPoleSymplectic4Pass'
         Length: 1
          Class: 'Quadrupole'
              K: 0.200000000000000
       PolynomB: [0 0.200000000000000]
       PolynomA: [0 0]
       MaxOrder: 1
    NumIntSteps: 10
```

```matlab
fodocell{:}
```

```matlabTextOutput
ans = struct with fields:
        FamName: 'QF'
     PassMethod: 'StrMPoleSymplectic4Pass'
         Length: 1
          Class: 'Quadrupole'
              K: 0.200000000000000
       PolynomB: [0 0.200000000000000]
       PolynomA: [0 0]
       MaxOrder: 1
    NumIntSteps: 10
ans = struct with fields:
       FamName: 'DR01'
    PassMethod: 'DriftPass'
        Length: 3
         Class: 'Drift'
ans = struct with fields:
        FamName: 'QD'
     PassMethod: 'StrMPoleSymplectic4Pass'
         Length: 1
          Class: 'Quadrupole'
              K: -0.400000000000000
       PolynomB: [0 -0.400000000000000]
       PolynomA: [0 0]
       MaxOrder: 1
    NumIntSteps: 10
ans = struct with fields:
       FamName: 'DR02'
    PassMethod: 'DriftPass'
        Length: 2
         Class: 'Drift'
ans = struct with fields:
        FamName: 'QF'
     PassMethod: 'StrMPoleSymplectic4Pass'
         Length: 1
          Class: 'Quadrupole'
              K: 0.200000000000000
       PolynomB: [0 0.200000000000000]
       PolynomA: [0 0]
       MaxOrder: 1
    NumIntSteps: 10
```


Let's build a cell array `fodoring` that represents a closed ring with 10 periods of `fodocell` the same way we would build any other array in MATLAB from the command line:

```matlab
fodoring = [fodocell fodocell fodocell fodocell fodocell...
           fodocell fodocell fodocell fodocell fodocell]; 
whos fodoring
```

```matlabTextOutput
  Name          Size            Bytes  Class    Attributes

  fodoring      1x50            58980  cell               
```


The first element in fodoring is:

```matlab
fodoring{1}
```

```matlabTextOutput
ans = struct with fields:
        FamName: 'QF'
     PassMethod: 'StrMPoleSymplectic4Pass'
         Length: 1
          Class: 'Quadrupole'
              K: 0.200000000000000
       PolynomB: [0 0.200000000000000]
       PolynomA: [0 0]
       MaxOrder: 1
    NumIntSteps: 10
```


To inspect or change the value of a specific field we can use MATLAB syntax for accessing cells in cell arrays and field in structures

```matlab
oldK = fodoring{1}.K
```

```matlabTextOutput
oldK =    0.200000000000000
```

```matlab
fodoring{1}.K = 0.25;
newK = fodoring{1}.K
```

```matlabTextOutput
newK =    0.250000000000000
```


The lattice `fodoring` is a variable in MATLAB workspace.


We can use it in accelerator physics functions and scripts.


For example: the function `findm44` finds 4\-by\-4 transverse transfer matrix

```matlab
M = findm44(fodoring,0)
```

```matlabTextOutput
M = 4x4
  -0.635168208131458  11.030527257692606                   0                   0
  -0.076281158176798  -0.249663952717313                   0                   0
                   0                   0  -0.970624346209193   0.311457414865072
                   0                   0  -0.047116691539558  -1.015145726463917

```

<a id="H_51B373F1"></a>

# Summary
1.  Individual elements are represented by 1\-by\-1  MATLAB structures,
2. Sequences of elements (lattices) are represented by 1\-dimensional MATLAB cell arrays of structures,
3. MATLAB syntax for handling structures and cell arrays applies.

 No special language is required to define a lattice.

