# electrochem
This project is used for analysis of electrochemistry data.
That includes importing raw data from MPR files created by ECLab Software as well as subsequent analysis and fitting to models

![Electrochem](https://raw.githubusercontent.com/boblansdorp/electrochem/main/screenshot.png)

Data was included from a Biologic DC3 Dummy Cell and best fit parameters extracted and compared to actual.

## Getting Started
The easiest way to get started is to download the files and then use Matlab to run:
```Matlab
runtests('tests')
```

One of the tests is testZFit. This test will take a EIS DC3.mpr file in the testData directory, convert it into Matlab usable data, and use ZFit to find the best fit paramaters. 
Since it's a Biologic DC3 Dummy cell, we know that is should be            
```Matlab
nominal_p = [499, 1000, 10e-9, 3.57e3, 2.2e-6];
```
which is a 499 Ω resistor in series with:
a parallel combo of (1000 Ω resistor and 10 nF capacitor) which is in series with
a 3.57 kΩ resistor in parallel with a 2.2 μF capacitor
we can see the definition of the parallel p and series parts of the circuit s in this part of the test: 
```
circuit='s(R1,s(p(R1,C1),p(R1,C1)))';  
```
We need to give it some starting guesses for the resistors and capacitors:
```
param=[100, 1000, 1e-9,1e3, 1e-6];  
```
and a lower bound for everything so it doesn't return negative resistors or anything weird:
 ```
LB=[0,0,0,0,0];
```

it should show:
```
Best Fit Parameter values vs Expected:
Fit(1) = 4.990e+02 Ω, Expected(1) = 4.990e+02 Ω 
Fit(2) = 9.962e+02 Ω, Expected(2) = 1.000e+03 Ω 
Fit(3) = 9.059e-09 F, Expected(3) = 1.000e-08 F 
Fit(4) = 3.575e+03 Ω, Expected(4) = 3.570e+03 Ω 
Fit(5) = 2.249e-06 F, Expected(5) = 2.200e-06 F 
```
the resulting best-fit resistor and capacitor values should be within 10% of the values expectd based on the labels on the DC3 cell. Neat!


## Licensing

This repository contains MATLAB scripts with different licensing terms:

- **Most of the code** (`src/`) is released under the **MIT License**.
- **The file `third_party/parseBiologicMPR.m` is derived from a GPL-licensed Python project** and is licensed under **GPLv3**.
- To comply with GPL, `parseBiologicMPR.m` must be distributed under the same license.
- **The file `third_party/Zfit.m` is derived from a custom license 

For details, see:
- [`LICENSE`](./LICENSE) for MIT-licensed files.
- [`third_party/LICENSE_GPL.txt`](./third_party/LICENSE_GPL.txt) for GPL-covered code.
- [`third_party/license_ZFit.txt`](./third_party/LICENSE_GPL.txt)
