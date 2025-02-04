# electrochem
This project is used for analysis of electrochemistry data.
That includes importing raw data from MPR files created by ECLab Software as well as subsequent analysis and fitting to models

![Electrochem](https://raw.githubusercontent.com/boblansdorp/electrochem/main/screenshot.png)

Data was included from a Biologic DC3 Dummy Cell and best fit parameters extracted and compared to actual.

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
