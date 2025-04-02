# Make_QE_file
This is a collection of Python scripts that partially automate the creation of input files required for Quantum Espresso calculations.   
The tool provides the following functionalities:   
 - Generate input files from CIF files
 - Read structures from relax and vc-relax output files and create the next input file
 - Read structures from relax and vc-relax output files and convert them to CIF format
 - Plot band structures, PDOS, and other related data

The intended use case is when a calculation has already converged for one material, and there are many similar materials for which you want to run calculations under the same conditions.   
Although it can also be used for the initial calculation of a single material, the script is primarily designed to reuse successful conditions efficiently.   
I was able to use it for my application, but I think there are cases where it may not work.   
Please modify it by yourself.   

## Getting Started

### Prerequisites

Requires pymatgen,numpy,pandas. If not, installation is automatic.

### Installing

First, activate the virtual environment if it is separated by conda.

```
#examples
conda activate myenv
```

Download and Installation

```
 pip install git+https://github.com/souno1218/Make_QE_file.git
```

pymatgen, numpy, pandas, matplotlib will be installed automatically when not present

## Running

Explain how to run the automated tests for this system



## Built With

* [moyopy](https://spglib.github.io/moyo/python/index.html) - Used to obtain the Hermann-Mauguin notation from the space group number.
* [numpy](https://numpy.org) - Used for various calculations.
* [pandas](https://pandas.pydata.org) - Used to .
* [matplotlib](https://matplotlib.org) - Used to plot.



## Authors

* **河野 颯之介(Sonosuke Kono)**

## License

This project is licensed under Apache License, Version 2.0 - see the [LICENSE](LICENSE) file for details

## Finally.

I am Japanese and had never used GitHub until I wrote this.   
I use Deepl because I am not good at English.   
This ReedMe is also written with reference to the following page.   
https://gist.github.com/PurpleBooth/109311bb0361f32d87a2

I would like to ask you to support us continuously.
