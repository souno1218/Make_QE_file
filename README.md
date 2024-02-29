# Make_QE_file

Quantum Espresso calculations first require an input file.   
Since it is time-consuming to create such a file by copy-and-paste from a crystal information file (CiF), etc., we have scripted it in python.

- cif2qe_in
  - Create QE input file from CiF (Only relax and vc-relax are currently confirmed.)
- qe_out2in
  - Create QE input file from QE output file (Only relax and vc-relax are currently confirmed.)
  - Used for relax->vc-relax, etc.
- qe_out2cif
  - Create CiF from QE output file (Only relax and vc-relax are currently confirmed.)

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

pymatgen,numpy,pandas will be installed automatically when not present

## Running

Explain how to run the automated tests for this system



## Built With

* [Pymatgen](https://pymatgen.org) - Used to read CiF
* [numpy](https://numpy.org) - Used for various calculations
* [pandas](https://pandas.pydata.org) - Used to read csv



## Authors

* **河野 颯之介(Sonosuke Kono)**

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Finally.

I am Japanese and had never used GitHub until I wrote this.   
I use Deepl because I am not good at English.   
This ReedMe is also written with reference to the following page.   
https://gist.github.com/PurpleBooth/109311bb0361f32d87a2

I would like to ask you to support us continuously.
