# lammps-hrex
This package implements Hamiltonian replica exchange (HREX) into the molecular dynamics code LAMMPS. It also features a well-tempered metadynamics (WTMD) implementation with support for the hybrid WTMD-HREX method.

[Accelerating Solvent Dynamics with Replica Exchange for Improved Free Energy Sampling](https://doi.org/10.1021/acs.jctc.3c00786) <br />
Robert Darkins<sup>1</sup>, Dorothy M. Duffy<sup>1</sup>, Ian J. Ford<sup>1</sup> <br />
<sup>1</sup>University College London

## Installation

* Download the package

```
git clone https://github.com/RDarkins/lammps-hrex.git
```

* Copy the contents of the hrex/ folder into the LAMMPS src/ folder
* Build LAMMPS

## Usage

Please consult [documentation.pdf](documentation.pdf).

## Citation

If this package helps you in your research, please cite:

```
   @article{darkins2023accelerating,
     title={Accelerating Solvent Dynamics with Replica Exchange for Improved Free Energy Sampling},
     author={Darkins, R. and Duffy, D. M. and Ford, I. J.},
     journal={J. Chem. Theory Comput.},
     volume={19},
     number={21},
     pages={7527--7532},
     year={2023},
     publisher={ACS Publications},
     doi={10.1021/acs.jctc.3c00786}
   }
```
