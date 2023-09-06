# plind.py Documentation

## Introduction
plind.py (Picard-Lefschetz Integration in N-Dimensions) is a Python package that implements the Picard-Lefschetz method for integrating highly oscillatory functions in any number of dimensions. 

## Description
This code is for people interested in computing integrals of the form
$$I=\int_\Omega \text{d}^n\mathbf{x}\exp i S(\mathbf{x}\;\mathbf{\mu}).$$ 
These integrals occur in quantum theory and wave optics. 
Such integrals are not possible to do using traditional numerical methods due to their highly oscillatory nature. 

However, Picard-Lefschetz (PL) theory gives a recipe for computing them by analytically continuing the integration domain to $\mathbb{C}^n$, deforming the integration domain according to certain rules, then integrating on the new domain. 
This package is an implementation of the PL algorithm in Python. 

For details, see papers in the [additional information](#additional-information) section. 

## Installation
To use plind.py, run `pip install plind`.

## Usage
To use the code, follow the examples in `Example Notebook.ipynb`. 

## Additional Information
For more information on the code and its usage, refer to the docstrings and comments within the code files, or email the authors. See also [arXiv:2103.08687](https://arxiv.org/abs/2103.08687) and [arXiv:1909.04632](https://arxiv.org/abs/1909.04632) for details on the method.

## Citation
Please cite [doi:10.1093/mnras/stab2337](https://doi.org/10.1093/mnras/stab2337) if you use this code in your project.
