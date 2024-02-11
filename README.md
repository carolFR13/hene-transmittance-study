# Analysis of optical transmittance for FTIR using a HeNe laser

The main goal of this repository is to provide a python package with methods and functions for the analysis of
measurements of the trasmitted and reflected power for a given optical system.

We will work with a system consisting of two triangular prisms coupled as it follows: 

<p align="center">
  <img src="https://github.com/carolFR13/hene-transmission-study/blob/main/data/img/prisms.png" width="450">
</p>

measuring simoultaneously the trasmitted and reflected power with 2 different detectors. Each measurement
was performed with the 2 possible configurations of the detectors ( 1 : 1st detector measuring T and 2nd measuring R and 2 : 1st detector with R and 2nd with T)

The optics package provided here can be installed by following the steps below:

1. Clone this repository on your computer
> git clone https://github.com/carolFR13/hene-transmission-study.git
2. From the new directory install the package from the terminal 
> pip install -e .

Once installed, you can do the following: 

- Read the measurement files and obtain the transmittance and reflectance from the transmitted and reflected power,
  taking into account the constant associated with the detector used.
  > An explanation on how to do the above can be found at `data/doc/hene-measurements.ipynb`
- Analyse the data from the graphs, where we have a region in FTIR and a region outside the same (working as a Fabry-Perot), to calculate
  the angle of the prism ($\alpha$ in the figure) and the distance between the prisms ($d$).
  > See `data/doc/hene-parameters.ipynb` for an explanation of the above.

