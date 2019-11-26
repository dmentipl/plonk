---
title: 'Plonk: Smoothed particle hydrodynamics analysis and visualization with Python'
tags:
  - Python
  - astronomy
  - smoothed particle hydrodynamics
  - visualization
authors:
  - name: Daniel Mentiplay
    orcid: 0000-0002-5526-8798
    affiliation: 1
affiliations:
 - name: School of Physics and Astronomy, Monash University
   index: 1
date: 08 November 2019
bibliography: paper.bib
---

# Summary

``Plonk`` is a Python package for analysis and visualization of smoothed particle hydrodynamics data focussing on astrophysical applications.

Hydrodynamical simulations are a standard part of the astronomer’s toolkit. The knowledge required to run these simulations and analyze the output can be specialized. Many of the codes written to run and analyze these simulations are written in compiled languages such as Fortran and C++. The advantage of using compiled languages is performance. One disadvantage is that expertise in these languages is becoming less widely available than other, typically interpreted, languages such as Python. Access to analysis and visualization of hydrodynamical simulations should be available to a wide audience. Much of the modern astronomy pipeline is written in Python.

Smoothed particle hydrodynamics (SPH) is a method for solving the equations of hydrodynamics [@Lucy:1977; @Gingold:1977; @Price:2012]. The fluid is discretized as a set of massive particles. The density field is computed as a sum over particles weighted by a smoothing kernel. The particles are evolved by solving the equations of hydrodynamics written in Lagrangian form.

Visualization of SPH data is complicated by the fact that the data is represented on a set of particles arranged irregularly in space, as opposed to a regular grid. One visualization method is to plot the particles as a scatter plot, possibly with a color to represent some quantity defined on the particles. An alternative is to interpolate quantities defined on the particles to a two-dimensional pixel array. The latter approach provides an accurate representation of the data, in that it uses the same underlying approach to computing density as in the SPH method.

SPH simulations typically output files which contain quantities defined on the particles at a given snapshot in time. The size of these files scales with the number of particles. These snapshot files may also include metadata, such as physical units, time, and other physical and numerical parameters. In addition, output may also include auxiliary files containing globally-averaged quantities output more frequently than snapshot files.

``Plonk`` can read data from these various data sets and provide access to the underlying data as ``NumPy`` arrays. Large arrays are loaded only as required. The entire data set does not need to be read into memory at once. After the array is loaded from file it is cached in memory for later use. ``Plonk`` can read particles representing gas, dust (including multiple species), or other fluids, as well as sink particles (representing compact objects like stars or planets).

``Plonk`` has tools to perform typical tasks associated with the analysis and visualization of smoothed particle hydrodynamics data, such as: calculating extra quantities on the particles, calculating and plotting radial profiles, and visualization of any quantity defined on the particles via kernel density estimation. The interpolation component of the visualization package makes use of ``KDEpy``, a Python package for weighted kernel density estimation using fast Fourier transforms [@Odland:2018]. As such it has performance on par with other (compiled) SPH visualization programs such as Splash [@Price:2007].

Given that ``Plonk`` is a Python package it can be used via a Jupyter notebook [@Kluyver:2016]. The Jupyter notebook is a platform which allows for exploratory data analysis and visualization, combined with the ability to share the notebook, as a record of analysis, with collaborators. A user can install Plonk, using a package manager such as pip or Conda, on the supercomputer where their data resides, and analyze and visualize their data via a web browser running on their local machine. The expensive computations take place on the remote machine.

Similar existing Python packages for analysis and visualization of SPH (or N-body) data such as ``pynbody`` [@Pontzen:2013] and ``Py-SPHViewer`` [@Benitez-Llambay:2015] focus on cosmological simulations. The, more general, scientific analysis and visualization package ``yt`` [@Turk:2011] was originally designed to work with gridded simulation output. The next major version of ``yt`` (version 4.0) will support analyzing SPH data directly on the particles, rather than via interpolating onto an octree data structure, as in the current version (3.5.1). In contrast, ``Plonk`` uses the particle data directly and is focussed on astrophysical fluid dynamics simulations, such as protoplanetary and black hole accretion disks, star formation, and molecular clouds. ``Plonk`` was originally designed with users of the Phantom SPH code [@Price:2018] in mind. At the time of writing, ``Plonk`` is the only open source Python package for analysis and visualization that can read Phantom SPH data.

# Acknowledgements

DM is funded by a Research Training Program Stipend from the Australian government.

Plonk relies on the following scientific Python packages: ``NumPy`` [@Oliphant:2006; @van-der-Walt:2011], ``SciPy`` [@Virtanen:2019], ``matplotlib`` [@Hunter:2007], ``h5py`` [@Collette:2013], ``KDEpy`` [@Odland:2018], ``pandas`` [@McKinney:2010], ``Astropy`` [@Astropy:2013], ``Pint``, ``scikit-image`` [@van-der-Walt:2014].

# References
