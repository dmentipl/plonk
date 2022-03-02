# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Types of changes:

- `Added` for new features.
- `Changed` for changes in existing functionality.
- `Deprecated` for soon-to-be removed features.
- `Removed` for now removed features.
- `Fixed` for any bug fixes.
- `Security` in case of vulnerabilities.

## [Unreleased]

## [0.7.4] - 2021-10-13

### Changed

- Moved from Travis CI to GitHub actions for tests/CI.
- Moved to src layout.
- Version is now set in `src/plonk/__init__.py` which is read in `setup.cfg`.
- Update MANIFEST.in.
- Renamed master branch to main. Changes reflected in docs.
- Restructure docs, including changing from reStructuredText to Markdown using myst-parser.

## [0.7.3] - 2020-08-28

### Added

- Add add_alias method to Profile.
- Add disc_viscosity, alpha_shakura_sunyaev, epicyclic_frequency, midplane_stokes_number profiles.
- Add utils vector_array_names like dust_array_names.
- Add plotting of multiple dust or vector profiles.
- Add extra test data to complement the test suite.
- Add function to set missing header values on converted Phantom-HDF5 files.

### Changed

- Renamed toomre_Q profile to toomre_q.
- Use cache context manager to generate some Snap attributes so as to not cache some arrays in memory.
- Some profiles have been removed when making non-radial profiles.
- The analysis.total module now uses sink particles in computing summed quantities.
- Moved some quantities from particles to discs as they are not generic but appropriate for discs simulations.
- Rename interpolate arguments: number_of_pixels -> num_pixels and density_weighted -> weighted. And weighted is a named argument (i.e. not just caught with **kwargs) to image, vector, and interpolate
- Rename `load_ev` to `load_time_series`.
- Rename `load_sim` to `load_simulation`.
- Replaced Snap.set_gravitational_parameter to Snap.set_central_body. A central body is required to calculate orbital dynamics quantities, e.g. eccentricity, on the particles.
- Refactored plonk.load_snap. Much of the code that was in the phantom reader module was, in fact, more general. Now, some of that code lives in Snap.load_snap which requires defining three functions in the reader modules: snap_properties_and_units, snap_array_registry, and snap_sink_registry
- Use figshare to host sample datasets not Anaconda Cloud.

### Deprecated

- `load_ev` is deprecated in favour of `load_time_series`.
- `load_sim` is deprecated in favour of `load_simulation`.

### Fixed

- Fix units labels in Profile.plot.
- Fix units of internal energy.
- Fix bugs in rotation to face-on and edge-on.
- Fix bugs in disc position and inclination angles.

## [0.7.2] - 2020-08-23

### Added

- Add label argument to Profile.plot.
- Add dust_array_names utility function. This makes a list of array names broken into sub-species.
- Add plonk.animate function to provide a common interface to animation functions.
- Always cache arrays during plotting using the caching context manager.
- Added function name to logging messages.

### Changed

- Analysis particles functions with dust quantites are now for the dust "mixture method" (i.e. 1-fluid).
- Remove '_tot' as an array suffix for dust arrays. This was causing confusion as it makes no sense to sum stopping times.
- Rename base_array_name to base_profile_name.
- Change interpolate to return dimensionful quantity.
- Moved utils functions from each sub-package to modules in utils sub-package.
- Moved animation, animation_particles, animation_profiles from plonk namespace to plonk.visualize namespace.
- Renamed animation to animation_images.

### Fixed

- Fix bug in setting default units for dimensionless arrays.
- Fix bugs in Profile related to units.
- Fix bug in std dev shading in Profile.plot.

## [0.7.1] - 2020-08-21

### Added

- Use a TOML config file to configure options.
- Add default_units and set_units on Snap and Profile.
- Sinks are iterable.
- Plot error bars on profiles.

### Changed

- Renamed subsnaps_by_type to subsnaps_as_dict.
- Renamed units_defaults to array_units.
- Name mapping, units, aliases are no longer hard coded and are now in config.toml.
- If no units specified in image/plot/vector functions then use the default units on the Snap if available.
- Renamed some analysis.particles functions.
- Sinks analysis functions take a Sinks object as arguments.
- Profile.plot "std_dev_shading" argument changed to "std", and now is a string not a bool.

### Fixed

- Fixed bug in using Snap for snaps with no sinks.
- Fixed bug in accessing a single Sink if a np.int rather than int was passed in.
- Fixed bug in reading Phantom datasets without knowing units.

## [0.7.0] - 2020-08-17

### Added

- Added plonk.image to make image plots (with interpolation and then matplotlib imshow).
- Added plonk.vector to make vector plots (with interpolation and then matplotlib quiver).
- Added plot_smoothing_length function to plot the smoothing length on particles, or accretion radius on sink particles.
- Added pretty_array_name function to prettify array names.
- Added visualize_sim as a method on Simulation objects.
- Allow getting subsets of Sinks.
- Added ax_kwargs to plotting functions for passing to ax.set.
- Added xlim, ylim on visualize.plot.
- Added units_dict function to return a units dictionary for passing to plotting functions.
- Added bulk load and unload functions on Snap for loading/unloading multiple arrays into/out-of memory.
- Add context manager for caching arrays on Snap.
- Added public methods on Snap: family for accessing particle families, and array for accessing particle arrays. These are already accessible via __getitem__ but this makes the underlying methods available.
- Add function to add missing units to array on Snap.

### Changed

- Removed plonk.particle_plot in favour of plonk.plot.
- Changed plonk.plot to produce particle plots
- Renamed MultiVisualization and plot_snaps to VisualizeSimulation and visualize_sim.
- Changed units system from cgs to SI.
- Simplified animation functions by only allowing one axes per animation.
- Changed default units to more conventional SI units, e.g. Pascal for pressure and Joule for energy.
- Simplified tree and neighbours functions on Snap. Now there is only one tree for a Snap. If you want a tree for just, say, dust particles, then first create a SubSnap and get the tree on that.
- Changed _Sinks into Sinks, i.e. a public class.
- All plotting functions/methods use the same argument for setting units.
- Renamed Snap.available_arrays argument "all" to "verbose".
- Changed Snap.units to Snap.code_units.
- Use pretty_array_name for plots labels.
- Rename Snap.unset to Snap.reset and allow for unloading cached arrays.
- When setting Snap.cache_arrays to False, no longer unload currently cached arrays.

### Fixed

- Fixed writing Snap.to_dataframe units.

## [0.6.2] - 2020-08-11

### Changed

- Use setup.cfg for setuptools, and pyproject.toml (and setup.cfg) for config of tools.
- Version is set in setup.cfg and imported into plonk via importlib_metadata.
- Changed API documentation.
- Moved sph module from utils sub-package to analysis.

### Deprecated

- plonk.particle_plot will be removed.
- plonk.plot will change from image plots to particle plots, and plonk.image and plonk.vector will be added to replace plonk.plot.
- Default units will change from cgs to SI.

### Fixed

- Fixed bug in Profile with getting number of mixture dust species.
- Fixed bugs in animation functions (due to making physical units on by default).
- Fixed issues with colorbar size matching height of plots.

## [0.6.1] - 2020-08-09

### Added

- Snap.sinks attribute has more features.
- Cross section interpolation in a non-xy plane specified by a normal vector to the plane.
- Snap.rotate can be set by an axis vector and angle as opposed to a scipy Rotation object.
- discs module to analysis.
- filters module to analysis to set SubSnaps easily.
- 'id' array on Snap to help track particles.
- Function to plot the smoothing length as a circle.
- Profile method to generate a function from a profile to help create particle filters, for example.
- Simulation method to create a particle array over the whole simulation.

### Changed

- Snap.available_arrays does not reference sink particles; see Snap.sinks.available_arrays.
- Profile.plot units are now consistent with visualize functions.
- Dust profiles in Profile are now distinguished by whether they are mixture (dust/gas) particles or dust-only particles.

### Fixed

- Setting origin in extra quantities.
- All analysis functions have better physical units support.
- Bug in Snap.num_particles_of_type.

## [0.6.0] - 2020-08-05

### Added

- Added plot and particle_plot as methods of the Snap class. This allows for plotting with `snap.plot(quantity='quantity')` as opposed to `plonk.visualize.plot(snap=snap, quantity='quantity)`.
- Axis and colorbars have labels by default now, including units.
- The to_dataframe Snap method now indicates units in the column names, e.g. `position [au]`.
- The available_arrays method of Snap has additional arguments to see all sub-arrays on particles, e.g. `velocity_x` and `dust_fraction_001`.
- Added to examples and quick-start in documentation.
- Added method to re-open a closed Snap file.

### Changed

- Physical units are turned on by default on Snap objects. All particle and sink arrays have units (provided by Pint).
- The units attribute of Snap and Simulation now only has core units, i.e. length, time, mass, and magnetic field.
- Some extra quantities have been renamed.
- Extra quantities are available on Snap objects by default.
- The arguments radius_min and radius_max in Profile have been renamed cmin and cmax to reflect that profiles are not just radial.

### Fixed

- Fixed setting pressure from Phantom equation of states.

## [0.5.3] - 2020-07-28

### Fixed

- Fixed major bug in setting extra dust arrays on a Snap.

## [0.5.2] - 2020-07-28

### Added

- Change log.
- Cartesian profiles in y- and z-direction, in addition to x-direction which was already implemented.

### Changed

- Do not raise exception in extra_quantities and physical_units if already set.
- Scikit-image and tqdm are no longer required dependencies.
- Conda environment renamed from plonk-dev to plonk.
- Refactor Plonk namespace. Fewer modules are directly imported.

## [0.5.1] - 2020-07-11

### Added

- Analysis functions for sink particles.
- Function to animate particle plots.
- Different colours per particle type in particle plots.
- Tqdm progress bar for animation.

### Changed

- Use Matplotlib consistent argument names in particle_plot.

### Fixed

- Fix bug in standard deviation shading in Profile.

## [0.5.0] - 2020-04-20

### Added

- Neighbour finding via kd-tree.
- Compute SPH derivatives using kd-tree.
- IPython tab completion for Snap arrays and Profile profiles.
- Profile can have ndim==1 which gives a linear profile, useful for box calculations.
- Option to turn off caching of particle arrays, so that they are always read from file.
- Write derived arrays to HDF5 file, and read arrays from that file onto a Snap.
- Added logging of warning and other information messages.

### Changed

- Generalize sub-types: dust_type → sub_type: this allows for Phantom boundary particle sub-types.

### Removed

- Remove `Visualization` class in favour of just returning matplotlib's Axes and Figure objects.

## [0.4.1] - 2020-03-24

### Added

- Add scatter plots, i.e. particle plots with variable color and size markers.
- Add `extra_quantities` method to easily calculate extra quantities on the snapshot.
- Allow for setting array units, whether the array is rotatable, and whether it is a dust on derived arrays.
- Profiles are automatically generated from any 1d Snap quantity.
- Access individual dust quantities on profiles via '_001', etc.
- Read Phantom equation of state information to get pressure, sound speed, temperature.
- Add extra Snap quantities, e.g. Stokes number.
- Add extra profiles, e.g. Toomre Q.
- Allow accessing components of Snap quantities via '_x', etc.
- Calculate standard deviation on profiles.

### Changed

- Use verbose names for all snapshot quantities, e.g. 'dust_fraction' not 'dustfrac' and 'velocity_divergence' not 'divv'.

### Removed

- Remove `Evolution` object in favour of pandas DataFrame.

## [0.4.0] - 2020-03-15

### Added

- Add physical units on the `Snap` object.
- Physical units are compatible with all visualization and analysis modules.

## [0.3.1] - 2020-03-06

### Added

- Add many analysis functions.
- Animations of visualizations.

### Changed

- Make it easier to add profiles to Profile
- Make `plonk.visualize.plot` easier to use.

### Fixed

- Fix bug in `Snap.rotate` not rotating sink particles.

## [0.3.0] - 2019-12-07

### Changed

- Interpolation functions are now written in Python and JIT-compiled with Numba.

## [0.2.1] - 2019-11-27

### Added

- Add the JOSS paper describing the code.

## [0.2.0] - 2019-11-06

### Changed

- Use KDEpy for interpolation.

## [0.1.0] - 2019-06-28

- Initial release.
