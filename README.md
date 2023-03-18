# spatial-filters-2d
Two-dimensional spatial filters, developed for application to ocean and atmosphere model output and gridded observations

## Introduction
The functions in this repository apply various two-dimensional spatial filters to gridded data. These data can be output from a model simulation, or observations in a gridded layout (e.g., from satellite-based products, reanalyses, Argo). There is no requirement that the grid points be evenly spaced or aligned with longitude or latitude lines; however, the filters will work better if any curvature in the grid occurs on scales larger than the filter threshold scale. For input fields being filtered, non-finite (NaN, -Inf/Inf) or zero values are assumed to be land areas, and buffering and/or adjustments are implemented to reduce (but not entirely eliminate) variance loss at the land boundaries. Currently all of these functions are in Matlab, and should also work in [GNU Octave](https://octave.org) freeware--please let me know if there are issues implementing functions in Octave). Python versions will hopefully be developed soon.

## Low-pass filters

### General description
The most widely applicable filtering functions in this repository implement 2-D low-pass filters. The residual of applying these functions can also be used to isolate smaller scales, and these filters can also be combined to effectively "bandpass" for a specific range of spatial scales. The core filter in these functions is applied as a window in 2-D wavenumber space, after a planar fit to the current tile has been removed and the tile is nested in a larger 2x2 tile of zeros to buffer the edges. The transfer function for the 2-D spectral coefficients is

$$
F_f(k,l) = \left[0.5 - 0.5\rm{erf}\left( s \, \rm{ln}\frac{K}{K_0}\right)\right] F(k,l)
$$

where $F(k,l)$ is the input wavenumber coefficient for wavenumbers $k$ and $l$ with wavenumber magnitude $K = (k^2 + l^2)^{1/2}$, $K_0$ is the wavenumber representing the half-amplitude cutoff (with a slight adjustment from the half-power cutoff wavelength specified in `deg_threshold`), $s$ is `steepness_factor` representing how abrupt the cutoff is, and $F_f(k,l)$ is the filtered wavenumber coefficient. This filter is a two-dimensional extension of the one-dimensional filter used in Delman and Lee ([2020](https://doi.org/10.5194/os-16-979-2020); [2021](https://doi.org/10.5194/os-17-1031-2021)).

The low-pass filters here should yield broadly similar results to those implemented in the Python package [gcm-filters](https://github.com/ocean-eddy-cpt/gcm-filters). The spatial-filters-2d low-pass filters are most similar to the `TAPER` filters in gcm-filters. Note that the [definition](https://gcm-filters.readthedocs.io/en/latest/theory.html) of wavenumber/wavelength in gcm-filters is different from the half-power thresholds used as inputs in the functions here. Another key difference from gcm-filters is that spatial-filters-2d functions take longitude and latitude 2-D arrays as inputs and compute the relevant grid spacings automatically. However, if conservation of area- or transect-integrated quantities is required, then it is necessary for the user to weight the input field with the relevant grid cell area or face length.

The low-pass filters in spatial-filters-2d take as input half-power wavelength thresholds defined in the `deg_threshold` vector (in units of degrees latitude) with a corresponding `lat_deg_threshold` vector indicating the latitude of each entry in `deg_threshold`. The threshold cutoff wavelengths are linearly interpolated between the latitudes defined in `lat_deg_threshold`. For example, if `lat_deg_threshold` = [-90 0 90] and `deg_threshold` = [2 20 2], then the half-power cutoff wavelength is 2 deg lat $\approx$ 222.2 km at the poles, 20 deg lat $\approx$ 2222 km at the equator, and 11 deg lat $\approx$ 1222 km at 45 deg S and N.

Other aspects of the low-pass filtering process are described here:

### Overlapping tiles
In order to apply the filter at different latitudes with different threshold scales, the domain is first divided into a set of overlapping squares or tiles.  There are two parameters `ratio_side_deg_threshold` and `ratio_side_spacing` that set the size of the tile squares and the spacing between them, respectively. `ratio_side_deg_threshold` has a recommended range between 2 and 4, with 2.5 or 3 usually an optimal choice. If `ratio_side_deg_threshold` = 2.5 and the local threshold wavelength is 10 degrees latitude, then the side of each square is 25 degrees latitude wide (and the equivalent distance in longitude).  It is important that this value be somewhat larger than 1 (so that the threshold scale can be resolved) but not too large (otherwise you can get strange artifacts from the filtering especially near coastlines). `ratio_side_spacing` has a recommended range between 2 and 5, with 3 often an optimal choice. If `ratio_side_spacing` = 3 with squares of 25 degrees latitude on a side, then the centers of adjacent squares are separated by $\approx$ 8.3 degrees latitude (and equivalent distance in longitude). If this parameter is set to 1 then the squares do not overlap, but this is not recommended, because the result will have filtering artifacts near the boundaries of each square. The overlapping squares are blended together using a tapered weighting function of the form $[1 + \rm{cos}(r/R)]/2$, where $r$ is the radius from the center of the square scaled by a scale radius $R$. This function weights values close to the center of the square highest, and tapers to zero in both the weighting value and its first derivative near the edge of the square. Thus first-order spatial derivatives should not have abrupt transitions due to the tiling scheme (higher-order derivatives could though). At the end of each function the filtered field is normalized by the sum of the weights to generate the blended low-passed field.

### Coastal buffering and horizontal velocity field adjustments
To improve the filtered fields near coastlines (as determined by the presence of NaN, -Inf/Inf, or zero values), a coastal "buffer" is prepared. Before the final application of the filter, the field minus planar fit is filtered a number of times as specified by `n_iter_coastbuffer`, and after each application, the unfiltered values (minus planar fit) are restored to water grid cells. The effect is that after this process, the water grid cells are unchanged, but the land cells near the coastlines are buffered, so that the transitions in the low-pass filtered field are not abrupt near the coastlines, and attenuation of the filtered field near land areas is much reduced (though not entirely eliminated). Once the planar fit is restored, a new planar fit is computed and removed again prior to applying the filter one final time, with the filtered values over water grid cells retained. The fully low-passed field then has the planar fit restored, the original tile is recovered from the larger nested tile, and land areas are re-masked.

For the functions specifically intended to handle horizontal velocities (which contain `vel` in the name), an additional adjustment is applied. Namely, the low-pass filtered "wet" mask (with values of 1 for water or air grid cells and 0 for land) is used to compute a local angle for the low-passed coastline, and the coast-perpendicular component of the velocity field is set to zero before the final application of the filter. Note that the objective of this adjustment is *not* to eliminate low-passed flow across the real coastline, which would be impossible while retaining the scale properties of the filtered field. This correction merely reduces the low-passed cross-shore flow while retaining the low-passed alongshore component, a useful adjustment for improving large-scale transport representation in the presence of strong near-coastal flows (e.g., western boundary currents).

### Which functions to use when
spat_filter_2d_lp.m: For filtering of most scalar quantities (including vertical velocity) in a non-periodic, locally Cartesian grid. This is best used for regional subsets of model output or data, though it could be applied to more complex global grids if the user manually stitches the connections of the grid together before calling the function. It can also be used to filter horizontal velocities if the user does not want the near-coastal velocity field adjustment described above.

spat_filter_vel_2d_lp.m: For filtering of horizontal velocity components $u$, $v$ in a non-periodic, locally Cartesian grid. The near-coastal velocity field adjustment is implemented to reduce flow across the low-passed (smoothed) coastline.
