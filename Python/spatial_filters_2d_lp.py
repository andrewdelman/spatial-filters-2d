def spat_filter_2d_lp(curr_var,lon_array,lat_array,lat_deg_threshold,deg_threshold,ratio_side_deg_threshold=2.5,ratio_side_spacing=3,steepness_factor=2,n_iter_coastbuffer=4):
    
    """
    This function applies spatial 2-D low-pass filter in wavenumber space, to a given field curr_var in a rectangular grid domain that is not periodic at the boundaries. If it is desired for the field to connect laterally across the boundaries of the domain (e.g., global tripole grid or LLC tiles), this stitching must be done manually by the user before calling the function.
    
    Inputs:
    
    lon_array, lat_array: 2-D Numpy arrays; contain longitude and latitude coordinates respectively. lon_array may contain abrupt ~360 degree jumps due to the branch cut (e.g., 360 to 0 deg, or +180 to -180 deg). Filters should propagate across these jumps, but not across array boundaries.
    
    curr_var: 2-D or 3-D Numpy array; first 2 dimensions correspond to horizontal dimensions defined by lon_array and lat_array
    
    lat_deg_threshold: list, tuple, or 1-D Numpy array; latitude base points to define wavelength thresholds for filtering. If a consistent filter threshold is desired across the domain, set lat_deg_threshold to [np.nanmin(lat_array),np.nanmax(lat_array)].
    
    deg_threshold: list, tuple, or 1-D Numpy array; wavelength thresholds for filtering, to be interpolated between latitude base points. A value of 1 corresponds to a half-power cutoff wavelength of 1 degree latitude, or approx. 111.1 km. If a consistent filter threshold is desired across the domain, set deg_threshold to [threshold,threshold] or threshold*np.ones((2,)).
    
    ratio_side_deg_threshold: float; ratio of side length of filtering tiles to threshold wavelength (must be somewhat greater than 1; recommended range is 2 to 4). Default value is 2.5.
    
    ratio_side_spacing: float; ratio of side length of filtering tiles to spacing between tile centers (must be somewhat greater than 1; recommended range is 2 to 5). Default value is 3.
    
    steepness factor: float; sets the steepness of the wavenumber taper near the wavelength cutoff defined in deg_threshold. Recommended values between 1 and 5; lower values have a more gradual (less precise) cutoff, while higher values are more precise but may have spectral ringing especially near sharp peaks. Default value is 2.
    
    n_iter_coastbuffer: non-negative integer; number of iterations to apply filter in order to smooth transitions across coastlines (in curr_var, NaN, Inf, and zero values are assumed to be land areas). This value can be set to zero, which will save a bit of computing time but this is only recommended in open ocean areas where the edges of the domain are not the focus of the analysis. Otherwise, a value of at least 4 is recommended for decent conservation of large-scale integrals near coastlines or the edges of the domain. Default value is 4.
    
    
    Outputs:
    
    curr_var_lp: array the same size as curr_var, consisting of the values of curr_var spatially low-pass filtered in the last 2 dimensions of the array
    
    wet_mask_spat_lp: array the same size as curr_var, consisting of the spatially low-passed "wet" mask (open ocean values near 1, land values near 0, with a transition in values corresponding to the filtering scale near coastlines)
    """
    
    #---------------------------------------------------------------------------
    
    import numpy as np
    from scipy import special
    
    
    # define NaN mask array
    
    orig_var = curr_var
    
    if len(curr_var.shape) >= 3:
        size_tup = curr_var.shape
    elif len(curr_var.shape) == 2:
        size_tup = (1,) + curr_var.shape
        curr_var = np.expand_dims(curr_var,axis=0)
    
    
    wet_mask_array = np.ones(curr_var.shape)
    max_abs_curr_var = np.nanmax(np.abs(curr_var[np.logical_and(~np.isnan(curr_var),~np.isinf(curr_var))]))
    wet_mask_array[np.logical_or(np.logical_or(np.isnan(curr_var) == 1,np.isinf(curr_var) == 1),np.abs(curr_var) < (1.e-10)*max_abs_curr_var)] = 0
    
    curr_var[wet_mask_array < 1.e-5] = np.nan
        
    
    
    # define center locations for tapered bins
    
    diff_x_lat_array = np.abs(np.diff(lat_array,axis=-1))
    diff_y_lat_array = np.abs(np.diff(lat_array,axis=-2))
    sorted_diff_lat = np.sort(np.hstack((diff_x_lat_array[~np.isnan(diff_x_lat_array)],diff_y_lat_array[~np.isnan(diff_y_lat_array)])))
    typical_diff_lat = sorted_diff_lat[int(np.round(0.75*len(sorted_diff_lat)))]    
    
    lat_deg_threshold = np.array(lat_deg_threshold)
    deg_threshold = np.array(deg_threshold)
    
    lat_center_decomp = []
    deg_threshold_center = []
    ratio_spacing_deg_threshold = ratio_side_deg_threshold/ratio_side_spacing
    curr_lat_center_spacing = 360/(np.round(360/(ratio_spacing_deg_threshold*deg_threshold[0])))
    curr_lat_center = lat_deg_threshold[0] + curr_lat_center_spacing
    while curr_lat_center <= np.max(lat_deg_threshold):
        closest_south_ind = (lat_deg_threshold < curr_lat_center).nonzero()[0][-1]
        north_weight = (curr_lat_center - lat_deg_threshold[closest_south_ind])/np.diff(lat_deg_threshold[closest_south_ind + [0,1]])[0]
        curr_deg_threshold = np.sum(np.array([(1 - north_weight),north_weight])*deg_threshold[closest_south_ind + [0,1]])
        lat_center_decomp.append(curr_lat_center)
        deg_threshold_center.append(curr_deg_threshold)
        
        curr_lat_center_spacing = 360/(np.round(360/(ratio_spacing_deg_threshold*curr_deg_threshold)))
        curr_lat_center = curr_lat_center + curr_lat_center_spacing
    
    in_lat_range_ind = np.logical_and(np.array(lat_center_decomp) > np.nanmin(lat_array),\
                                      np.array(lat_center_decomp) < np.nanmax(lat_array))
    lat_center_decomp = list(np.array(lat_center_decomp)[in_lat_range_ind])
    deg_threshold_center = list(np.array(deg_threshold_center)[in_lat_range_ind])
    
    lon_center_decomp_list = []
    for curr_lat_decomp,curr_deg_threshold in zip(lat_center_decomp,deg_threshold_center):
        curr_lon_center_spacing = 360/(round(360/(ratio_spacing_deg_threshold*curr_deg_threshold/(np.cos((np.pi/180)*curr_lat_decomp)))))
        near_lat_ind = (abs(lat_array - curr_lat_decomp) < (2*typical_diff_lat));
        lon_center_decomp_list.append(np.arange(np.nanmin(lon_array[near_lat_ind]),np.nanmax(lon_array[near_lat_ind]),curr_lon_center_spacing))
    
    
    
    
    # tapered spatial weighting (in x and y)
    
    waven_x_low_bound = 1/(4*360)
    waven_y_low_bound = waven_x_low_bound
    
    
    cum_spat_weighting_array = np.zeros(size_tup[-2:])
    curr_var_lp = np.zeros(curr_var.shape)
    wet_mask_spat_lp = np.zeros(curr_var.shape)
    for count,(curr_center_lat,curr_deg_threshold,lon_center_decomp) in \
      enumerate(zip(lat_center_decomp,deg_threshold_center,lon_center_decomp_list)):
        curr_lat_center_spacing = 360/(round(360/(ratio_spacing_deg_threshold*curr_deg_threshold)))
        
        curr_min_lat = curr_center_lat - ((ratio_side_spacing/2)*curr_lat_center_spacing)
        curr_max_lat = curr_center_lat + ((ratio_side_spacing/2)*curr_lat_center_spacing)
        
        for curr_center_lon in lon_center_decomp:
            curr_min_lon = curr_center_lon - ((ratio_side_spacing/2)*(curr_lat_center_spacing/(np.cos((np.pi/180)*curr_center_lat))))
            curr_max_lon = curr_center_lon + ((ratio_side_spacing/2)*(curr_lat_center_spacing/(np.cos((np.pi/180)*curr_center_lat))))
            
            in_box_ind = (np.logical_and(np.logical_and(((lon_array - curr_min_lon + 180) % 360) - 180 >= 0,((lon_array - curr_max_lon + 180) % 360) - 180 <= 0),np.logical_and(lat_array >= curr_min_lat,lat_array <= curr_max_lat))).nonzero()
            in_box_i_ind = in_box_ind[-1]
            in_box_j_ind = in_box_ind[-2]
            
            if len(in_box_i_ind) < 5:
                continue
            
            in_box_i_ind = np.arange(np.min(in_box_i_ind),np.max(in_box_i_ind) + 1)
            in_box_j_ind = np.reshape(np.arange(np.min(in_box_j_ind),np.max(in_box_j_ind) + 1),(-1,1))

            
            if np.logical_or(len(in_box_i_ind) < 5,len(in_box_j_ind) < 5):
                continue
            
            
            curr_wet_mask_in_box = wet_mask_array[:,in_box_j_ind,in_box_i_ind]
            
            if np.sum(np.abs(curr_wet_mask_in_box[0,:,:]) < 1.e-5) > 0.8*(curr_wet_mask_in_box[0,:,:].size):
                continue
    
            curr_lon_in_box = lon_array[in_box_j_ind,in_box_i_ind]
            curr_lat_in_box = lat_array[in_box_j_ind,in_box_i_ind]
            unfilt_var_in_box = curr_var[:,in_box_j_ind,in_box_i_ind]
            
            
            # define weighting array within current tile, tapering to zero towards edges of tile
            if np.abs(curr_center_lat) >= 80:
                dist_from_center = np.abs(((90 - np.abs(curr_lat_in_box))*np.exp(1j*(np.pi/180)*curr_lon_in_box)) - ((90 - np.abs(curr_center_lat))*np.exp(1j*(np.pi/180)*curr_center_lon)))
                max_dist_from_center = (ratio_side_spacing/2)*curr_lat_center_spacing
                dist_norm_from_center = dist_from_center/((2**(-1/2))*max_dist_from_center)
                curr_spat_weighting_array = (np.cos((np.pi/2)*dist_norm_from_center/(2.5/ratio_side_spacing)))**2;
                curr_spat_weighting_array[dist_norm_from_center/(2.5/ratio_side_spacing) > 1] = 0
            else:
                below_center_dist_x_norm = (((curr_lon_in_box - curr_center_lon + 180) % 360) - 180)/(curr_min_lon - curr_center_lon)
                above_center_dist_x_norm = (((curr_lon_in_box - curr_center_lon + 180) % 360) - 180)/(curr_max_lon - curr_center_lon)
                below_center_dist_y_norm = (curr_lat_in_box - curr_center_lat)/(curr_min_lat - curr_center_lat)
                above_center_dist_y_norm = (curr_lat_in_box - curr_center_lat)/(curr_max_lat - curr_center_lat)
                dist_x_norm_from_center = np.fmax(below_center_dist_x_norm,above_center_dist_x_norm)
                dist_y_norm_from_center = np.fmax(below_center_dist_y_norm,above_center_dist_y_norm)
                dist_from_center_renorm = np.abs(dist_x_norm_from_center + (1j*dist_y_norm_from_center))/(2.5/ratio_side_spacing)
                curr_spat_weighting_array = (np.cos((np.pi/2)*dist_from_center_renorm))**2
                curr_spat_weighting_array[dist_from_center_renorm > 1] = 0
    
            
            waven_y_high_bound = 1/curr_deg_threshold
            waven_x_high_bound = waven_y_high_bound
            
            diff_x_dist = np.abs((np.cos((np.pi/180)*(curr_lat_in_box[:,1:] - (np.diff(curr_lat_in_box,axis=-1)/2)))*(((np.diff(curr_lon_in_box,axis=-1) + 180) % 360) - 180)) + (1j*np.diff(curr_lat_in_box,axis=-1)))
            diff_y_dist = np.abs((np.cos((np.pi/180)*(curr_lat_in_box[1:,:] - (np.diff(curr_lat_in_box,axis=-2)/2)))*(((np.diff(curr_lon_in_box,axis=-2) + 180) % 360) - 180)) + (1j*np.diff(curr_lat_in_box,axis=-2)))
            delta_x = np.nanmean(diff_x_dist)
            delta_y = np.nanmean(diff_y_dist)
                    
            
            if len(unfilt_var_in_box.shape) > 2:
                size_tup_nested = ()
                for num in range(len(unfilt_var_in_box.shape) - 2):
                    size_tup_nested += (1,)
                size_tup_nested = tuple(np.array(size_tup_nested + (2,2))*np.array(unfilt_var_in_box.shape))
            else:
                size_tup_nested = tuple(np.array([2,2])*np.array(unfilt_var_in_box.shape))
            
            if n_iter_coastbuffer < 1:
                curr_var_in_box = unfilt_var_in_box
            else:
                # define buffer over land areas near coast, using 2-D low-pass filter
                
                steepness_factor_buffer = np.fmin(steepness_factor,2)
                
                curr_half_power_adj = np.exp(special.erfinv((2**(1/2)) - 1)/steepness_factor_buffer)   # adjustment factor to set bounds at half-power (rather than half-amplitude)
                
                
                # find planar fit for each level
                curr_var_planar_fit = planar_fit(unfilt_var_in_box)
                
                
                # iterate to buffer values in masked areas near boundaries
                
                curr_var_in_box = unfilt_var_in_box - curr_var_planar_fit;
                curr_var_in_box[np.logical_or(np.isnan(curr_var_in_box),np.abs(curr_wet_mask_in_box) < 1e-5)] = 0
                
                bandpass_filter = error_fcn_bp_filter_2d(size_tup_nested,waven_x_low_bound,waven_x_high_bound,delta_x,delta_y,steepness_factor_buffer,power_cutoff_opt=1)
                for iter in range(n_iter_coastbuffer):
                    curr_var_in_box = bp_filter_2d_apply(curr_var_in_box,bandpass_filter)
                    
                    curr_var_in_box[np.abs(curr_wet_mask_in_box - 1) < 1.e-5] = unfilt_var_in_box[np.abs(curr_wet_mask_in_box - 1) < 1.e-5] - curr_var_planar_fit[np.abs(curr_wet_mask_in_box - 1) < 1.e-5]
                
                
                curr_var_in_box += curr_var_planar_fit
            
            
            
            # apply 2-D low pass filter after buffering
            
            bandpass_filter = error_fcn_bp_filter_2d(size_tup_nested,waven_x_low_bound,waven_x_high_bound,delta_x,delta_y,steepness_factor,power_cutoff_opt=1)
            
            curr_array_planar_fit = planar_fit(curr_var_in_box)
            
            curr_array = curr_var_in_box - curr_array_planar_fit
            curr_array[np.logical_or(np.isnan(curr_array),np.abs(curr_wet_mask_in_box) < 1.e-5)] = 0
            curr_var_in_box_lp = bp_filter_2d_apply(curr_array,bandpass_filter) + curr_array_planar_fit
            
            
            curr_array_planar_fit = planar_fit(curr_wet_mask_in_box)
            
            curr_array = curr_wet_mask_in_box - curr_array_planar_fit
            curr_array[np.logical_or(np.isnan(curr_array),np.abs(curr_wet_mask_in_box) < 1.e-5)] = 0
            curr_wet_mask_spat_lp = bp_filter_2d_apply(curr_array,bandpass_filter) + curr_array_planar_fit
            
            
            # add filtered current tile to full-size array
            cum_spat_weighting_array[in_box_j_ind,in_box_i_ind] += curr_spat_weighting_array
            curr_var_lp[:,in_box_j_ind,in_box_i_ind] += (curr_spat_weighting_array*curr_var_in_box_lp)
            wet_mask_spat_lp[:,in_box_j_ind,in_box_i_ind] += (curr_spat_weighting_array*curr_wet_mask_spat_lp)
        
        
        print('Completed latitude row number ' + str(count) + ' (out of ' + str(len(lat_center_decomp)) + ')')
        
    
    # normalize by cumulative weighting array
    curr_var_lp = curr_var_lp/cum_spat_weighting_array
    wet_mask_spat_lp = wet_mask_spat_lp/cum_spat_weighting_array
    
    
    # reshape arrays to size of original array
    curr_var_lp = np.reshape(curr_var_lp,orig_var.shape)
    wet_mask_spat_lp = np.reshape(wet_mask_spat_lp,orig_var.shape)
    
    # reset masked areas to original values
    mask = (np.abs(np.reshape(wet_mask_array,orig_var.shape)) < 1.e-5)
    curr_var_lp[mask] = orig_var[mask]
    
        
    return curr_var_lp,wet_mask_spat_lp


#===============================================================================


def planar_fit(input_array):
    "Compute linear regression planar fit each 2-D level (last 2 dimensions) of a 3-D Numpy array"
    
    #---------------------------------------------------------------------------
    
    import numpy as np
    
    planar_fit_array = np.empty(input_array.shape)
    planar_fit_array.fill(np.nan)
    for z in range(planar_fit_array.shape[-3]):
        curr_array = input_array[z,:,:]
        good_ind = np.logical_and(~np.isnan(curr_array),~np.isinf(curr_array))
        if np.sum(good_ind) < 0.2*(curr_array.size):
            G = np.ones((curr_array.size,))
            G_good = G[good_ind]
            m = np.dot(np.linalg.inv(np.dot(G_good,G_good)),np.dot(G_good,curr_array[good_ind]))
        else:
            G = np.hstack((np.ones((curr_array.size,1)),\
                           np.reshape(np.tile(np.reshape(np.arange(curr_array.shape[-2]),(-1,1)),(1,curr_array.shape[-1])),(-1,1)),\
                           np.tile(np.reshape(np.arange(curr_array.shape[-1]),(-1,1)),(curr_array.shape[-2],1))))
            G_good = np.hstack((np.ones((int(np.sum(good_ind)),1)),\
                                np.reshape(good_ind.nonzero()[-2],(-1,1)),\
                                np.reshape(good_ind.nonzero()[-1],(-1,1))))
            m = np.dot(np.linalg.inv(np.matmul(G_good.transpose(),G_good)),np.dot(G_good.transpose(),curr_array[good_ind]))
        
        planar_fit_array[z,:,:] = np.reshape(np.dot(G,m),curr_array.shape)
    
    return planar_fit_array


#===============================================================================


def error_fcn_bp_filter_2d(size_tup_nested,low_bound,high_bound,delta_x,delta_y,steepness_factor=2,power_cutoff_opt=1):
    """
    Compute 2-D bandpass filter in spectral (wavenumber) space. This function does not actually apply a filter, but outputs the transfer /coefficients in 2-D wavenumber space for an array of size specified by tuple array_shape.
    """
    
    #---------------------------------------------------------------------------
    
    import numpy as np
    from scipy import special
    
    if power_cutoff_opt == 1:
        half_power_adj = np.exp(special.erfinv((2**(1/2)) - 1)/steepness_factor)   # adjustment factor to set bounds at half-power (rather than half-amplitude)
    else:
        half_power_adj = 1
    
    k_array = np.tile((1/(delta_x*size_tup_nested[-1]))*(((np.arange(size_tup_nested[-1]) + (size_tup_nested[-1]/2)) % size_tup_nested[-1]) - (size_tup_nested[-1]/2)),(size_tup_nested[-2],1))
    l_array = np.tile((1/(delta_y*size_tup_nested[-2]))*(((np.reshape(np.arange(size_tup_nested[-2]),(-1,1)) + (size_tup_nested[-2]/2)) % size_tup_nested[-2]) - (size_tup_nested[-2]/2)),(1,size_tup_nested[-1]))
    mag_waven_array = ((k_array**2) + (l_array**2))**(1/2)
    low_bound = low_bound/half_power_adj
    high_bound = high_bound*half_power_adj
    bandpass_filter = 0.5*(special.erf(steepness_factor*(np.log(np.abs(mag_waven_array)) - np.log(low_bound))) - special.erf(steepness_factor*(np.log(np.abs(mag_waven_array)) - np.log(high_bound))))
    
    return bandpass_filter
    
    
#===============================================================================


def bp_filter_2d_apply(array_to_filter,bandpass_filter):
    """
    Apply 2-D bandpass filter in spectral (wavenumber) space.
    
    Inputs:
    array_to_filter: 2-D or 3-D Numpy array
    bandpass_filter: 2-D Numpy array of filter transfer coefficients in spectral (wavenumber) space, as computed by error_fcn_bp_filter_2d
    
    Output:
    array_filtered: 2-D or 3-D Numpy array
    
    """
    
    #---------------------------------------------------------------------------
    
    import numpy as np
    from scipy import special

    if len(array_to_filter.shape) > 2:
        size_tup_nested = tuple(np.array([1,2,2])*np.array(array_to_filter.shape))
    else:
        size_tup_nested = (1,) + tuple(np.array([2,2])*np.array(array_to_filter.shape))
    
    i_ind_nested = int(np.floor((array_to_filter.shape[-1])/2)) + np.arange(array_to_filter.shape[-1])
    j_ind_nested = np.reshape(int(np.floor((array_to_filter.shape[-2])/2)) + np.arange(array_to_filter.shape[-2]),(-1,1))
    curr_array_nested = np.zeros(size_tup_nested)
    curr_array_nested[:,j_ind_nested,i_ind_nested] = array_to_filter
    fft_2d_curr_array = np.fft.fft(np.fft.fft(curr_array_nested,axis=-2),axis=-1)
    curr_array_filtered = np.real(np.fft.ifft(np.fft.ifft(np.expand_dims(bandpass_filter,axis=0)*fft_2d_curr_array,axis=-2),axis=-1))
    array_filtered = np.reshape(curr_array_filtered[:,j_ind_nested,i_ind_nested],array_to_filter.shape)
    
    return array_filtered