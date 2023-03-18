function [curr_var_lp,wet_mask_spat_lp] = spat_filter_2d_lp(lon_array,lat_array,curr_var,lat_deg_threshold,deg_threshold,ratio_side_deg_threshold,ratio_side_spacing,steepness_factor,n_iter_coastbuffer)

% This function applies spatial 2-D low-pass filter in wavenumber space, to a given field curr_var in a rectangular grid domain that is not periodic at the boundaries. If it is desired for the field to connect laterally across the boundaries of the domain (e.g., global tripole grid or LLC tiles), this stitching must be done manually by the user before calling the function. A related function for the global tripole grid (for POP/CESM) is under development (available upon request) and will be released soon.
% 
% Inputs:
% 
% lon_array, lat_array: 2-D arrays with lon/lat coordinates
% 
% curr_var: 2-D or 3-D array; first 2 dimensions correspond to horizontal dimensions defined by lon_array and lat_array
% 
% lat_deg_threshold: vector, latitude base points to define wavelength thresholds for filtering. If a consistent filter threshold is desired across the domain, set lat_deg_threshold to [min(min(lat_array,'omitnan'),'omitnan') max(max(lat_array,'omitnan'),'omitnan')].
% 
% deg_threshold: vector, wavelength thresholds for filtering, to be interpolated between latitude base points. A value of 1 corresponds to a half-power cutoff wavelength of 1 degree latitude, or approx. 111.1 km. If a consistent filter threshold is desired across the domain, set deg_threshold to threshold*ones([1 2]).
% 
% ratio_side_deg_threshold: scalar, ratio of side length of filtering tiles to threshold wavelength (must be somewhat greater than 1; recommended range is 2 to 4).
% 
% ratio_side_spacing: scalar, ratio of side length of filtering tiles to spacing between tile centers (must be somewhat greater than 1; recommended range is 2 to 5).
% 
% steepness factor: scalar, sets the steepness of the wavenumber taper near the wavelength cutoff defined in deg_threshold. Recommended values between 1 and 5; lower values have a more gradual (less precise) cutoff, while higher values are more precise but may have spectral ringing especially near sharp peaks.
% 
% n_iter_coastbuffer: non-negative integer scalar, number of iterations to apply filter in order to smooth transitions across coastlines (in curr_var, NaN, Inf, and zero values are assumed to be land areas). This value can be set to zero, which will save a bit of computing time but this is only recommended in open ocean areas where the edges of the domain are not the focus of the analysis. Otherwise, a value of at least 4 is recommended for decent conservation of large-scale integrals near coastlines or the edges of the domain.
% 
% 
% Outputs:
% 
% curr_var_lp: array the same size as curr_var, consisting of the values of curr_var spatially low-pass filtered in the first 2 dimensions of the array
% 
% wet_mask_spat_lp: array the same size as curr_var, consisting of the spatially low-passed "wet" mask (open ocean values near 1, land values near 0, with a transition in values corresponding to the filtering scale near coastlines)
% 
% ===========================================================================
    
    
    
    % define NaN mask array
    
    orig_var = curr_var;
    
    if length(size(curr_var)) >= 3
        size_array = size(curr_var);
    elseif length(size(curr_var)) == 2
        size_array = [size(curr_var) 1];
    end
    wet_mask_array = ones(size_array(1:3));
    max_abs_curr_var = max(abs(curr_var((~isnan(curr_var)) & (~isinf(curr_var)))));
    wet_mask_array(((isnan(curr_var) == 1) | (isinf(curr_var) == 1)) | (abs(curr_var) < ((1e-10)*max_abs_curr_var))) = 0;
    
    curr_var(wet_mask_array < 1e-5) = NaN;
    
    
    
    % define center locations for tapered bins
    
    diff_x_lat_array = abs(diff(lat_array,1,1));
    diff_y_lat_array = abs(diff(lat_array,1,2));
    sorted_diff_lat = sort([diff_x_lat_array(isnan(diff_x_lat_array) == 0); diff_y_lat_array(isnan(diff_y_lat_array) == 0)],1,'ascend');
    typical_diff_lat = sorted_diff_lat(round(0.75*numel(sorted_diff_lat)));
    
    
    lat_center_decomp = [];
    deg_threshold_center = [];
    ratio_spacing_deg_threshold = ratio_side_deg_threshold/ratio_side_spacing;
    curr_lat_center_spacing = 360/(round(360/(ratio_spacing_deg_threshold*deg_threshold(1))));
    curr_lat_center = lat_deg_threshold(1) + curr_lat_center_spacing;
    while curr_lat_center <= max(lat_deg_threshold)
        closest_south_ind = find(lat_deg_threshold < curr_lat_center,1,'last');
        north_weight = (curr_lat_center - lat_deg_threshold(closest_south_ind))/diff(lat_deg_threshold(closest_south_ind + [0 1]));
        curr_deg_threshold = sum([(1 - north_weight) north_weight].*deg_threshold(closest_south_ind + [0 1]));
        lat_center_decomp = [lat_center_decomp curr_lat_center];
        deg_threshold_center = [deg_threshold_center curr_deg_threshold];
        
        curr_lat_center_spacing = 360/(round(360/(ratio_spacing_deg_threshold*curr_deg_threshold)));
        curr_lat_center = curr_lat_center + curr_lat_center_spacing;
    end
    
    in_lat_range_ind = find((lat_center_decomp > min(min(lat_array,[],'omitnan'),[],'omitnan')) & (lat_center_decomp < max(max(lat_array,[],'omitnan'),[],'omitnan')));
    lat_center_decomp = lat_center_decomp(in_lat_range_ind);
    deg_threshold_center = deg_threshold_center(in_lat_range_ind);
    
    lon_center_decomp_cellarray = cell(1,length(deg_threshold_center));
    for lat_center_ind = 1:length(lat_center_decomp)
        curr_lon_center_spacing = 360/(round(360/(ratio_spacing_deg_threshold*deg_threshold_center(lat_center_ind)/(cosd(lat_center_decomp(lat_center_ind))))));
        near_lat_ind = find(abs(lat_array - lat_center_decomp(lat_center_ind)) < (2*typical_diff_lat));
        lon_center_decomp_cellarray{lat_center_ind} = (min(lon_array(near_lat_ind)):curr_lon_center_spacing:max(lon_array(near_lat_ind)));
    end
    
    
    
    % tapered spatial weighting (in x and y)
    
    waven_x_low_bound = 1/(4*360);
    waven_y_low_bound = waven_x_low_bound;
    
    
    cum_spat_weighting_array = zeros(size_array(1:2));
    curr_var_lp = zeros(size(curr_var));
    wet_mask_spat_lp = zeros(size_array(1:3));
    for lat_center_ind = 1:length(lat_center_decomp)
        curr_center_lat = lat_center_decomp(lat_center_ind);
        curr_deg_threshold = deg_threshold_center(lat_center_ind);
        curr_lat_center_spacing = 360/(round(360/(ratio_spacing_deg_threshold*curr_deg_threshold)));
        
        curr_min_lat = curr_center_lat - ((ratio_side_spacing/2)*curr_lat_center_spacing);
        curr_max_lat = curr_center_lat + ((ratio_side_spacing/2)*curr_lat_center_spacing);
    
        lon_center_decomp = lon_center_decomp_cellarray{lat_center_ind};
        
        for lon_center_ind = 1:length(lon_center_decomp)
            curr_center_lon = lon_center_decomp(lon_center_ind);
            
            curr_min_lon = curr_center_lon - ((ratio_side_spacing/2)*(curr_lat_center_spacing/(cosd(curr_center_lat))));
            curr_max_lon = curr_center_lon + ((ratio_side_spacing/2)*(curr_lat_center_spacing/(cosd(curr_center_lat))));
            
            in_box_ind = find((mod(lon_array - curr_min_lon + 180,360) - 180 >= 0) & (mod(lon_array - curr_max_lon + 180,360) - 180 <= 0) & (lat_array >= curr_min_lat) & (lat_array <= curr_max_lat));
            in_box_i_ind = mod(in_box_ind - 1,size_array(1)) + 1;
            in_box_j_ind = ceil(in_box_ind/size_array(1));
            
            in_box_i_ind = (min(in_box_i_ind):1:max(in_box_i_ind))';
            in_box_j_ind = (min(in_box_j_ind):1:max(in_box_j_ind))';
            
            if ((isempty(in_box_ind) == 1) || (length(in_box_i_ind) < 5) || (length(in_box_j_ind) < 5))
                continue
            end
            
            
            curr_wet_mask_in_box = wet_mask_array(in_box_i_ind,in_box_j_ind,:);
            
            if max(sum(sum(abs(curr_wet_mask_in_box(:,:,1)) < 1e-5))) > 0.8*numel(curr_wet_mask_in_box(:,:,1))
                continue
            end
    
            curr_lon_in_box = lon_array(in_box_i_ind,in_box_j_ind);
            curr_lat_in_box = lat_array(in_box_i_ind,in_box_j_ind);
            unfilt_var_in_box = curr_var(in_box_i_ind,in_box_j_ind,:);
            
            
            % define weighting array within current tile, tapering to zero towards edges of tile
            if abs(curr_center_lat) >= 80
                dist_from_center = abs(((90 - abs(curr_lat_in_box)).*exp(1i*(pi/180)*curr_lon_in_box)) - ((90 - abs(curr_center_lat)).*exp(1i*(pi/180)*curr_center_lon)));
                max_dist_from_center = (ratio_side_spacing/2)*curr_lat_center_spacing;
                dist_norm_from_center = dist_from_center/((2^(-1/2))*max_dist_from_center);
                curr_spat_weighting_array = (cos((pi/2)*dist_norm_from_center/(2.5/ratio_side_spacing))).^2;
                curr_spat_weighting_array(dist_norm_from_center/(2.5/ratio_side_spacing) > 1) = 0;
            else
                below_center_dist_x_norm = (mod(curr_lon_in_box - curr_center_lon + 180,360) - 180)/(curr_min_lon - curr_center_lon);
                above_center_dist_x_norm = (mod(curr_lon_in_box - curr_center_lon + 180,360) - 180)/(curr_max_lon - curr_center_lon);
                below_center_dist_y_norm = (curr_lat_in_box - curr_center_lat)/(curr_min_lat - curr_center_lat);
                above_center_dist_y_norm = (curr_lat_in_box - curr_center_lat)/(curr_max_lat - curr_center_lat);
                n_in_box = length(in_box_i_ind)*length(in_box_j_ind);
                dist_x_norm_from_center = reshape(max([reshape(below_center_dist_x_norm,[n_in_box 1]) reshape(above_center_dist_x_norm,[n_in_box 1])],[],2),[length(in_box_i_ind) length(in_box_j_ind)]);
                dist_y_norm_from_center = reshape(max([reshape(below_center_dist_y_norm,[n_in_box 1]) reshape(above_center_dist_y_norm,[n_in_box 1])],[],2),[length(in_box_i_ind) length(in_box_j_ind)]);
                curr_spat_weighting_array = (cos((pi/2)*(abs(dist_x_norm_from_center + (1i*dist_y_norm_from_center)))/(2.5/ratio_side_spacing))).^2;
                curr_spat_weighting_array(abs(dist_x_norm_from_center + (1i*dist_y_norm_from_center))/(2.5/ratio_side_spacing) > 1) = 0;
            end
    
            
            waven_y_high_bound = 1/curr_deg_threshold;
            waven_x_high_bound = waven_y_high_bound;
            
            diff_x_dist = abs((cosd(curr_lat_in_box(2:length(in_box_i_ind),:) - (diff(curr_lat_in_box,1,1)/2)).*(mod(diff(curr_lon_in_box,1,1) + 180,360) - 180)) + (1i*diff(curr_lat_in_box,1,1)));
            diff_y_dist = abs((cosd(curr_lat_in_box(:,2:length(in_box_j_ind)) - (diff(curr_lat_in_box,1,2)/2)).*(mod(diff(curr_lon_in_box,1,2) + 180,360) - 180)) + (1i*diff(curr_lat_in_box,1,2)));
            delta_x = mean(diff_x_dist(isnan(diff_x_dist) == 0));
            delta_y = mean(diff_y_dist(isnan(diff_y_dist) == 0));
                    
            
            if length(size(unfilt_var_in_box)) > 2
                size_array_nested = [2 2 1].*size(unfilt_var_in_box);
            else
                size_array_nested = [2 2].*size(unfilt_var_in_box);
            end
            
            if n_iter_coastbuffer < 1
                curr_var_in_box = unfilt_var_in_box;
            else
                % define buffer over land areas near coast, using 2-D low-pass filter
                
                steepness_factor_buffer = min([steepness_factor 2]);
                
                curr_half_power_adj = exp(erfinv((2^(1/2)) - 1)/steepness_factor_buffer);   % adjustment factor to set bounds at half-power (rather than half-amplitude)
                
                % find planar fit for each level
                
                curr_var_planar_fit = NaN(size(unfilt_var_in_box));
                for z = 1:size(curr_var_planar_fit,3)
                    curr_array = unfilt_var_in_box(:,:,z);
                    good_ind = find(~isnan(curr_array));
                    if length(good_ind) < 0.2*numel(curr_array)
                        G = ones([numel(curr_array) 1]);
                        G_good = G(good_ind);
                        m = (((G_good')*G_good)^(-1))*((G_good')*curr_array(good_ind));
                    else
                        G = [ones([numel(curr_array) 1]) repmat((0:1:(size(curr_array,1) - 1))',[size(curr_array,2) 1]) reshape(repmat(0:1:(size(curr_array,2) - 1),[size(curr_array,1) 1]),[numel(curr_array) 1])];
                        G_good = G(good_ind,:);
                        m = (((G_good')*G_good)^(-1))*((G_good')*curr_array(good_ind));
                    end
                    curr_var_planar_fit(:,:,z) = reshape(G*m,[size(curr_array) 1]);
                end
                
                % iterate to produce non-tapered values near edges of data
                
                curr_var_in_box = unfilt_var_in_box - curr_var_planar_fit;
                curr_var_in_box((isnan(curr_var_in_box)) | (abs(curr_wet_mask_in_box) < 1e-5)) = 0;
                
                i_ind_nested = floor(size(curr_var_in_box,1)/2) + (1:1:size(curr_var_in_box,1))';
                j_ind_nested = floor(size(curr_var_in_box,2)/2) + (1:1:size(curr_var_in_box,2))';
                k_array = repmat((1/(delta_x*size_array_nested(1)))*(mod((0:1:(size_array_nested(1) - 1))' + (size_array_nested(1)/2),size_array_nested(1)) - (size_array_nested(1)/2)),[1 size_array_nested(2)]);
                l_array = repmat((1/(delta_y*size_array_nested(2)))*(mod((0:1:(size_array_nested(2) - 1)) + (size_array_nested(2)/2),size_array_nested(2)) - (size_array_nested(2)/2)),[size_array_nested(1) 1]);
                mag_waven_array = ((k_array.^2) + (l_array.^2)).^(1/2);
                low_bound = waven_x_low_bound/curr_half_power_adj;
                high_bound = waven_x_high_bound*curr_half_power_adj;
                bandpass_filter = 0.5*(erf(steepness_factor_buffer*(log(abs(mag_waven_array)) - log(low_bound))) - erf(steepness_factor_buffer*(log(abs(mag_waven_array)) - log(high_bound))));
                for iter = 1:n_iter_coastbuffer
                    curr_array = curr_var_in_box;
                    curr_array_nested = zeros(size_array_nested);
                    curr_array_nested(i_ind_nested,j_ind_nested,:) = curr_array;
                    fft_2d_curr_array = fft(fft(curr_array_nested,[],2),[],1);
                    
                    curr_array_filtered = real(ifft(ifft(repmat(bandpass_filter,[1 1 size(curr_array,3)]).*fft_2d_curr_array,[],2),[],1));
                    curr_var_in_box = curr_array_filtered(i_ind_nested,j_ind_nested,:);
        
                    curr_var_in_box(abs(curr_wet_mask_in_box - 1) < 1e-5) = unfilt_var_in_box(abs(curr_wet_mask_in_box - 1) < 1e-5) - curr_var_planar_fit(abs(curr_wet_mask_in_box - 1) < 1e-5);
                end
                
                curr_var_in_box = curr_var_in_box + curr_var_planar_fit;
            end
            
            
            % apply 2-D low pass filter after buffering
            
            half_power_adj = exp(erfinv((2^(1/2)) - 1)/steepness_factor);   % adjustment factor to set bounds at half-power (rather than half-amplitude)
            
            i_ind_nested = floor(size(curr_var_in_box,1)/2) + (1:1:size(curr_var_in_box,1))';
            j_ind_nested = floor(size(curr_var_in_box,2)/2) + (1:1:size(curr_var_in_box,2))';
            k_array = repmat((1/(delta_x*size_array_nested(1)))*(mod((0:1:(size_array_nested(1) - 1))' + (size_array_nested(1)/2),size_array_nested(1)) - (size_array_nested(1)/2)),[1 size_array_nested(2)]);
            l_array = repmat((1/(delta_y*size_array_nested(2)))*(mod((0:1:(size_array_nested(2) - 1)) + (size_array_nested(2)/2),size_array_nested(2)) - (size_array_nested(2)/2)),[size_array_nested(1) 1]);
            mag_waven_array = ((k_array.^2) + (l_array.^2)).^(1/2);
            low_bound = waven_x_low_bound/half_power_adj;
            high_bound = waven_x_high_bound*half_power_adj;
            bandpass_filter = 0.5*(erf(steepness_factor*(log(abs(mag_waven_array)) - log(low_bound))) - erf(steepness_factor*(log(abs(mag_waven_array)) - log(high_bound))));
            
            
            curr_array_planar_fit = NaN(size(curr_var_in_box));
            for z = 1:size(curr_array_planar_fit,3)
                curr_array = curr_var_in_box(:,:,z);
                good_ind = find(~isnan(curr_array));
                if length(good_ind) < 0.2*numel(curr_array)
                    G = ones([numel(curr_array) 1]);
                    G_good = G(good_ind);
                    m = (((G_good')*G_good)^(-1))*((G_good')*curr_array(good_ind));
                else
                    G = [ones([numel(curr_array) 1]) repmat((0:1:(size(curr_array,1) - 1))',[size(curr_array,2) 1]) reshape(repmat(0:1:(size(curr_array,2) - 1),[size(curr_array,1) 1]),[numel(curr_array) 1])];
                    G_good = G(good_ind,:);
                    m = (((G_good')*G_good)^(-1))*((G_good')*curr_array(good_ind));
                end
                curr_array_planar_fit(:,:,z) = reshape(G*m,[size(curr_array) 1]);
            end
            curr_array = curr_var_in_box - curr_array_planar_fit;
            curr_array((isnan(curr_array)) | (abs(curr_wet_mask_in_box) < 1e-5)) = 0;
            curr_array_nested = zeros(size_array_nested);
            curr_array_nested(i_ind_nested,j_ind_nested,:) = curr_array;
            fft_2d_curr_array = fft(fft(curr_array_nested,[],2),[],1);
            curr_array_filtered = real(ifft(ifft(repmat(bandpass_filter,[1 1 size(curr_array,3)]).*fft_2d_curr_array,[],2),[],1));
            curr_var_in_box_lp = curr_array_filtered(i_ind_nested,j_ind_nested,:) + curr_array_planar_fit;
            
            
            curr_array_planar_fit = NaN(size(curr_wet_mask_in_box));
            for z = 1:size(curr_array_planar_fit,3)
                curr_array = curr_wet_mask_in_box(:,:,z);
                good_ind = find(~isnan(curr_array));
                if length(good_ind) < 0.2*numel(curr_array)
                    G = ones([numel(curr_array) 1]);
                    G_good = G(good_ind);
                    m = (((G_good')*G_good)^(-1))*((G_good')*curr_array(good_ind));
                else
                    G = [ones([numel(curr_array) 1]) repmat((0:1:(size(curr_array,1) - 1))',[size(curr_array,2) 1]) reshape(repmat(0:1:(size(curr_array,2) - 1),[size(curr_array,1) 1]),[numel(curr_array) 1])];
                    G_good = G(good_ind,:);
                    m = (((G_good')*G_good)^(-1))*((G_good')*curr_array(good_ind));
                end
                curr_array_planar_fit(:,:,z) = reshape(G*m,[size(curr_array) 1]);
            end
            curr_array = curr_wet_mask_in_box - curr_array_planar_fit;
            curr_array((isnan(curr_array)) | (abs(curr_wet_mask_in_box) < 1e-5)) = 0;
            curr_array_nested = zeros(size_array_nested);
            curr_array_nested(i_ind_nested,j_ind_nested,:) = curr_array;
            fft_2d_curr_array = fft(fft(curr_array_nested,[],2),[],1);
            curr_array_filtered = real(ifft(ifft(repmat(bandpass_filter,[1 1 size(curr_array,3)]).*fft_2d_curr_array,[],2),[],1));
            curr_wet_mask_spat_lp = curr_array_filtered(i_ind_nested,j_ind_nested,:) + curr_array_planar_fit;
            
    
            cum_spat_weighting_array(in_box_i_ind,in_box_j_ind) = cum_spat_weighting_array(in_box_i_ind,in_box_j_ind) + curr_spat_weighting_array;
            curr_var_lp(in_box_i_ind,in_box_j_ind,:) = curr_var_lp(in_box_i_ind,in_box_j_ind,:) + ((repmat(curr_spat_weighting_array,[1 1 size(curr_var,3)])).*curr_var_in_box_lp);
            wet_mask_spat_lp(in_box_i_ind,in_box_j_ind,:) = wet_mask_spat_lp(in_box_i_ind,in_box_j_ind,:) + ((repmat(curr_spat_weighting_array,[1 1 size(wet_mask_array,3)])).*curr_wet_mask_spat_lp);
            
        end
    
        disp(['Completed latitude row number ',num2str(lat_center_ind),' (out of ',num2str(length(lat_center_decomp)),')'])
        
    end
    
    curr_var_lp = curr_var_lp./(repmat(cum_spat_weighting_array,[1 1 size(curr_var_lp,3)]));
    wet_mask_spat_lp = wet_mask_spat_lp./(repmat(cum_spat_weighting_array,[1 1 size(wet_mask_spat_lp,3)]));
    
    curr_var_lp(abs(wet_mask_array) < 1e-5) = orig_var(abs(wet_mask_array) < 1e-5);
    
end
