function [cost_neig, conf_score, dispar_int_val] = calc_init_dispar(im_L,im_R,params,x_range,y_range,reduce_fact,y_flag)
% Calculate interger disparity based on two images:
% Input: im_L,im_R - images
%        params - params used for disparity cost
%        x_range + y_range - disparity ranges for x+y
%        reduce_fact - reduction factor of image size.
%        y_flag - calculate y-disparity
% Output: cost_neig - cost in neighbhoorhood of minmum
%         conf_score - confidence score
%         dispar_int_val - disparity integer value


if y_flag == 1
    params.dispRangeY = round(y_range*(1/reduce_fact)); % Disparity range - (this is for first level of pyramid) for DP this is symmetric 
    params.dispRangeX = round(x_range*(1/reduce_fact)); % Disparity range - (this is for first level of pyramid) for DP this is symmetric 
    params.dispar_vals_x = params.dispRangeX(1)-1:params.dispRangeX(end)+1; % Add -+1 values to disparity range for calculating parabolas    
    params.dispar_vals_y = params.dispRangeY(1)-1:params.dispRangeY(end)+1; % Add -+1 values to disparity range for calculating parabolas 
    
    [cost_neig, conf_score, dispar_int_val] = disparCostYdir(im_L, im_R, params);
else
    
    params.dispRangeX = round(x_range*(1/reduce_fact)); % Disparity range - (this is for first level of pyramid) for DP this is symmetric 
    params.dispar_vals_x = params.dispRangeX(1)-1:params.dispRangeX(end)+1; % Add -+1 values to disparity range for calculating parabolas    
    params.dispar_vals_y = 0;

    [cost_neig, conf_score, dispar_int_val] = disparCostYdir(im_L, im_R, params);
end

end

