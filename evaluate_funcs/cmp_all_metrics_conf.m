function [aiwe1val, aiwe2val, srccval, geomean] = cmp_all_metrics_conf(depth_estimate, depth_GT, conf,name_save,flag_save,idx_im)
% Calculate metrix to estimate disparity. 
crop_size = [512,384]; % Crop because of GT of DSLR

% Convert to 8-bit image
est_crop = double(cropCenter(depth_estimate, crop_size));
est_crop = est_crop-min(est_crop(:));
est_crop = est_crop/max(est_crop(:));
est_crop = double(uint8(est_crop*255))/255;

gt_crop = double(cropCenter(depth_GT, crop_size))/255;
conf_crop = double(cropCenter(conf, crop_size));

aiwe1val = aiwe1_calc_conf(est_crop,gt_crop,conf_crop);

aiwe2val = aiwe2_calc_conf(est_crop,gt_crop,conf_crop);

srccval = srcc_calc_conf(est_crop,gt_crop,conf_crop);
geomean = (srccval*aiwe1val*aiwe2val)^(1/3);

% fprintf('%.3f %.3f %.3f %.3f \n',aiwe1val, aiwe2val, srccval, geomean); 


if flag_save == 1
    name_save = [name_save,'\im',num2str(idx_im),'.mat'];
    est_int = depth_estimate;
    est_int = est_int-min(est_int(:));
    est_int = est_int/max(est_int(:));
    est_int = uint8(est_int*255);
    
    
%     est_crop = uint8(est_crop*255);
    save(name_save,'est_int');
end
end

