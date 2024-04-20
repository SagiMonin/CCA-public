function [aiwe1val, aiwe2val, srccval, geomean] = cmp_all_metrics(depth_estimate, depth_GT)
% Calculate metrix to estimate disparity. 
crop_size = 100; % Crop because of GT of DSLR

gt_crop = depth_GT(crop_size:end-crop_size,crop_size:end-crop_size);
est_crop = depth_estimate - min(depth_estimate(:));
est_crop = floor(est_crop/max(est_crop(:)) * 255)/255;  % This step quanitzes to 256 levels so we can compare to spearman 
est_crop = est_crop(crop_size:end-crop_size,crop_size:end-crop_size);

% Note I checked quantization step - quantized to 512 and 1024 and it
% didn't effect scores

aiwe1val = aiwe1_calc(est_crop,gt_crop);

[aiwe2val,a,b] = aiwe2_calc(est_crop,gt_crop);

srccval = srcc_calc(est_crop,gt_crop);

geomean = (srccval*aiwe1val*aiwe2val)^(1/3);
fprintf('%.3f %.3f %.3f %.3f \n',aiwe1val, aiwe2val, srccval, geomean); 
end

