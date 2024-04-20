function [im_l, im_r, im_guide, gt_depth, conf_depth] = loadPhone(pathData,imname)
%% Read all images and crop them
crop_l = [1344,1008];
crop_s = [672,504];

file_path_right = [pathData,'right_pd\',imname,'result_rightPd_center.png'];
im_r = imread(file_path_right);
im_r = double(cropCenter(im_r, crop_l));

file_path_left = [pathData,'left_pd\',imname,'result_leftPd_center.png'];
im_l = imread(file_path_left);
im_l = double(cropCenter(im_l, crop_l));

file_path_guide = [pathData,'scaled_images\',imname,'result_scaled_image_center.jpg'];
im_guide = imread(file_path_guide);
im_guide = double(cropCenter(im_guide, crop_s));
im_guide = imresize(im_guide,2);

file_path_depth = [pathData,'merged_depth\',imname,'result_merged_depth_center.png'];
gt_depth = imread(file_path_depth);
gt_depth = double(cropCenter(gt_depth, crop_s));

file_path_conf = [pathData,'merged_conf\',imname,'result_merged_conf_center.exr'];
conf_depth = cv.imread(file_path_conf, 'Unchanged', true, 'Color', false, 'Grayscale', true); conf_depth = conf_depth(:,:,1);
conf_depth = double(cropCenter(conf_depth, crop_s));

if 0
    figure(999);
    ax(1) = subplot(151); imagesc(im_r);
    ax(2) = subplot(152); imagesc(im_l);
    ax(3) = subplot(153); imagesc(im_guide/255);
    ax(4) = subplot(154); imagesc(gt_depth);
    ax(5) = subplot(155); imagesc(conf_depth);
    linkaxes(ax);

end

end

