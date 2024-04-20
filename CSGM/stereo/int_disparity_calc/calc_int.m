% Choose data set
imgset = 'training';
%imgset = 'test';

% Specify which resolution you are using for the stereo image set (F, H, or Q?)
imgsize = 'Q';
% imgsize = 'H';
% imgsize = 'F';

%
path_save = 'C:\Users\local_admin\Documents\GitHub\C-SGM\CSGM\stereo\int_disparity_calc';


for im_num = 1
    [im_L,im_R,GT,mask,ndisp] = load_midbury(im_num,imgsize,imgset);
%     im_L = imread(['MiddEval3/',imgset,imgsize,'/',image_names{im_num},'/im0.png']);
%     im_R = imread(['MiddEval3/',imgset,imgsize,'/',image_names{im_num},'/im1.png']);
    
    im_L = double(im_L);
    im_R = double(im_R);
    
    im_guide = im_L;
    im_L = rgb2gray(im_L/255)*255;
    im_R = rgb2gray(im_R/255)*255;
    
    % Adjust the range of disparities to the chosen resolution
    if imgsize == 'Q'
        DisparityRange = [-round(ndisp(im_num)/4),-1];
    elseif imgsize == 'H'
        DisparityRange = [-round(ndisp(im_num)/2),-1];
    else
        DisparityRange = [-round(ndisp(im_num)),-1];
    end
    
     im_L = im_L - cv.bilateralFilter(im_L ,'Diameter',15, 'SigmaSpace', 3, 'SigmaColor', 20);
     im_R = im_R - cv.bilateralFilter(im_R ,'Diameter',15, 'SigmaSpace', 3, 'SigmaColor', 20);
%      im_L = im_L -min(im_L(:));
%      im_R = im_R -min(im_R(:));
%     params.normKerSize = 7;5
%     [im_L,~]=normalizeImg(im_L,params.normKerSize);
%     [im_R,~]=normalizeImg(im_R,params.normKerSize);
    %% Calculate disparity 
    params.cost = 'BT';
    params.gaussKerSigma = 5;
    params.levels = 3;
    
    % Calculate pyramids
    [im_blur_L, im_diff_L, im_pyr_L] = pyrImg(im_L, params.levels); 
    [im_blur_R, im_diff_R, im_pyr_R] = pyrImg(im_R, params.levels);
    [~,~,im_pyr_guide] = pyrImg(im_guide, params.levels);
    tic;
    for ii = 1:params.levels
        params.dispar_vals = (floor(DisparityRange(1)/2.^(params.levels-ii))-1):(DisparityRange(2)+1);
%         [cost_neig{ii}, conf_score{ii}, dispar_int_val{ii},conf_score_no_suprress] = disparCost(im_pyr_L{ii},im_pyr_R{ii}, params);
        [cost_neig{ii}, conf_score{ii}, dispar_int_val{ii}] = disparCostAltConf(im_pyr_L{ii},im_pyr_R{ii}, params);
        params.dispar_vals = -fliplr(params.dispar_vals); 
%         [cost_neig_R{ii}, conf_score_R{ii}, dispar_int_val_R{ii},conf_score_no_suprress] = disparCost(im_pyr_R{ii},im_pyr_L{ii}, params);
        [cost_neig_R{ii}, conf_score_R{ii}, dispar_int_val_R{ii}] = disparCostAltConf(im_pyr_R{ii},im_pyr_L{ii}, params);

    end
    toc
    curr_save = [path_save,'\tgt',num2str(im_num),'_res',imgsize,'.mat'];

    save(curr_save, 'im_L', 'im_R', 'im_guide',...
        'cost_neig', 'conf_score', 'dispar_int_val', ...
        'cost_neig_R', 'conf_score_R', 'dispar_int_val_R');

end
