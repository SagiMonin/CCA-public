% Choose data set
imgset = 'training';
%imgset = 'test';

% Specify which resolution you are using for the stereo image set (F, H, or Q?)
imgsize = 'Q';
% imgsize = 'H';
% imgsize = 'F';

% Parameters to choose 
params.cost = 'BT'; % Cost - BT or SAD 
params.gaussKerSigma = 5; % Size of window of cost
path_save = 'C:\Users\local_admin\CSGM\CSGM_Monin\stereo\int_disparity_calc'; % Save path location


for im_num = 1:15
    [im_L,im_R,GT,mask,ndisp] = load_midbury(im_num,imgsize,imgset);
    
    im_L = double(im_L);
    im_R = double(im_R);
    
    im_guide = im_L;
    im_guide_R = im_R;

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
%     
     im_L = im_L - cv.bilateralFilter(im_L ,'Diameter',15, 'SigmaSpace', 3, 'SigmaColor', 20);
     im_R = im_R - cv.bilateralFilter(im_R ,'Diameter',15, 'SigmaSpace', 3, 'SigmaColor', 20);
%      
     
%     params.normKerSize = 6;
%     [im_L,~]=normalizeImg(im_L,params.normKerSize);
%     [im_R,~]=normalizeImg(im_R,params.normKerSize);

    %% Calculate disparity 
    tic
    params.dispar_vals = DisparityRange(1)-1:DisparityRange(2)+1;
    [cost_neig, conf_score, dispar_int_val] = disparCostAltConf(im_L, im_R, params, im_guide);
    params.dispar_vals = -fliplr(params.dispar_vals); 

    [cost_neig_R, conf_score_R, dispar_int_val_R] = disparCostAltConf(im_R,im_L, params, im_guide_R);

    
    toc
    curr_save = [path_save,'\tgt',num2str(im_num),'_res',imgsize,'.mat'];

    save(curr_save, 'im_L', 'im_R', 'im_guide',...
        'cost_neig', 'conf_score', 'dispar_int_val', ...
        'cost_neig_R', 'conf_score_R', 'dispar_int_val_R');

end
