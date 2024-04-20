function [cost_neig, conf_score, dispar_int_val,conf_score_no_suprress] = disparCost(im_L, im_R, params)
% Calculate minus dispairty cost - this means we are calculating max
% instead of min
%% Output:
% cost_neig - cost containing min cost and its two neighboors
% conf_score - confidence score used to scale parabolas.


gaussKerSigma = params.gaussKerSigma;

% Some parameters/hyper-params
defval = -1e-4; % Make sure we don't divide by zero
mx = 2.2; % Found emperically
mn = 1; % Found empeariaclly 
conf_max = 1; % Limit confidence score to max of 1
conf_min = 0.01; % Limit confidence score to min of 0.01

% Calculate scores
idx_score = 1;
if strcmp(params.cost, 'SAD')  % Weighted SAD 
    for idx_dispar = params.dispar_vals
        % Note: Check if weighted can be optimized - instead of Gaussian
        score(:,:,idx_score) = imgaussfilt(-abs(im_L - imtranslate(im_R,[-idx_dispar,0],'linear','OutputView' ,'same','FillValues',0)),params.gaussKerSigma,'Padding','symmetric');
        idx_score = idx_score + 1;
    end
elseif strcmp(params.cost, 'SAD_win')  % Weighted SAD 
    h = ones(gaussKerSigma)/gaussKerSigma^2;

    for idx_dispar = params.dispar_vals
        % Note: Check if weighted can be optimized - instead of Gaussian
%         score(:,:,idx_score) = imgaussfilt(-abs(im_L - imtranslate(im_R,[-idx_dispar,0],'linear','OutputView' ,'same','FillValues',0)),params.gaussKerSigma,'Padding','symmetric');
        score(:,:,idx_score) = imfilter(-abs(im_L - imtranslate(im_R,[-idx_dispar,0],'linear','OutputView' ,'same','FillValues',0)),h,'same');

        idx_score = idx_score + 1;
    end
elseif strcmp(params.cost, 'SSD')  % Weighted Correlation 
    for idx_dispar = params.dispar_vals
        score(:,:,idx_score) = imgaussfilt(-abs(im_L - imtranslate(im_R,[-idx_dispar,0],'linear','OutputView' ,'same','FillValues',0)).^2,params.gaussKerSigma,'Padding','symmetric');
        idx_score = idx_score + 1;
    end
elseif strcmp(params.cost, 'CC')  % Weighted Correlation 
    for idx_dispar = params.dispar_vals
        score(:,:,idx_score) = imgaussfilt(im_L.*imtranslate(im_R,[-idx_dispar,0],'linear','OutputView' ,'same','FillValues',0),params.gaussKerSigma,'Padding','symmetric');
        idx_score = idx_score + 1;
    end
elseif strcmp(params.cost, 'NCC')  % Weighted Correlation 
    for idx_dispar = params.dispar_vals
        im_R_shift = imtranslate(im_R,[-idx_dispar,0],'linear','OutputView' ,'same','FillValues',0);
        ene_norm = sqrt(imgaussfilt(im_R_shift.*im_R_shift,params.gaussKerSigma, 'Padding','symmetric').* imgaussfilt(im_L.*im_L,params.gaussKerSigma, 'Padding','symmetric'));
        score(:,:,idx_score) = imgaussfilt(im_L.*im_R_shift,params.gaussKerSigma,'Padding','symmetric')./ene_norm;
        idx_score = idx_score + 1;
    end
    
% elseif strcmp(params.cost, 'NCC')  % Weighted - Normalized cross-correlation 
%     sizeKer = 2*gaussKerSigma;
%     avgFilt = ones(sizeKer)/sizeKer^2;
%     
%     im_L_ene = imfilter(im_L.*im_L,avgFilt,'same','symmetric');
%     for idx_dispar = params.dispar_vals
%         % Note: Check if weighted can be optimized - instead of Gaussian
%         im_R_shift = imtranslate(im_R,[-idx_dispar,0],'linear','OutputView' ,'same','FillValues',0);
% %         im_R_shift_ene = im_R_shift.*im_R_shift;
%         im_R_shift_ene = imfilter(im_R_shift.*im_R_shift, avgFilt,'same','symmetric');
%         normal_val = max(sqrt(im_R_shift_ene.*im_L_ene),abs(defval));
%         
% %         score(:,:,idx_score) = imgaussfilt(im_L.*imtranslate(im_R,[-idx_dispar,0],'linear','OutputView' ,'same','FillValues',0),params.gaussKerSigma,'Padding','symmetric')./normal_val;
%         score(:,:,idx_score) = imfilter(im_L.*im_R_shift,avgFilt,'same','symmetric')./normal_val;
%         idx_score = idx_score + 1;
%     end
 elseif strcmp(params.cost, 'ZNCC')  % Weighted - Normalized cross-correlation 
    im_L_mean = imgaussfilt(im_L, gaussKerSigma,'Padding','symmetric');
    im_L_ene = imgaussfilt((im_L-im_L_mean).*(im_L-im_L_mean),gaussKerSigma,'Padding','symmetric');
    for idx_dispar = params.dispar_vals
        % Note: Check if weighted can be optimized - instead of Gaussian
        im_R_shift = imtranslate(im_R,[-idx_dispar,0],'linear','OutputView' ,'same','FillValues',0);
%         im_R_shift_ene = im_R_shift.*im_R_shift;
        im_R_shift_mean = imgaussfilt(im_R, gaussKerSigma,'Padding','symmetric');   
        im_R_shift_ene = imgaussfilt((im_R_shift-im_R_shift_mean).*(im_R_shift-im_R_shift_mean), params.gaussKerSigma,'Padding','symmetric');
        normal_val = max(sqrt(im_R_shift_ene.*im_L_ene),abs(defval));
        
        score(:,:,idx_score) = imgaussfilt((im_L-im_L_ene).*imtranslate((im_R-im_R_shift_mean),[-idx_dispar,0],'linear','OutputView' ,'same','FillValues',0),params.gaussKerSigma,'Padding','symmetric')./normal_val;
        idx_score = idx_score + 1;
    end   
elseif strcmp(params.cost, 'ENCC')  % Weighted SAD 
    for idx_dispar = params.dispar_vals
        % Note: Check if weighted can be optimized - instead of Gaussian
        score(:,:,idx_score) = imgaussfilt(-abs(im_L - imtranslate(im_R,[-idx_dispar,0],'linear','OutputView' ,'same','FillValues',0)),params.gaussKerSigma,'Padding','symmetric');
        idx_score = idx_score + 1;
    end
elseif strcmp(params.cost, 'BT')
    h = ones(gaussKerSigma)/gaussKerSigma^2;
    im_L_minus = (im_L+imtranslate(im_L,[-1,0],'linear','OutputView' ,'same','FillValues',0))/2;
    im_L_plus = (im_L+imtranslate(im_L,[1,0],'linear','OutputView' ,'same','FillValues',0))/2;
    im_L_min = min(cat(3,im_L_minus,im_L,im_L_plus),[],3);
    im_L_max = max(cat(3,im_L_minus,im_L,im_L_plus),[],3);
    
    im_R_minus = (im_R+imtranslate(im_R,[-1,0],'linear','OutputView' ,'same','FillValues',0))/2;
    im_R_plus = (im_R+imtranslate(im_R,[1,0],'linear','OutputView' ,'same','FillValues',0))/2;
    im_R_min = min(cat(3,im_R_minus,im_R,im_R_plus),[],3);
    im_R_max = max(cat(3,im_R_minus,im_R,im_R_plus),[],3);
    
    for idx_dispar = params.dispar_vals
        shift_R_min = imtranslate(im_R_min,[-idx_dispar,0],'linear','OutputView' ,'same','FillValues',0);
        shift_R_max = imtranslate(im_R_max,[-idx_dispar,0],'linear','OutputView' ,'same','FillValues',0);
        shift_R = imtranslate(im_R,[-idx_dispar,0],'linear','OutputView' ,'same','FillValues',0);

        A1 = im_L-shift_R_max;
        A2 = shift_R_min-im_L;
        A = max(cat(3,zeros(size(im_R_min)),A1,A2),[],3);
        
        B1 = shift_R-im_L_max;
        B2 = im_L_min-shift_R;
        B = max(cat(3,zeros(size(im_R_min)),B1,B2),[],3);

        score(:,:,idx_score) = imfilter(min(cat(3,A,B),[],3),h,'replicate','same');
        idx_score = idx_score + 1;

    end
    score = -score;
    
end

% Get maximum score (remove disparity values from edges)
[max_score, max_idx] = maxk(score(:,:,2:end-1),2,3); 
max_idx = max_idx+1; % Add 1 back to compinsate for removing edges 

% Scale parabaolas depending on values of second maximum
rat_score = max_score(:,:,2) ./ min(max_score(:,:,1), defval);
a = 1/(mx-mn); b = - a*mn;
conf_score = max(min(a*rat_score + b, conf_max), conf_min).^2;
conf_score_no_suprress{2} = conf_score;
conf_score(abs(diff(max_idx,[],3))<=1) = 1; % If second max is close to first max don't reduce confidence
conf_score_no_suprress{1} = conf_score;
% 
if 1 % This is an alternative method - as described by google - check if there is second close score
    r_0 = 0.6;
    r_1 = 0.8;
    eps_d = 0.5;
    w_r = 3;
    d_i1 = -max_score(:,:,1); %Add here "-" because I computed -SSD
    d_i2 = -max_score(:,:,2);

    in1 = max(d_i1,eps_d);
    in1 = (in1-d_i2*r_0).^2;
    in1 = in1./(d_i2.^2*(r_1-r_0)^2);
    in1 = max(0,in1);
    in2 = min(in1,1);
    conf_score_no_suprress{3} = exp(-w_r*in2);
end


%% Look for places where minimum is not distinct 
% p1 = 0.25;
% p2 = 1;
% score_diff = abs(max_score(:,:,2) - max_score(:,:,1));
% disp_diff = abs(max_idx(:,:,1) - max_idx(:,:,2));
% % disp_diff_1 = disp_diff==1;
% disp_diff_2 = disp_diff>1;
% 
% % idx_penalty_1 = disp_diff_1&(score_diff<p1);
% idx_penalty_2 = disp_diff_2&(score_diff<p2);
% 
% % conf_score = ones(size(disp_diff));
% % conf_score(idx_penalty_1) = 1/100;
% conf_score(idx_penalty_2) = conf_score(idx_penalty_2)*1/1000;


%cost_neig = score(:,:, max_idx(:,:,1)-1:max_idx(:,:,1)+1); % This doesn't work 

% Only need to pass costs of max+2 neighboors
[xs,ys]=meshgrid(1:size(score,2),1:size(score,1));
cost_neig = zeros(size(score,1),size(score,2),3);
cost_neig(:,:,1) = score(sub2ind(size(score),ys,xs,max_idx(:,:,1)-1));
cost_neig(:,:,2) = score(sub2ind(size(score),ys,xs,max_idx(:,:,1)));
cost_neig(:,:,3) = score(sub2ind(size(score),ys,xs,max_idx(:,:,1)+1));


dispar_int_val = params.dispar_vals(max_idx(:,:,1));

end


