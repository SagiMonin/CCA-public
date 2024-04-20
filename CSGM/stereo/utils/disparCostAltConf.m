function [cost_neig, conf_score, dispar_int_val] = disparCostAltConf(im_L, im_R, params,im_guide)
% Calculate minus dispairty cost - this means we are calculating max
% instead of min
%% Output:
% cost_neig - cost containing min cost and its two neighboors
% conf_score - confidence score used to scale parabolas.


gaussKerSigma = params.gaussKerSigma;
% Some parameters/hyper-params
defval = -1e-10; % Make sure we don't divide by zero
mx = 2.2; % Found emperically
mn = 1; % Found empeariaclly 
a = 1/(mx-mn); b = - a*mn;
conf_max = 1; % Limit confidence score to max of 1
conf_min = 1e-2; % Limit confidence score to min of 0.01
h =  ones(params.gaussKerSigma)/params.gaussKerSigma.^2;
% h = fspecial('average',params.gaussKerSigma);

% Calculate scores
idx_score = 1;
if strcmp(params.cost, 'SAD')  % Weighted SAD 
    for idx_dispar = params.dispar_vals
        % Note: Check if weighted can be optimized - instead of Gaussian
%         score(:,:,idx_score) = imgaussfilt(-abs(im_L - imtranslate(im_R,[-idx_dispar,0],'linear','OutputView' ,'same','FillValues',0)),params.gaussKerSigma,'Padding','symmetric');
        score(:,:,idx_score) = imfilter(-abs(im_L - imtranslate(im_R,[-idx_dispar,0],'linear','OutputView' ,'same','FillValues',0)),h,'replicate','same');
        idx_score = idx_score + 1;
    end
    
elseif strcmp(params.cost, 'BT')
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
[max_score, max_idx] = maxk(score(:,:,2:end-1),5,3); 
max_idx = max_idx+1; % Add 1 back to compinsate for removing edges 

dispar_int_val = params.dispar_vals(max_idx(:,:,1));


GT = evalin('base', 'GT');

% Scale parabola if second minimum is not adjacent - this score is not good for large disparity (if I have a second minimum which is far and not that large it can penalize)
rat_score1 = max_score(:,:,2) ./ min(max_score(:,:,1), defval);
conf_score1 = max(min(a*rat_score1 + b, conf_max), conf_min).^2;
conf_score1(abs(max_idx(:,:,1)-max_idx(:,:,2))<=1) = 1; % If second max is close to first max don't reduce confidence

rat_score2 = max_score(:,:,3) ./ min(max_score(:,:,1), defval);
conf_score2 = max(min(a*rat_score2 + b, conf_max), conf_min).^2;
conf_score2(abs(max_idx(:,:,1)-max_idx(:,:,3))<=1) = 1; % If second max is close to first max don't reduce confidence

rat_score3 = max_score(:,:,4) ./ min(max_score(:,:,1), defval);
conf_score3 = max(min(a*rat_score3 + b, conf_max), conf_min).^2;
conf_score3(abs(max_idx(:,:,1)-max_idx(:,:,4))<=1) = 1; % If second max is close to first max don't reduce confidence

rat_score4 = max_score(:,:,5) ./ min(max_score(:,:,1), defval);
conf_score4 = max(min(a*rat_score4 + b, conf_max), conf_min).^2;
conf_score4(abs(max_idx(:,:,1)-max_idx(:,:,5))<=1) = 1; % If second max is close to first max don't reduce confidence


conf_score = min(cat(3,conf_score1,conf_score2,conf_score3,conf_score4),[],3);



% conf_score(abs(diff(max_idx,[],3))<=1) = 1; % If second max is close to first max don't reduce confidence
% 
%% Look for places where minimum is not distinct 
% look_val = gaussKerSigma; % Distance to look at from minimum in order to say if confident in score
% rat_val = 0.5; % Ratio value to look if ratio is smaller than. Higher will result in less confident pixels
[score_sort,idx_sort ]= sort(score,3,'descend'); % Sort scores
% 
val_disp = 0.1;
rat_disp = 2;
idx_disp_large = abs(idx_sort(:,:,2:5)-idx_sort(:,:,1)) > rat_disp;
idx_disp_score = abs(score_sort(:,:,2:5)-score_sort(:,:,1)) < val_disp;
idx_remove = any(idx_disp_large&idx_disp_score,3);
conf_score(idx_remove) = conf_min.^2;

% figure; imagesc(idx_disp_large&idx_disp_score)
% rat_score_other = (score_sort(:,:,2:end)./score_sort(:,:,1))<rat_val; % Check if ratio is is less than rat_val
% idx_to_close = (idx_sort(:,:,2:end).*rat_score_other)>look_val; % If dispariy value is under rat_val and it is at least look_val pixels away, mark it as not confident
% tmp_conf = sum(idx_to_close,3);

% conf_score = conf_score .* exp(-tmp_conf);
% 
% conf_score = conf_score + 1e-100;
% 
% %% Look for places where minimum is not distinct 
% % p1 = 0.25;
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
% conf_score(idx_penalty_2) = conf_score(idx_penalty_2)*1/100;
% 
% 
% 
% % 
%cost_neig = score(:,:, max_idx(:,:,1)-1:max_idx(:,:,1)+1); % This doesn't work 

if 0
    figure(999);
    ax(1) = subplot(231); imagesc(conf_score1);
    ax(2) = subplot(232); imagesc(conf_score2);
    ax(3) = subplot(233); imagesc(-dispar_int_val);
    ax(4) = subplot(234); imagesc(conf_score);
    ax(5) = subplot(235); imagesc(GT);
    ax(6) = subplot(236); imagesc(abs(GT+dispar_int_val)<1);
    linkaxes(ax);
end

% Only need to pass costs of max+2 neighboors
[xs,ys]=meshgrid(1:size(score,2),1:size(score,1));
cost_neig = zeros(size(score,1),size(score,2),3);
cost_neig(:,:,1) = score(sub2ind(size(score),ys,xs,max_idx(:,:,1)-1));
cost_neig(:,:,2) = score(sub2ind(size(score),ys,xs,max_idx(:,:,1)));
cost_neig(:,:,3) = score(sub2ind(size(score),ys,xs,max_idx(:,:,1)+1));



end



