function conf = calc_conf_dslr(sum_parab,im_guide,dispar_total,params)
%% Calculate Confidecnec 
% Step 1 based in quadratic value and image value
sig_a = 5;
sig_c = 512;
% a = (left_val+right_val)/2; % - on both scores since they are for max

% c1 = exp( log(abs(sum_parab.a))/sig_a^2 + sum_parab.c/sig_c^2); % - on score since it is for max : exp(
c1 = exp( log(abs(sum_parab.a))/sig_a^2); % - on score since it is for max : exp(

% Step 2: Check neighbors if there is a large difference remove
sig_u = 1;

idx_shift = 0;
for idx_x = -1:1
    for idx_y = -1:1
        idx_shift = idx_shift + 1;
        dispar_total_shift(:,:,idx_shift) = imtranslate(dispar_total,[idx_x,idx_y]);
    end
end
dispar_total_shift(:,:,5) = inf; % Dont consider zero translation

in2 = min((dispar_total-dispar_total_shift).^2/sig_u.^2 ,[], 3);
c2= exp(-in2);

% Step 3: Check horizontal neighbors 
w_v = 125;
e_v = 25;

sobel_filt = [-1 0 1; -2 0 2; -1 0 1];
grads = abs(imfilter(rgb2gray(im_guide/255)*255,sobel_filt,'same','replicate'));    
grads_filt = imgaussfilt(grads,params.gaussKerSigma);
in3 = max(0,w_v./grads_filt-e_v);

c3 = exp(-in3);

% conf = c1.*c2.*c3;
conf = c1.*c2.*c3;

end

