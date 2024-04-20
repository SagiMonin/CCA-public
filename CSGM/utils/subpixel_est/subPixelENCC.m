function [tau_clean, invalid_mask,tau_dirty] = subPixelENCC(im_L, im_R, int_dispar_val, param, next_val_mat)
% Assume right image is translated to the left. Disparity is positive.
% Ouput: tau_clean - subpixel map
%        invalid_mask - invalid pixel masks.
% Input: Image left, image right, integer disparity value, param - window size


window_size = param.encc_window;
win_ceil = ceil(window_size/2); win_floor = floor(window_size/2);
[M,N] = size(im_L); % Assume gray-scale image

% Padding so we can calculate borders as-well.
max_dispar = max(abs(int_dispar_val(:)));

zero_pad_im_L = [zeros(M+2*win_floor,1+max_dispar) ,padarray(im_L, [win_floor, win_floor],0),zeros(M+2*win_floor,1)]; % Add extra column to left for disparity correction, add to the right because of interpolation to the right instead of to the left
zero_pad_im_R = [zeros(M+2*win_floor,1+max_dispar) ,padarray(im_R, [win_floor, win_floor],0),zeros(M+2*win_floor,1)]; % Add extra column to left
zero_pad_int_dispar_val = [zeros(M+2*win_floor,1+max_dispar) ,padarray(int_dispar_val, [win_floor, win_floor],0),zeros(M+2*win_floor,1)]; % Add extra column to left
zero_pad_next_val_mat = [zeros(M+2*win_floor,1+max_dispar) ,padarray(next_val_mat, [win_floor, win_floor],0),zeros(M+2*win_floor,1)]; % Add extra column to left

[M_pad,N_pad] = size(zero_pad_im_L);

for idx_x = max_dispar+win_ceil+1:N_pad-win_floor-1 %This is shifted by 
    tic
    for idx_y = win_ceil:M_pad-win_floor
        curr_dipar = zero_pad_int_dispar_val(idx_y, idx_x);
        curr_next_val = zero_pad_next_val_mat(idx_y, idx_x);
        
        % Get left window 
        win_L = zero_pad_im_L(idx_y-win_floor:idx_y+win_floor, idx_x-win_floor:idx_x+win_floor);
        win_L = win_L - mean(win_L(:)); % Subtract mean
        win_L_norm = win_L/sqrt((win_L(:)'*win_L(:))); % Normalize window 
        % Get right image window
        win_R = zero_pad_im_R(idx_y-win_floor:idx_y+win_floor, idx_x-win_floor-curr_dipar+curr_next_val:idx_x+win_floor-curr_dipar+curr_next_val);
        win_R = win_R - mean(win_R(:));
        win_R_norm = win_R/sqrt((win_R(:)'*win_R(:)));
        % Get right inage window to the left
        win_R_left = zero_pad_im_R(idx_y-win_floor:idx_y+win_floor, idx_x-win_floor-1-curr_dipar+curr_next_val:idx_x+win_floor-1-curr_dipar+curr_next_val); % Taking window to the left of one comparing 
        win_R_left = win_R_left - mean(win_R_left(:));
        win_R_left_norm = win_R_left/sqrt((win_R_left(:)'*win_R_left(:)));
        
        lambda = sqrt((win_R_left(:)'*win_R_left(:)))/sqrt((win_R(:)'*win_R(:)));
        r = win_R_norm(:)'*win_R_left_norm(:);
        rho_d = win_R_norm(:)'*win_L_norm(:);
        rho_d_1 = win_R_left_norm(:)'*win_L_norm(:);
        
        val = (rho_d_1 - r*rho_d)/(lambda*(r*rho_d_1 - rho_d) + r*rho_d - rho_d_1);
        tau(idx_y,idx_x) = val;
    end
    toc
end

% Remove border
tau_clean  = tau(win_ceil:end, win_ceil+1+max_dispar:end); 
tau_dirty = tau_clean;
% Get invalid values (nan, larger than 0, smaller than -1)
nan_mask = isnan(tau_clean);
too_small = tau_clean<-1;
too_big = tau_clean>0;
invalid_mask = nan_mask|too_small|too_big;

tau_clean(invalid_mask) = 0; % This might happen in the borders if disparity is larger then window


% % Main loop - Not taking care of borders 
% for idx_x = win_ceil+1:N-win_floor %This is shifted by 
%     for idx_y = win_ceil:M-win_floor
%         curr_dipar = int_dispar_val(idx_y, idx_x);
%         
%         % Get left window 
%         win_L = im_L(idx_y-win_floor:idx_y+win_floor, idx_x-win_floor:idx_x+win_floor);
%         win_L = win_L - mean(win_L(:)); % Subtract mean
%         win_L = win_L/(win_L(:)'*win_L(:)); % Normalize window 
%         % Get left window
%         win_R = im_R(idx_y-win_floor:idx_y+win_floor, idx_x-win_floor-curr_dipar:idx_x+win_floor-curr_dipar);
%         win_R_left = im_R(idx_y-win_floor:idx_y+win_floor, idx_x-win_floor-1-curr_dipar :idx_x+win_floor-1-curr_dipar); % Taking window to the left of one comparing 
%         win_R = win_R - mean(win_R(:));
%         win_R_left = win_R_left - mean(win_R_left(:));
%         win_R = win_R/(win_R(:)'*win_R(:));
%         win_R_left = win_R_left/(win_R_left(:)'*win_R_left(:));
%         
%         lambda = (win_R(:)'*win_R(:))/(win_R_left(:)'*win_R_left(:));
%         r = win_R(:)'*win_R_left(:);
%         rho_d = win_R(:)'*win_L(:);
%         rho_d_1 = win_R_left(:)'*win_L(:);
%         
%         tau(idx_y,idx_x) = (rho_d_1 - r*rho_d)/(lambda*(r*rho_d_1 - rho_d) + r*rho_d - rho_d_1);
%     end
% end


end

