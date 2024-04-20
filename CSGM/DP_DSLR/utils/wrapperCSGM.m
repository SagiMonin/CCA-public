function [dispar_map_pyr, sum_parab,conf_score_no_suprress] = wrapperCSGM(im_L, im_R, im_guide, params )
% Wrapper function to run CSGM
debug_flag = 1;

if params.levels == 1
    [dispar_map, sum_parab,conf_score_no_suprress] = currLvlCSGM(im_L, im_R, im_guide, params, [], []);
    dispar_map_pyr{1} = dispar_map;
elseif params.levels > 1
    % Generate pyramid images - coarse to fine (level 1 is low res) 
    [im_blur_L, im_diff_L, im_pyr_L] = pyrImg(im_L, params.levels); 
    [im_blur_R, im_diff_R, im_pyr_R] = pyrImg(im_R, params.levels);
    [~,~,im_pyr_guide] = pyrImg(im_guide, params.levels);
    params.levels = size(im_pyr_R,2); % If levels is larger than possible number of levels.

    % First step on lowest res
    [dispar_map_pyr{1}, sum_parab,conf_score_no_suprress] = currLvlCSGM(im_pyr_L{1}, im_pyr_R{1}, im_pyr_guide{1}, params, [], []);
    if params.plt_pyr == 1
        figure(201);
        subplot(121); imshow(dispar_map_pyr{1},[]); title(['Disparity map - pyramid lvl: 1']); colorbar;
        subplot(122); imshow(im_pyr_guide{1},[]); title(['Guide image - pyramid lvl: 1']); colorbar;
        pause(1e-6);
    end
    for idx_pyr = 2:params.levels
        params.idx_pyr = idx_pyr;
        curr_lvl_size = size(im_pyr_L{idx_pyr});
        prev_lvl_size = size(im_pyr_L{idx_pyr-1});
        
        % Calculate rescale factor: this is usually 2, but maybe because of
        % rounding a bit differemt
        rescale_fact = curr_lvl_size(1)./prev_lvl_size(1);

        % Bilinear im-resizing to current size
        dispar_init = imresize(dispar_map_pyr{idx_pyr-1}, curr_lvl_size, 'bilinear'); 
        dispar_init = dispar_init*rescale_fact;
        
        % resize the parabolas and use as prior
        % Update the parabolas according to a.*(d/rescaleFac).^2+b.*(d/rescaleFac)+c
        sum_parab.a = sum_parab.a .* (1/rescale_fact)^2;
        sum_parab.b = sum_parab.b .* (1/rescale_fact);
        
        sum_parab.a = imresize(sum_parab.a, curr_lvl_size, 'bilinear');
        sum_parab.b = imresize(sum_parab.b, curr_lvl_size, 'bilinear');
        sum_parab.c = imresize(sum_parab.c, curr_lvl_size, 'bilinear');       
        
        % Perform next step
        [dispar_map_pyr{idx_pyr}, sum_parab,conf_score_no_suprress] = currLvlCSGM(im_pyr_L{idx_pyr}, im_pyr_R{idx_pyr}, im_pyr_guide{idx_pyr}, params, sum_parab, dispar_init);
        
        if params.plt_pyr == 1
            figure(200+idx_pyr);
            subplot(121); imshow(dispar_map_pyr{idx_pyr},[]); title(['Disparity map - pyramid lvl: ', num2str(idx_pyr)]); colorbar;
            subplot(122); imshow(im_pyr_guide{idx_pyr},[]); title(['Guide image - pyramid lvl: ', num2str(idx_pyr)]); colorbar;
            pause(1e-6);
        end

        
    end
end 






% %% Pre-processing - generate pyramids for multi-level CSGM
% 
% % Create an image pyramid




end

