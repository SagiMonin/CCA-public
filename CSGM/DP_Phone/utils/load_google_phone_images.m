function imname = load_google_phone_images(path_im, num_im)

folder_info = dir(path_im);
perm_idxs = randperm(length(folder_info));
perm_idxs(perm_idxs == 1) = [];
perm_idxs(perm_idxs == 2) = [];
% perm_idxs = 3:52;
if num_im == 684
    perm_idxs = 3:num_im+2;
end
for idx_im = 1:num_im
    curr_idx = perm_idxs(idx_im);
    
    imname{idx_im} = [folder_info(curr_idx).name,'\'];
end
    
end

