%% ENNC
params = paramsCSGM();
params.encc_window = 3;

im_L = kron(rand(64),ones(4));
% im_L = im2double(imread('cameraman.tif'));
im_R = imtranslate(im_L, [-0.7,0],'linear','OutputView' ,'same','FillValues',0);
im_guide = im_L;
dispar_map = zeros(size(im_L))-1;

[D2, T2, W2] = enccMEX(im_L, im_R, [15,15], 3, 1, 0);


% [D2, T2, W2] = encc(im_L, im_R, [15,15], 3, 1, 0);
[ T3, W3,max_map] = encc_subpixel_alter2(im_L, im_R, [15,15], max(abs(dispar_map(:)))+1, 1, 0, dispar_map);
% figure; imagesc(D)
% tic
% [ T3, W3,max_map1] = encc_subpixel_alter3(im_L, im_R, [15,15], max(abs(dispar_map(:)))+1, 1, 0, dispar_map);
% toc
tic
[ T4,max_map2] = encc_subpixel_alter3_2(im_L, im_R, [15,15], max(abs(dispar_map(:)))+1, 1, 0, dispar_map);

toc

% final_dispar = dispar_map;
% final_dispar = zeros(size(max_map));
% final_dispar = dispar_map+double(max_map==2).*T3+double(max_map==3).*(T3+1);
final_dispar = dispar_map+(T3+1);

% final_dispar = dispar_map+double(max_map==3).*(T3+1);
figure(1);
subplot(221); imagesc(final_dispar+D2,[-1,1]);
subplot(222); imagesc(T4-T3);
subplot(223); imagesc(-D2);
subplot(224); imagesc(final_dispar)


%% This is how to compensate if the disparity is positive
im_L = kron(rand(64),ones(4));
% im_L = im2double(imread('cameraman.tif'));
im_R = imtranslate(im_L, [1.2,0],'linear','OutputView' ,'same','FillValues',0);
im_guide = im_L;
dispar_map = zeros(size(im_L))+1;

max_disp = max(dispar_map(:))+1;
dispar_map = dispar_map - max_disp;
im_R_shift = imtranslate(im_R, [-max_disp,0],'linear','OutputView' ,'same','FillValues',0);

[D2, T2, W2] = enccMEX(im_L, im_R_shift, [15,15], max_disp, 1, 0);
D2 = D2-max_disp;
% [D2, T2, W2] = enccMEX(im_R, im_L, [15,15], 4, 1, 0);

% [D2, T2, W2] = encc(im_L, im_R, [15,15], 4, 1, 0);
[ T3, W3,max_map] = encc_subpixel_alter2(im_L, im_R_shift, [15,15], max(abs(dispar_map(:)))+1, 1, 0, dispar_map);

% final_dispar = dispar_map;
% final_dispar = zeros(size(max_map));
final_dispar = dispar_map+double(max_map==2).*T3+double(max_map==3).*(T3+1) + max_disp;

figure(1);
subplot(221); imagesc(final_dispar+D2,[-1,1]);
subplot(222); imagesc(T2-T3);
subplot(223); imagesc(-D2);
subplot(224); imagesc(final_dispar)