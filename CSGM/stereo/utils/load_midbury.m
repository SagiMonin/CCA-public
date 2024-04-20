function [im_L,im_R,GT,mask,ndisp] = load_midbury(im_idx,imgsize,imgset)
% Choose data set
% imgset = 'training';
%imgset = 'test';

% Specify which resolution you are using for the stereo image set (F, H, or Q?)
% imgsize = 'Q';
% imgsize = 'H';
%imgsize = 'F';
% imgsize = res;
%
% path_save = 'C:\Users\local_admin\Documents\GitHub\C-SGM\Stereo_data_example\int_disparity_calc';

if ~exist('MiddEval3results','dir')
    mkdir('MiddEval3results');
end
if ~exist(['MiddEval3results/',imgset,imgsize],'dir')
    mkdir(['MiddEval3results/',imgset,imgsize]);
end

if strcmp(imgset,'training')
    image_names{1} = 'Adirondack';
    image_names{2} = 'ArtL';
    image_names{3} = 'Jadeplant';
    image_names{4} = 'Motorcycle';
    image_names{5} = 'MotorcycleE';
    image_names{6} = 'Piano';
    image_names{7} = 'PianoL';
    image_names{8} = 'Pipes';
    image_names{9} = 'Playroom';
    image_names{10} = 'Playtable';
    image_names{11} = 'PlaytableP';
    image_names{12} = 'Recycle';
    image_names{13} = 'Shelves';
    image_names{14} = 'Teddy';
    image_names{15} = 'Vintage';
    ndisp = [290, 256, 640, 280, 280, 260, 260, 300, 330, 290, 290, 260, 240, 256, 760];
else
    image_names{1} = 'Australia';
    image_names{2} = 'AustraliaP';
    image_names{3} = 'Bicycle2';
    image_names{4} = 'Classroom2';
    image_names{5} = 'Classroom2E';
    image_names{6} = 'Computer';
    image_names{7} = 'Crusade';
    image_names{8} = 'CrusadeP';
    image_names{9} = 'Djembe';
    image_names{10} = 'DjembeL';
    image_names{11} = 'Hoops';
    image_names{12} = 'Livingroom';
    image_names{13} = 'Newkuba';
    image_names{14} = 'Plants';
    image_names{15} = 'Staircase';
    ndisp = [290, 290, 250, 610, 610, 256, 800, 800, 320, 320, 410, 320, 570, 320, 450];
end

% for im_num = 1:15
im_L = imread(['Stereo_data_example/MiddEval3/',imgset,imgsize,'/',image_names{im_idx},'/im0.png']);
im_R = imread(['Stereo_data_example/MiddEval3/',imgset,imgsize,'/',image_names{im_idx},'/im1.png']);

GT =0; mask=0;
if strcmp(imgset,'training')
    GT = readpfm(['Stereo_data_example/MiddEval3_GT/training',imgsize,'/',image_names{im_idx},'/disp0GT.pfm']);
    mask = imread(['Stereo_data_example/MiddEval3_GT/training',imgsize,'/',image_names{im_idx},'/mask0nocc.png']);
    mask = mask == 255;
%     Error = abs(DisparityMap{1} - GT) > 1;
%     Error(~mask) = 0;
%     ErrorRate(im_num) = sum(Error(:))/sum(mask(:));
%     fprintf('%s = %f\n', image_names{im_num}, ErrorRate(im_num));
end
% ndisp_out = ndisp(im_idx);
% im_GT = imread(['MiddEval3/',imgset,imgsize,'/',image_names{im_idx},'/im0.png']);
% GT = readpfm(['MiddEval3_GT/training',imgsize,'/',image_names{im_idx},'/disp0GT.pfm']);

% end

