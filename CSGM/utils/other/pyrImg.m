function [blurred,diffImg,imgs]=pyrImg(img, num_lvl)


max_num_lvl = floor(log2(min(size(img,1),size(img,2)))); % Calculate maximum possible levels of pyramid 
max_num_lvl = min(max_num_lvl, num_lvl);

imgs{max_num_lvl} = img;
% [blurred{ii},diffImg{ii},imgs{ii-1}] = doStep(imgs{ii});
for ii = max_num_lvl:-1:1
    [blurred{ii},diffImg{ii}] = doStep(imgs{ii});
    if ii > 1
        imgs{ii - 1} = impyramid(imgs{ii}, 'reduce');
    end
end


% Single pyramid reduction
function [blurred,diffImg]=doStep(img)
    sigma=1;
    
    blurred=imgaussfilt(img,sigma,'FilterSize' ,[5 5]); % This is not the filter used in the pyramid algo in matlab
    diffImg=img-blurred;

%     smallImg = impyramid(img, 'reduce');
    %diffImg=imfilter(blurred,fspecial('laplacian',0.2));

    %smallImg=imresize(blurred,size(img)/2);
end

end

