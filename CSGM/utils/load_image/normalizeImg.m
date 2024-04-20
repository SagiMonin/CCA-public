function [nrmImg,stdVal]=normalizeImg(img,kerSize)

k=ones(kerSize)./sum(sum(ones(kerSize)));
meanVal=imfilter(img,k,'symmetric');
img=img-meanVal;
meanValSqr=imfilter(img.*img,k,'symmetric');
stdVal=sqrt(meanValSqr);
nrmImg=img./(max(stdVal,1));
nrmImg=imguidedfilter(nrmImg,img,'NeighborhoodSize',[5 5],'DegreeOfSmoothing' ,0.001*256*256); 
end