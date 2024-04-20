function [good05,good1,good2,good4,rmse] = metrics_stereo(GT,res,mask,shift_val	,idx_fig)

if sum(sum(mask - ones(size(mask)))) == 0
    mask(GT==0) = 0;
end

GT = GT.*mask; 
res = res.*mask; 

if shift_val>0
    GT = GT(:,shift_val:end);
    res = res(:,shift_val:end);
end
num_pix = numel(res);



good05 = 100-sum(abs(GT(:) - res(:))<0.5)/num_pix*100;
good1 = 100-sum(abs(GT(:) - res(:))<1)/num_pix*100;
good2 = 100-sum(abs(GT(:) - res(:))<2)/num_pix*100;
good4 = 100-sum(abs(GT(:) - res(:))<4)/num_pix*100;
rmse = sqrt(sum((GT(:)-res(:)).^2)/num_pix);
if 1
    if idx_fig>0
        figure(idx_fig);
        subplot(141); imagesc(res.*(abs(GT - res)<0.5)); title(['0.5 pix bad- ', num2str(good05),'%']);
        set(gca,'FontSize',18) % Creates an axes and sets its FontSize to 18

        subplot(142); imagesc(res.*(abs(GT - res)<1)); title(['1 pix bad- ', num2str(good1),'%']);
        set(gca,'FontSize',18) % Creates an axes and sets its FontSize to 18

        subplot(143); imagesc(res.*(abs(GT - res)<2)); title(['2 pix bad- ', num2str(good2),'%']);
        set(gca,'FontSize',18) % Creates an axes and sets its FontSize to 18

        subplot(144); imagesc(res.*(abs(GT - res)<4)); title(['4 pix bad- ', num2str(good4),'%']);
        set(gca,'FontSize',18) % Creates an axes and sets its FontSize to 18
    end
else 
    figure;
    imagesc((abs(GT - res)<0.5)); title(['0.5 pix bad- ', num2str(good05),'%']); set(gca,'FontSize',18) % Creates an axes and sets its FontSize to 18
%     subplot(122);imagesc(res); title('disparity'); set(gca,'FontSize',18) % Creates an axes and sets its FontSize to 18
end

end

