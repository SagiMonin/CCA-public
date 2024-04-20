function sum_parab = propCSGM_test_stereo(parab, im_guide, params)

if size(im_guide,3)==3
    rot_im_guide=permute(im_guide,[2 1 3]);
else
    rot_im_guide=im_guide';
end

% Save the initial "a" values as they represent confidence and the
% propagation "ruins" this
aInit=parab.a./mean(parab.a(:)); % Sagi - try to use 1, and see what happens with larger iterations
% aInit = ones(size(parab.a));
% aInit = parab.a/max(abs(parab.a(:)));
% norma = max(abs(parab.a(:)));
% parab.a = parab.a/norma;
% parab.b = parab.b/norma;
% parab.c = parab.c/norma;

if params.numIter(params.idx_pyr) > 0
    sum_parab.a = zeros(size(parab.a));
    sum_parab.b = zeros(size(parab.a));
    sum_parab.c = zeros(size(parab.a));
    init_parab.a = parab.a;
    init_parab.b = parab.b;
    init_parab.c = parab.c;
    
elseif params.numIter(params.idx_pyr) == 0
    sum_parab.a = parab.a;
    sum_parab.b = parab.b;
    sum_parab.c = parab.c;
    init_parab.a = parab.a;
    init_parab.b = parab.b;
    init_parab.c = parab.c;
end

for iter=1:params.numIter(params.idx_pyr)
    params.thr_dis_curr = params.thr_dis(iter);
    params.thr_dis_curr2 = params.thr_dis2(iter);

    % Horizontal
    parab_h = propLine_test_stereo_P2(im_guide, parab, params, 'H');
%         parab_h = propLine_test_stereo(im_guide, parab, params, 'H');

    % Vertical
    parab_v = propLine_test_stereo_P2(rot_im_guide, parab, params ,'V');    
%         parab_v = propLine_test_stereo(rot_im_guide, parab, params ,'V');    

    % First pair of diagonals
    parab_MD = propDiag_test_stereo_P2(im_guide, parab, params, 'MD');    
%         parab_MD = propDiag_test_stereo(im_guide, parab, params, 'MD');    

    % Second pair of diagonals
    parab_OD = propDiag_test_stereo_P2(im_guide, parab, params, 'OD');
%         parab_OD = propDiag_test_stereo(im_guide, parab, params, 'OD');

    
    
    % Sum all parabolas
    sum_parab.a = parab_h.a + parab_v.a + parab_MD.a + parab_OD.a;
    sum_parab.b = parab_h.b + parab_v.b + parab_MD.b + parab_OD.b;
    sum_parab.c = parab_h.c + parab_v.c + parab_MD.c + parab_OD.c;
%     parab_cat.a = cat(3,parab_h.a, parab_v.a ,parab_MD.a , parab_OD.a);
%     parab_cat.b = cat(3,parab_h.b, parab_v.b ,parab_MD.b , parab_OD.b);
%     parab_cat.c = cat(3,parab_h.c, parab_v.c ,parab_MD.c , parab_OD.c);
%     [val,idx ] = min(parab_cat.a,[],3);
%     sum_parab.a = zeros(size(sum_parab.a));
%     sum_parab.c = zeros(size(sum_parab.c));
%     sum_parab.b = zeros(size(sum_parab.b));
%     for idx_iter = 1:4
%         sum_parab.a = sum_parab.a + parab_cat.a(:,:,idx_iter)/2 .* double(idx==idx_iter);
%         sum_parab.b = sum_parab.b + parab_cat.b(:,:,idx_iter)/2 .* double(idx==idx_iter);
%         sum_parab.c = sum_parab.c + parab_cat.c(:,:,idx_iter)/2 .* double(idx==idx_iter);
%     end
%     sum_parab = refineParab2(sum_parab, params);

    
    % Prepare the next iteration
    if 1 %iter~=params.numIter(params.idx_pyr)
        numPaths=8;
        parab.a = sum_parab.a./numPaths.*aInit;
        parab.b = sum_parab.b./numPaths.*aInit;
        parab.c = sum_parab.c./numPaths.*aInit;
        parab = refineParab2(parab, params);
        
    
    end
   
        
    if 0
        GT = evalin('base', 'GT');

        figure(55);
        subplot(231); imagesc(init_parab.b./(2*init_parab.a));
        subplot(232); imagesc(parab.b./(2*parab.a));
        subplot(233); imagesc(abs(-parab.b./(2*parab.a)+GT)<1); title(num2str(sum(sum(abs(-parab.b./(2*parab.a)+GT)<1))/numel(GT(:))*100));
        subplot(234); imagesc(init_parab.a); 
        subplot(235); imagesc(parab.a); 

%             figure(56);
%             ax(iter+params.numIter(params.idx_pyr)) =subplot(3,4,iter);  imagesc(parab.a);
    end
    
end

end


