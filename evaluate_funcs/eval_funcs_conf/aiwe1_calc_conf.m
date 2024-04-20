% # =============================================================================
% # @INPROCEEDINGS{punnappurath2020modeling,
% # author={Abhijith Punnappurath and Abdullah Abuolaim and Mahmoud Afifi and Michael S. Brown},
% # booktitle={IEEE International Conference on Computational Photography (ICCP)}, 
% # title={Modeling Defocus-Disparity in Dual-Pixel Sensors}, 
% # year={2020}}
% #
% # Author: Abhijith Punnappurath (05/2020)
% # pabhijith@eecs.yorku.ca
% # jithuthatswho@gmail.com
% # https://abhijithpunnappurath.github.io
% # =============================================================================

function aiwe1 = aiwe1_calc_conf(pred,GT,conf)
%% Based on googles method
max_iter = 5;
eps = 1e-3;

pred_vec = pred(:);
GT_vec = GT(:);
conf_vec = conf(:);
irls_weight = ones(size(conf_vec));

for idx_iter = 1:max_iter
    sqrt_weight = sqrt(irls_weight.*conf_vec);
    lhs = [sqrt_weight,sqrt_weight].* [pred_vec,ones(size(pred_vec))];
    rhs = sqrt_weight.* GT_vec;
    [affine_est,~] = lsqr(lhs,rhs);
    
    pred_affine = pred_vec.*affine_est(1) + affine_est(2);
    resid = abs(pred_affine - GT_vec);
    irls_weight = 1./(max(eps,resid(:)));
end
aiwe1 = sum(conf_vec.*resid)/sum(conf_vec);
    

% 
% 
% 
% 
% 
% 
% 
% A=[img1(:) ones(numel(img1),1)];
% b=img2(:);
% c=conf(:);
% W=spdiags(ones(numel(img1),1),0,numel(img1),numel(img1));
% x=(A'*W*A)\(A'*W*b);
% 
% for i=1:5
%       e=abs(A*x-b);
%       e=max(e,0.001).^(1-2);                    
%       W=spdiags(e,0,numel(img1),numel(img1));
%       x=(A'*W*A)\(A'*W*b);
% end
% 
% img3=x(1)*img1+x(2);
% aiwe1=mae_calc(img3,img2);