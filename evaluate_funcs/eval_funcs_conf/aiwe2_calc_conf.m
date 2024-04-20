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

function [aiwe2] = aiwe2_calc_conf(pred,GT,conf)
eps = 1e-3;

pred_vec = pred(:);
GT_vec = GT(:);
conf_vec = conf(:);

lhs = sqrt([conf_vec,conf_vec]).*[pred_vec,ones(size(pred_vec))];
rhs = sqrt(conf_vec).*GT_vec;
[affine_est,~] = lsqr(lhs,rhs);
pred_affine = pred_vec.*affine_est(1) + affine_est(2);

resid_sq = min((pred_affine-GT_vec).^2,1e20);

aiwe2 = sqrt(sum(conf_vec.*resid_sq)/sum(conf_vec));


% A=[img1(:) ones(numel(img1),1)];
% b=img2(:);
% x=A\b;
% img3=x(1)*img1+x(2);
% aiwe2=rmse_calc(img3,img2);
% a=x(1);
% b=x(2);