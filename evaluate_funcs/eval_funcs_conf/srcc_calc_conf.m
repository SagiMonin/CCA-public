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

function srcc = srcc_calc_conf(d1,d2,conf)
d1_vec = d1(:);
d2_vec = d2(:);
conf_vec = conf(:);

% Rank
d1_rank = cast_rescale(rank(d1_vec));
d1_rank_neg = cast_rescale(rank(-d1_vec));
d2_rank = cast_rescale(rank(d2_vec));

p1 = person_corr(d1_rank,d2_rank,conf_vec);
p2 = person_corr(d1_rank_neg,d2_rank,conf_vec);

srcc = 1- max(p1,p2);


end

function y = cast_rescale(x)
y = (x-numel(x)/2)./(numel(x)/2);
end

function y = rank(x)
[~,I] = sort(x);
[~,y] = sort(I);
y = y-1;
end

function y = person_corr(d1_rank,d2_rank,conf_vec)
conf_sum = sum(conf_vec);
mu_d1 = expectation(d1_rank,conf_vec,conf_sum);
mu_d2 = expectation(d2_rank,conf_vec,conf_sum);
var_d1 = expectation(d1_rank.^2,conf_vec,conf_sum) - mu_d1.^2;
var_d2 = expectation(d2_rank.^2,conf_vec,conf_sum) - mu_d2.^2;
cov = expectation(d1_rank.*d2_rank,conf_vec,conf_sum) - mu_d1.*mu_d2;
y = cov./sqrt(var_d1.*var_d2);
end


function mu_x = expectation(x,conf,conf_sum)
mu_x = sum(x.*conf)/conf_sum;
end