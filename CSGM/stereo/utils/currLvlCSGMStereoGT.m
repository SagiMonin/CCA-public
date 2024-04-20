function [dispar_map, sum_parab] = currLvlCSGMStereoGT(im_L,im_R,im_guide, params, prior_lvl_parab, cost_neig, conf_score, dispar_int_val,GT )
%% Calculate current level C-SGM

% Mark which level we are
flag_lvl = 1; % Current level in upscale
if ~isempty(prior_lvl_parab)
    flag_lvl = 2; % This is not the first lvl
end

% Sub-pixel estimation
if ~strcmp(params.interpolant, 'ENCC') % Not ENCC
    % Generate parab after sub-pixel estimation
    parab = genParab(cost_neig, dispar_int_val, params);
elseif strcmp(params.interpolant, 'ENCC') % ENCC
    % Calculate cost and neighboors + confidence score + disparity integer value.
    [parab,invalid_mask] = estSubPixENCC(cost_neig, dispar_int_val, im_L, im_R, params);
    conf_score(invalid_mask == 1) = 0.01^2; % Invalid conf
end
%
% parab.a =(left_val+right_val)./2; % a - calculated with the 3 original cost values (ignoring sub-pixel refinment), so just take second derivative 
parab.b = GT.*(-2.*parab.a); % b - is adjusted so x is the minimum x = b/-2a -> b = -2ax
% parab.c = mid_val; % Don't really need this value
%


% Multiply with confidence scores
parab.a = parab.a.*conf_score;
parab.b = parab.b.*conf_score;
parab.c = parab.c.*conf_score;    

% If we have data from prior from previous level of pyramid sacle the
% parabolas
if flag_lvl > 1
    parab.a = parab.a + params.priorW * prior_lvl_parab.a;
    parab.b = parab.b + params.priorW * prior_lvl_parab.b;
    parab.c = parab.c + params.priorW * prior_lvl_parab.c;
end

% Refine parabaols - low confidence + at borders
parab = refineParab(parab, params);

% Cost prop 8-way
sum_parab = propCSGM(parab, im_guide, params);

% final disparity
dispar_map =-sum_parab.b./(2.*sum_parab.a);


end