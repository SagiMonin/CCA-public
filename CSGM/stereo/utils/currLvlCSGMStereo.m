function [dispar_map, sum_parab] = currLvlCSGMStereo(im_L,im_R,im_guide, params, prior_lvl_parab, cost_neig, conf_score, dispar_int_val )
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

% Multiply with confidence scores
parab.a = parab.a.*conf_score;
parab.b = parab.b.*conf_score;
parab.c = parab.c.*conf_score;    

if 0
    figure(99);
    subplot(121); imagesc(parab.b./(2.*parab.a));
end

% If we have data from prior from previous level of pyramid sacle the
% parabolas
if flag_lvl > 1
    parab.a = parab.a + params.priorW * prior_lvl_parab.a;
    parab.b = parab.b + params.priorW * prior_lvl_parab.b;
    parab.c = parab.c + params.priorW * prior_lvl_parab.c;
end

% Refine parabaols - low confidence + at borders
parab = refineParab_Stereo(parab, params);

if 0
    figure(99);
    subplot(121); imagesc(parab.a);

    subplot(122); imagesc(parab.b./(2.*parab.a));
end

% Cost prop 8-way
sum_parab = propCSGM_test_stereo(parab, im_guide, params);

% final disparity
dispar_map =-sum_parab.b./(2.*sum_parab.a);


end