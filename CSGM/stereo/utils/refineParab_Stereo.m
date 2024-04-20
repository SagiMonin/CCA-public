function parab = refineParab_Stereo(parab, params)
% Create some refiniment of paraboals.
% First - check if they have low confidence (shallow curviature), if so than
% "flatten" them.
% Second - we are not as sure in the borders, so "flatten" parabolas near
% borders (note this is border in X-direction)

% Some parameters/hyper-params
% params.confidenceThresh = -max(abs(parab.a(:)))/100;
def_val_a = params.confidenceThresh*(1e-4);

% First step - remove low confidence curves

idx_filter = parab.a > params.confidenceThresh;
if 0
    figure(222); imagesc(abs(idx_filter));
end
% parab.a(idx_filter) = def_val_a;
% parab.b(idx_filter) = 0;
% parab.c(idx_filter) = 0;
% 
parab.a(idx_filter) = parab.a(idx_filter).*abs(def_val_a);
parab.b(idx_filter) = parab.b(idx_filter).*abs(def_val_a);
parab.c(idx_filter) = 0;
% % 


%% Give very low min for places 
x = -parab.b./(2*parab.a);
idx_min = (parab.a>-1e-32);
parab.a(idx_min) = -1e-32;
parab.b(idx_min) = x(idx_min).*(-2.*parab.a(idx_min)); % b - is adjusted so x is the minimum x = b/-2a -> b = -2ax
parab.b(isnan(parab.b)) = 0;

% Second step - penalize borders
border_penalty = ones(size(parab.a));
penalty_val = linspace(params.penalty_border, 1, params.border_len);
border_penalty(:,1:params.border_len) = repmat(penalty_val,[size(parab.a,1),1]);
border_penalty(:,end-params.border_len+1:end) = fliplr(repmat(penalty_val,[size(parab.a,1),1]));

parab.a=parab.a./border_penalty;
parab.b=parab.b./border_penalty;
parab.c=parab.c./border_penalty;

end


