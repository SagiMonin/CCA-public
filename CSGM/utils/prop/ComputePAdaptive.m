function Padaptive=ComputePAdaptive(Pedges,a,prevA,P1param)

% The weight of the previous parabola depends on its "a" coefficient.
% We don't want to weight in bad previous parabolas

% 
weightPrev=-prevA;
% if sum(weightPrev>1)
%     display('hi');
% end
weightCur=-a;
weightPrev=min(max(weightPrev,0),1); % Give more weight to high level of confidence depednign to previous iter
weightCur=min(max(weightCur,0),1);


% % if the quality of parabolas are equal, then we still want to weight in
% % the previous parabola for smoothness, this is why we lower the weight of
% % the weightCur
% Padaptive=P1param.*Pedges.*weightCur;%*max(weightPrev-currentToPrevWeight.*weightCur,0);

Pedges = single(Pedges);
Padaptive=P1param.*Pedges.*weightPrev;%*max(weightPrev-currentToPrevWeight.*weightCur,0);