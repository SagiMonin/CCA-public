function parab_prop = propDiag_test_stereo(im_guide, parab, params, dir_flag)
% Extract values 
P1param = params.P1param;
sigmaEdges = params.sigmaEdges;
Thrrat = params.thr_rat;
ThrDis = params.thr_dis_curr;
ThrEdge = params.thr_edge ;

if dir_flag == 'OD' % off-diagonal
    im_guide = fliplr(im_guide);
end

df=sum(im_guide(2:end,2:end,:)-im_guide(1:end-1,1:end-1,:),3)./size(im_guide,3);
PedgesTopLeft=exp(-df.^2./sigmaEdges^2);
PedgesTopLeft=[ones([size(PedgesTopLeft,1) 1]) PedgesTopLeft];
PedgesTopLeft=[ones([1 size(PedgesTopLeft,2)]); PedgesTopLeft];

% BW = edge(rgb2gray(im_guide/255),'canny',[params.thr_canny]);
% BW = imdilate(BW,[1 1 1; 1 1 1; 1 1 1]);
% PedgesTopLeft = PedgesTopLeft.*abs(0.99999-BW);

if dir_flag == 'MD' % Main diagonal
    SaLR = parab.a;
    SbLR = parab.b;
    ScLR = parab.c;
    orig_a = parab.a;

elseif dir_flag == 'OD' % Off-diagional
    SaLR = fliplr(parab.a);
    SbLR = fliplr(parab.b);
    ScLR = fliplr(parab.c);
    orig_a = fliplr(parab.a);

end

for ii=2:size(SaLR,2)    
    % Top left
%     prevATopLeft=[ SaLR(1,ii); SaLR(1:end-1,ii-1)];
%     prevBTopLeft=[ SbLR(1,ii); SbLR(1:end-1,ii-1)];
%    
%     % Diagonal aggegation
%     Padaptive=ComputePAdaptive(PedgesTopLeft(:,ii),orig_a(:,ii),prevATopLeft,P1param);
%     Ep=ComputeExpectedValue(prevATopLeft,prevBTopLeft);
%     [SaLR(:,ii),SbLR(:,ii),ScLR(:,ii)]=ComputePropagatedParabola(SaLR(:,ii),SbLR(:,ii),ScLR(:,ii),Padaptive,Ep);
    
    prevA = [SaLR(1,ii); SaLR(1:end-1,ii-1)];
    prevB = [SbLR(1,ii); SbLR(1:end-1,ii-1)];
    prevC = [ScLR(1,ii); ScLR(1:end-1,ii-1)];
    currA = SaLR(:,ii);
    currB = SbLR(:,ii);
    currC = ScLR(:,ii);

%     P1=-P1param.*PedgesTopLeft(:,ii).*prevA;%*max(weightPrev-currentToPrevWeight.*weightCur,0);
    P1=ComputePAdaptive(PedgesTopLeft(:,ii),orig_a(:,ii),prevA,P1param);
    P2=ComputePAdaptive(PedgesTopLeft(:,ii),orig_a(:,ii),prevA,P2param);

    Ep_prev=ComputeExpectedValue(prevA,prevB);
    Ep_curr=ComputeExpectedValue(currA,currB);
    diff_val = abs(Ep_prev-Ep_curr);
    
    
        % Diff thr_1 - 
    idx_p1 = diff_val<Thr1; % Disparity is small so prop with p1 
    
    % Diff thr_2 - strong diffusion 
    idx_p2 = diff_val>=Thr1 & diff_val<Thr2; % Disparity is small so prop with p1 
    
    % Diff thr_3 
    idx_p3 = diff_val>=Thr2; % rest of pixels
    
    % Find if ratio between previous and current is large enough 
    rat_prev_curr = prevA./currA;       
    idx_large_rat = rat_prev_curr>Thrrat; % Ratio of 
    idx_edge = PedgesTopLeft(:,ii)>ThrEdge; % Find weak edge - so it is okay to prop
    
    idx_p3 = idx_p3.*idx_edge.*idx_large_rat;
    
    % Rest of index - don't change.
    rest_idx = ~(idx_p1|idx_p2|idx_p3);
    
    % Propagation - regular - both - Nothing, large disparity + not edge, regular
    SaLR(:,ii) = (currA - P1).*double(idx_p1) + (currA - P2).*double(idx_p2) + prevA.*double(idx_p3) + currA.*double(rest_idx);
    SbLR(:,ii) = (currB +  2.*P1.*Ep_prev).*double(idx_p1) + (currB +  2.*P2.*Ep_prev).*double(idx_p2)+ prevB.*double(idx_p3)+ currB.*double(rest_idx);
    ScLR(:,ii) = (currC - P1.*Ep_prev.^2).*double(idx_p1) + (currC - P2.*Ep_prev.^2).*double(idx_p2)+ prevC.*double(idx_p3)+ currC.*double(rest_idx);
    
%     idx_p2 = diff_val>ThrDis; % Dipsarity large than 1
%     rat_p2 = prevA./currA;
%     rat_large = rat_p2>Thrrat; % Ratio of disparity is large enough
%     idx_edge = PedgesTopLeft(:,ii)>ThrEdge; % weak edge
%     large_disp = rat_large&idx_edge;
%     
% %     % Propagation - regular
% %     SaLR(:,ii) = currA.*double(idx_p2) + (currA - P1).*double(~idx_p2);
% %     SbLR(:,ii) = currB.*double(idx_p2) + (currB +  2.*P1.*Ep_prev).*double(~idx_p2);
% %     ScLR(:,ii) = currC.*double(idx_p2) + (currC - P1.*Ep_prev.^2).*double(~idx_p2) ;
%     
% %     % Propagation -  large disparity - decide which disparity to take
% %     SaLR(:,ii) = currA.*double(idx_p2.*(~rat_large)) + prevA.*double(idx_p2.*rat_large);
% %     SbLR(:,ii) = currB.*double(idx_p2.*(~rat_large)) + prevB.*double(idx_p2.*rat_large);
% %     ScLR(:,ii) = currC.*double(idx_p2.*(~rat_large)) + prevC.*double(idx_p2.*rat_large);
%     
%     SaLR(:,ii) = currA.*double(idx_p2.*(~large_disp)) + prevA.*double(idx_p2.*large_disp) + (currA - P1).*double(~idx_p2);
%     SbLR(:,ii) = currB.*double(idx_p2.*(~large_disp)) + prevB.*double(idx_p2.*large_disp) + (currB +  2.*P1.*Ep_prev).*double(~idx_p2);
%     ScLR(:,ii) = currC.*double(idx_p2.*(~large_disp)) + prevC.*double(idx_p2.*large_disp) + (currC - P1.*Ep_prev.^2).*double(~idx_p2) ;
%     
%     figure(250); 
%     subplot(121); imagesc(SbLR./(2*SaLR))
%     subplot(122); imagesc(SaLR)
    
end

if dir_flag == 'MD' % Main diagonal
    SaRL = parab.a;
    SbRL = parab.b;
    ScRL = parab.c;
    orig_a = parab.a;
elseif dir_flag == 'OD' % Off-diagional
    SaRL = fliplr(parab.a);
    SbRL = fliplr(parab.b);
    ScRL = fliplr(parab.c);
    orig_a = fliplr(parab.a);
end

df=sum(im_guide(2:end,2:end,:)-im_guide(1:end-1,1:end-1,:),3)./size(im_guide,3);
PedgesBotLeft=exp(-df.^2./sigmaEdges^2);
PedgesBotLeft=[ PedgesBotLeft ones([size(PedgesBotLeft,1) 1])];
PedgesBotLeft=[ PedgesBotLeft;ones([1 size(PedgesBotLeft,2)])];

% PedgesBotLeft = PedgesBotLeft.*abs(0.999-BW);


for ii=2:size(SaLR,2)
    % Bot left
%     prevABotLeft=[ SaRL(2:end,ii-1);SaRL(end,ii) ];
%     prevBBotLeft=[ SbRL(2:end,ii-1);SbRL(end,ii) ];
%     
%     Padaptive=ComputePAdaptive(PedgesBotLeft(:,ii),orig_a(:,ii),prevABotLeft,P1param);
%     Ep=ComputeExpectedValue(prevABotLeft,prevBBotLeft);
%     [SaRL(:,ii),SbRL(:,ii),ScRL(:,ii)]=ComputePropagatedParabola(SaRL(:,ii),SbRL(:,ii),ScRL(:,ii),Padaptive,Ep);
    
    
    prevA = [ SaRL(2:end,ii-1);SaRL(end,ii) ];
    prevB = [ SbRL(2:end,ii-1);SbRL(end,ii) ];
    prevC = [ ScRL(2:end,ii-1);ScRL(end,ii) ];
    currA = SaRL(:,ii);
    currB = SbRL(:,ii);
    currC = ScRL(:,ii);

%     P1=-P1param.*PedgesBotLeft(:,ii).*prevA;%*max(weightPrev-currentToPrevWeight.*weightCur,0);
    P1=ComputePAdaptive(PedgesBotLeft(:,ii),orig_a(:,ii),prevA,P1param);
    P2=ComputePAdaptive(PedgesBotLeft(:,ii),orig_a(:,ii),prevA,P2param);

    Ep_prev=ComputeExpectedValue(prevA,prevB);
    Ep_curr=ComputeExpectedValue(currA,currB);
    diff_val = abs(Ep_prev-Ep_curr);
    
    
    % Diff thr_1 - 
    idx_p1 = diff_val<Thr1; % Disparity is small so prop with p1 
    
    % Diff thr_2 - strong diffusion 
    idx_p2 = diff_val>=Thr1 & diff_val<Thr2; % Disparity is small so prop with p1 
    
    % Diff thr_3 
    idx_p3 = diff_val>=Thr2; % rest of pixels
    
    % Find if ratio between previous and current is large enough 
    rat_prev_curr = prevA./currA;       
    idx_large_rat = rat_prev_curr>Thrrat; % Ratio of 
    idx_edge = PedgesBotLeft(:,ii)>ThrEdge; % Find weak edge - so it is okay to prop
    
    idx_p3 = idx_p3.*idx_edge.*idx_large_rat;
    
    % Rest of index - don't change.
    rest_idx = ~(idx_p1|idx_p2|idx_p3);
    
    % Propagation - regular - both - Nothing, large disparity + not edge, regular
    SaRL(:,ii) = (currA - P1).*double(idx_p1) + (currA - P2).*double(idx_p2) + prevA.*double(idx_p3) + currA.*double(rest_idx);
    SbRL(:,ii) = (currB +  2.*P1.*Ep_prev).*double(idx_p1) + (currB +  2.*P2.*Ep_prev).*double(idx_p2)+ prevB.*double(idx_p3)+ currB.*double(rest_idx);
    ScRL(:,ii) = (currC - P1.*Ep_prev.^2).*double(idx_p1) + (currC - P2.*Ep_prev.^2).*double(idx_p2)+ prevC.*double(idx_p3)+ currC.*double(rest_idx);
    
    
    
    
    
%     idx_p2 = diff_val>ThrDis; % Dipsarity large than 1
%     rat_p2 = prevA./currA;
%     rat_large = rat_p2>Thrrat; % Ratio of disparity is large enough
%     idx_edge = PedgesBotLeft(:,ii)>ThrEdge; % weak edge
%     large_disp = rat_large&idx_edge;
    
    
%     % Propagation - regular
%     SaRL(:,ii) = currA.*double(idx_p2) + (currA - P1).*double(~idx_p2);
%     SbRL(:,ii) = currB.*double(idx_p2) + (currB +  2.*P1.*Ep_prev).*double(~idx_p2);
%     ScRL(:,ii) = currC.*double(idx_p2) + (currC - P1.*Ep_prev.^2).*double(~idx_p2) ;
%     
%     % Propagation -  large disparity - decide which disparity to take
%     SaRL(:,ii) = currA.*double(idx_p2.*(~rat_large)) + prevA.*double(idx_p2.*rat_large);
%     SbRL(:,ii) = currB.*double(idx_p2.*(~rat_large)) + prevB.*double(idx_p2.*rat_large);
%     ScRL(:,ii) = currC.*double(idx_p2.*(~rat_large)) + prevC.*double(idx_p2.*rat_large);
    
%     SaRL(:,ii) = currA.*double(idx_p2.*(~large_disp)) + prevA.*double(idx_p2.*large_disp) + (currA - P1).*double(~idx_p2);
%     SbRL(:,ii) = currB.*double(idx_p2.*(~large_disp)) + prevB.*double(idx_p2.*large_disp) + (currB +  2.*P1.*Ep_prev).*double(~idx_p2);
%     ScRL(:,ii) = currC.*double(idx_p2.*(~large_disp)) + prevC.*double(idx_p2.*large_disp) + (currC - P1.*Ep_prev.^2).*double(~idx_p2) ;
    
    
end

if dir_flag == 'MD' % Main diagonal
    parab_prop.a = SaLR+SaRL;
    parab_prop.b = SbLR+SbRL;
    parab_prop.c = ScLR+ScRL;
elseif dir_flag == 'OD' % Off-diagional
    parab_prop.a = fliplr(SaLR+SaRL);
    parab_prop.b = fliplr(SbLR+SbRL);
    parab_prop.c = fliplr(ScLR+ScRL);
end



end