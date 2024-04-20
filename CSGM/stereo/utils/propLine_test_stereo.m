function parab_prop = propLine_test_stereo(im_guide, parab, params, dir_flag)
% Extract values 
P1param = params.P1param;
sigmaEdges = params.sigmaEdges;
Thrrat = params.thr_rat;
ThrDis = params.thr_dis_curr ;
ThrEdge = params.thr_edge ;
% currentToPrevWeight = params.currentToPrevWeight;

% left to right propagation
df = sum(im_guide(:,2:end,:)-im_guide(:,1:end-1,:),3)./size(im_guide,3); % Calculate difference (forward operation)
Pedges = [ones(size(im_guide(:,1))) exp(-df.^2./sigmaEdges^2)]; % Scale prior value for edges.

% BW = edge(rgb2gray(im_guide/255),'canny',[params.thr_canny]);
% BW = imdilate(BW,[1 1 1; 1 1 1; 1 1 1]);
% Pedges = Pedges.*abs(0.99999-BW);


if dir_flag == 'H'
    SaLR = parab.a;
    SbLR = parab.b;
    ScLR = parab.c;
    orig_a = parab.a;
elseif dir_flag == 'V'
    SaLR = parab.a';
    SbLR = parab.b';
    ScLR = parab.c';
    orig_a = parab.a';
end

for ii=2:size(SaLR,2)

    prevA=SaLR(:,ii-1);
    prevB=SbLR(:,ii-1);
    prevC=ScLR(:,ii-1);
    currA=SaLR(:,ii);
    currB=SbLR(:,ii);
    currC=ScLR(:,ii);

%     P1=-P1param.*Pedges(:,ii).*prevA;%*max(weightPrev-currentToPrevWeight.*weightCur,0);
    P1=ComputePAdaptive(Pedges(:,ii),orig_a(:,ii),prevA,P1param);
%      Ep=ComputeExpectedValue(prevA,prevB);
%     [SaLR(:,ii),SbLR(:,ii),ScLR(:,ii)]=ComputePropagatedParabola(SaLR(:,ii),SbLR(:,ii),ScLR(:,ii),Padaptive,Ep);

    Ep_prev=ComputeExpectedValue(prevA,prevB);
    Ep_curr=ComputeExpectedValue(currA,currB);
    diff_val = abs(Ep_prev-Ep_curr);
    idx_p2 = diff_val>ThrDis; % Dipsarity large than thr
    rat_p2 = prevA./currA;
    rat_large = rat_p2>Thrrat; % Ratio of disparity is large enough
    idx_edge = Pedges(:,ii)>ThrEdge; % weak edge
    large_disp = rat_large&idx_edge;
% 
    
%     % Propagation - regular
%     SaLR(:,ii) = currA.*double(idx_p2) + (currA - P1).*double(~idx_p2);
%     SbLR(:,ii) = currB.*double(idx_p2) + (currB +  2.*P1.*Ep_prev).*double(~idx_p2);
%     ScLR(:,ii) = currC.*double(idx_p2) + (currC - P1.*Ep_prev.^2).*double(~idx_p2) ;
%     
%     % Propagation -  large disparity - decide which disparity to take
%     SaLR(:,ii) = currA.*double(idx_p2.*(~rat_large)) + prevA.*double(idx_p2.*rat_large);
%     SbLR(:,ii) = currB.*double(idx_p2.*(~rat_large)) + prevB.*double(idx_p2.*rat_large);
%     ScLR(:,ii) = currC.*double(idx_p2.*(~rat_large)) + prevC.*double(idx_p2.*rat_large);
    
    
    % Propagation - regular - both - Nothing, large disparity + not edge, regular
    SaLR(:,ii) = currA.*double(idx_p2.*(~large_disp)) + prevA.*double(idx_p2.*large_disp) + (currA - P1).*double(~idx_p2);
    SbLR(:,ii) = currB.*double(idx_p2.*(~large_disp)) + prevB.*double(idx_p2.*large_disp) + (currB +  2.*P1.*Ep_prev).*double(~idx_p2);
    ScLR(:,ii) = currC.*double(idx_p2.*(~large_disp)) + prevC.*double(idx_p2.*large_disp) + (currC - P1.*Ep_prev.^2).*double(~idx_p2) ;
% % 
%     figure(250); 
%     subplot(131); imagesc(SbLR./(2*SaLR))
%     subplot(132); imagesc(SaLR)
%     subplot(133); imagesc(Pedges>ThrEdge)

end

% right to left propagation
if dir_flag == 'H'
    SaRL = parab.a;
    SbRL = parab.b;
    ScRL = parab.c;
    orig_a = parab.a;
elseif dir_flag == 'V'
    SaRL = parab.a';
    SbRL = parab.b';
    ScRL = parab.c';
    orig_a = parab.a';
end

df = sum(im_guide(:,1:end-1,:)-im_guide(:,2:end,:),3)./size(im_guide,3);
Pedges=[ exp(-df.^2./sigmaEdges.^2) ones(size(im_guide(:,1)))];
% Pedges = Pedges.*abs(0.999-BW);

for ii=size(SaRL,2)-1:-1:1

%     prevA=SaRL(:,ii+1);
%     prevB=SbRL(:,ii+1);
%     
%     Padaptive=ComputePAdaptive(Pedges(:,ii),orig_a(:,ii),prevA,P1param);
%     Ep=ComputeExpectedValue(prevA,prevB);
%     [SaRL(:,ii),SbRL(:,ii),ScRL(:,ii)]=ComputePropagatedParabola(SaRL(:,ii),SbRL(:,ii),ScRL(:,ii),Padaptive,Ep);
    
    prevA=SaRL(:,ii+1);
    prevB=SbRL(:,ii+1);
    prevC=ScRL(:,ii+1);
    currA=SaRL(:,ii);
    currB=SbRL(:,ii);
    currC=ScRL(:,ii);

%     P1=-P1param.*Pedges(:,ii).*prevA;
    P1=ComputePAdaptive(Pedges(:,ii),orig_a(:,ii),prevA,P1param);

    Ep_prev=ComputeExpectedValue(prevA,prevB);
    Ep_curr=ComputeExpectedValue(currA,currB);
    diff_val = abs(Ep_prev-Ep_curr);
    idx_p2 = diff_val>ThrDis; % Dipsarity large than 1
    rat_p2 = prevA./currA;
    rat_large = rat_p2>Thrrat; % Ratio of disparity is large enough
    idx_edge = Pedges(:,ii)<ThrEdge; % weak edge
    large_disp = rat_large&idx_edge;
    
%     % Propagation - regular
%     SaRL(:,ii) = currA.*double(idx_p2) + (currA - P1).*double(~idx_p2);
%     SbRL(:,ii) = currB.*double(idx_p2) + (currB +  2.*P1.*Ep_prev).*double(~idx_p2);
%     ScRL(:,ii) = currC.*double(idx_p2) + (currC - P1.*Ep_prev.^2).*double(~idx_p2) ;
%     
%     % Propagation -  large disparity - decide which disparity to take
%     SaRL(:,ii) = currA.*double(idx_p2.*(~rat_large)) + prevA.*double(idx_p2.*rat_large);
%     SbRL(:,ii) = currB.*double(idx_p2.*(~rat_large)) + prevB.*double(idx_p2.*rat_large);
%     ScRL(:,ii) = currC.*double(idx_p2.*(~rat_large)) + prevC.*double(idx_p2.*rat_large) ;    

    SaRL(:,ii) = currA.*double(idx_p2.*(~large_disp)) + prevA.*double(idx_p2.*large_disp) + (currA - P1).*double(~idx_p2);
    SbRL(:,ii) = currB.*double(idx_p2.*(~large_disp)) + prevB.*double(idx_p2.*large_disp) + (currB +  2.*P1.*Ep_prev).*double(~idx_p2);
    ScRL(:,ii) = currC.*double(idx_p2.*(~large_disp)) + prevC.*double(idx_p2.*large_disp) + (currC - P1.*Ep_prev.^2).*double(~idx_p2) ;

%     figure(250); 
%     subplot(121); imagesc(SbRL./(2*SaRL))
%     subplot(122); imagesc(SaRL)
    
end

if dir_flag == 'H'
    parab_prop.a = SaLR+SaRL;
    parab_prop.b = SbLR+SbRL;
    parab_prop.c = ScLR+ScRL;
elseif dir_flag == 'V'
    parab_prop.a = (SaLR+SaRL)';
    parab_prop.b = (SbLR+SbRL)';
    parab_prop.c = (ScLR+ScRL)';
end


end

