function parab_prop = propLine_test_stereo_P2(im_guide, parab, params, dir_flag)
% Extract values 
P1param = params.P1param;
P2param = params.P2param;
sigmaEdges = params.sigmaEdges;
Thrrat = params.thr_rat;
Thr1 = params.thr_dis_curr ;
Thr2 = params.thr_dis_curr2 ;
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
    P2=ComputePAdaptive(Pedges(:,ii),orig_a(:,ii),prevA,P2param);

%      Ep=ComputeExpectedValue(prevA,prevB);
%     [SaLR(:,ii),SbLR(:,ii),ScLR(:,ii)]=ComputePropagatedParabola(SaLR(:,ii),SbLR(:,ii),ScLR(:,ii),Padaptive,Ep);

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
    idx_edge = Pedges(:,ii)>ThrEdge; % Find weak edge - so it is okay to prop
    
    idx_p3 = idx_p3.*idx_edge.*idx_large_rat;
    
    % Rest of index - don't change.
    rest_idx = ~(idx_p1|idx_p2|idx_p3);
    
    % Propagation - regular - both - Nothing, large disparity + not edge, regular
    SaLR(:,ii) = (currA - P1).*double(idx_p1) + (currA - P2).*double(idx_p2) + prevA.*double(idx_p3) + currA.*double(rest_idx);
    SbLR(:,ii) = (currB +  2.*P1.*Ep_prev).*double(idx_p1) + (currB +  2.*P2.*Ep_prev).*double(idx_p2)+ prevB.*double(idx_p3)+ currB.*double(rest_idx);
    ScLR(:,ii) = (currC - P1.*Ep_prev.^2).*double(idx_p1) + (currC - P2.*Ep_prev.^2).*double(idx_p2)+ prevC.*double(idx_p3)+ currC.*double(rest_idx);

    
    % Propagation - regular - both - Nothing, large disparity + not edge, regular
%     SaLR(:,ii) = currA.*double(idx_p2.*(~large_disp)) + prevA.*double(idx_p2.*large_disp) + (currA - P1).*double(~idx_p2);
%     SbLR(:,ii) = currB.*double(idx_p2.*(~large_disp)) + prevB.*double(idx_p2.*large_disp) + (currB +  2.*P1.*Ep_prev).*double(~idx_p2);
%     ScLR(:,ii) = currC.*double(idx_p2.*(~large_disp)) + prevC.*double(idx_p2.*large_disp) + (currC - P1.*Ep_prev.^2).*double(~idx_p2) ;
% % % 
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

%     P1=-P1param.*Pedges(:,ii).*prevA;%*max(weightPrev-currentToPrevWeight.*weightCur,0);
    P1=ComputePAdaptive(Pedges(:,ii),orig_a(:,ii),prevA,P1param);
    P2=ComputePAdaptive(Pedges(:,ii),orig_a(:,ii),prevA,P2param);

%      Ep=ComputeExpectedValue(prevA,prevB);
%     [SaLR(:,ii),SbLR(:,ii),ScLR(:,ii)]=ComputePropagatedParabola(SaLR(:,ii),SbLR(:,ii),ScLR(:,ii),Padaptive,Ep);

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
    idx_edge = Pedges(:,ii)>ThrEdge; % Find weak edge - so it is okay to prop
    
    idx_p3 = idx_p3.*idx_edge.*idx_large_rat;
    
    % Rest of index - don't change.
    rest_idx = ~(idx_p1|idx_p2|idx_p3);
    
    % Propagation - regular - both - Nothing, large disparity + not edge, regular
    SaRL(:,ii) = (currA - P1).*double(idx_p1) + (currA - P2).*double(idx_p2) + prevA.*double(idx_p3) + currA.*double(rest_idx);
    SbRL(:,ii) = (currB +  2.*P1.*Ep_prev).*double(idx_p1) + (currB +  2.*P2.*Ep_prev).*double(idx_p2)+ prevB.*double(idx_p3)+ currB.*double(rest_idx);
    ScRL(:,ii) = (currC - P1.*Ep_prev.^2).*double(idx_p1) + (currC - P2.*Ep_prev.^2).*double(idx_p2)+ prevC.*double(idx_p3)+ currC.*double(rest_idx);
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

