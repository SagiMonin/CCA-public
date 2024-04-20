function parab_prop = propLine(im_guide, parab, params, dir_flag)
% Extract values 
P1param = params.P1param;
sigmaEdges = params.sigmaEdges;
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
    
    Padaptive=ComputePAdaptive(Pedges(:,ii),orig_a(:,ii),prevA,P1param);
    Ep=ComputeExpectedValue(prevA,prevB);
    [SaLR(:,ii),SbLR(:,ii),ScLR(:,ii)]=ComputePropagatedParabola(SaLR(:,ii),SbLR(:,ii),ScLR(:,ii),Padaptive,Ep);
    
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

    prevA=SaRL(:,ii+1);
    prevB=SbRL(:,ii+1);
    
    Padaptive=ComputePAdaptive(Pedges(:,ii),orig_a(:,ii),prevA,P1param);
    Ep=ComputeExpectedValue(prevA,prevB);
    [SaRL(:,ii),SbRL(:,ii),ScRL(:,ii)]=ComputePropagatedParabola(SaRL(:,ii),SbRL(:,ii),ScRL(:,ii),Padaptive,Ep);
    
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

