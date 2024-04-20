function parab_prop = propDiag(im_guide, parab, params, dir_flag)
% Extract values 
P1param = params.P1param;
sigmaEdges = params.sigmaEdges;

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
    prevATopLeft=[ SaLR(1,ii); SaLR(1:end-1,ii-1)];
    prevBTopLeft=[ SbLR(1,ii); SbLR(1:end-1,ii-1)];
   
    % Diagonal aggegation
    Padaptive=ComputePAdaptive(PedgesTopLeft(:,ii),orig_a(:,ii),prevATopLeft,P1param);
    Ep=ComputeExpectedValue(prevATopLeft,prevBTopLeft);
    [SaLR(:,ii),SbLR(:,ii),ScLR(:,ii)]=ComputePropagatedParabola(SaLR(:,ii),SbLR(:,ii),ScLR(:,ii),Padaptive,Ep);
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
    prevABotLeft=[ SaRL(2:end,ii-1);SaRL(end,ii) ];
    prevBBotLeft=[ SbRL(2:end,ii-1);SbRL(end,ii) ];
    
    Padaptive=ComputePAdaptive(PedgesBotLeft(:,ii),orig_a(:,ii),prevABotLeft,P1param);
    Ep=ComputeExpectedValue(prevABotLeft,prevBBotLeft);
    [SaRL(:,ii),SbRL(:,ii),ScRL(:,ii)]=ComputePropagatedParabola(SaRL(:,ii),SbRL(:,ii),ScRL(:,ii),Padaptive,Ep);
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