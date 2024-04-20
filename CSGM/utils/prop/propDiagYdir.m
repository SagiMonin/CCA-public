function parab_prop = propDiagYdir(im_guide, parab, params, dir_flag)
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

if dir_flag == 'MD' % Main diagonal
    SaLR = parab.a;
    SbLR = parab.b;
    ScLR = parab.c;
    SdLR = parab.d;
    SeLR = parab.e;
    SfLR = parab.f;
    orig_a = parab.a;

elseif dir_flag == 'OD' % Off-diagional
    SaLR = fliplr(parab.a);
    SbLR = fliplr(parab.b);
    ScLR = fliplr(parab.c);
    SdLR = fliplr(parab.d);
    SeLR = fliplr(parab.e);
    SfLR = fliplr(parab.f);
    orig_a = fliplr(parab.a);

end

for ii=2:size(SaLR,2)    
    % Top left
    prevATopLeft = [SaLR(1,ii); SaLR(1:end-1,ii-1)];
    prevBTopLeft = [SbLR(1,ii); SbLR(1:end-1,ii-1)];
    prevCTopLeft = [ScLR(1,ii); ScLR(1:end-1,ii-1)];
    prevDTopLeft = [SdLR(1,ii); SdLR(1:end-1,ii-1)];
    prevETopLeft = [SeLR(1,ii); SeLR(1:end-1,ii-1)];
    prevFTopLeft = [SfLR(1,ii); SfLR(1:end-1,ii-1)];
    
    
    [P1adaptive] = ComputePAdaptiveYdir(PedgesTopLeft(:,ii),P1param);
%     [Epx, Epy] = getXYPrabolodid(prevA,prevB,prevC,prevD,prevE,prevF);
    [SaLR(:,ii),SbLR(:,ii),ScLR(:,ii),SdLR(:,ii),SeLR(:,ii),SfLR(:,ii)] = ...
        ComputePropagatedParabolaYdir(SaLR(:,ii),SbLR(:,ii),ScLR(:,ii),SdLR(:,ii),SeLR(:,ii),SfLR(:,ii),...
                                      prevATopLeft,prevBTopLeft,prevCTopLeft,prevDTopLeft,prevETopLeft,prevFTopLeft,P1adaptive);
                                  


   
%     % Diagonal aggegation
% %     [P1adaptive, P2adaptive] = ComputePAdaptiveYdir(PedgesTopLeft(:,ii),prevATopLeft,prevDTopLeft,P1param,P2param);
%     [Epx, Epy] = getXYPraboloid(prevATopLeft,prevBTopLeft,prevCTopLeft,prevDTopLeft,prevETopLeft,prevFTopLeft);
%     [SaLR(:,ii),SbLR(:,ii),ScLR(:,ii),SdLR(:,ii),SeLR(:,ii),SfLR(:,ii)] = ...
%         ComputePropagatedParabolaYdir(SaLR(:,ii),SbLR(:,ii),ScLR(:,ii),SdLR(:,ii),SeLR(:,ii),SfLR(:,ii),P1adaptive,P2adaptive,Epx,Epy);
end

if dir_flag == 'MD' % Main diagonal
    SaRL = parab.a;
    SbRL = parab.b;
    ScRL = parab.c;
    SdRL = parab.d;
    SeRL = parab.e;
    SfRL = parab.f;
    orig_a = parab.a;
elseif dir_flag == 'OD' % Off-diagional
    SaRL = fliplr(parab.a);
    SbRL = fliplr(parab.b);
    ScRL = fliplr(parab.c);
    SdRL = fliplr(parab.d);
    SeRL = fliplr(parab.e);
    SfRL = fliplr(parab.f);
    orig_a = fliplr(parab.a);
end

% df=sum(im_guide(2:end,2:end,:)-im_guide(1:end-1,1:end-1,:),3)./size(im_guide,3);
df=sum(im_guide(1:end-1,1:end-1,:)-im_guide(2:end,2:end,:),3)./size(im_guide,3);
PedgesBotLeft=exp(-df.^2./sigmaEdges^2);
PedgesBotLeft=[ PedgesBotLeft ones([size(PedgesBotLeft,1) 1])];
PedgesBotLeft=[ PedgesBotLeft;ones([1 size(PedgesBotLeft,2)])];

for ii=size(SaRL,2)-1:-1:1
    % Bot left
    prevABotLeft=[ SaRL(2:end,ii+1);SaRL(end,ii) ];
    prevBBotLeft=[ SbRL(2:end,ii+1);SbRL(end,ii) ];
    prevCBotLeft=[ ScRL(2:end,ii+1);ScRL(end,ii) ];
    prevDBotLeft=[ SdRL(2:end,ii+1);SdRL(end,ii) ];
    prevEBotLeft=[ SeRL(2:end,ii+1);SeRL(end,ii) ];
    prevFBotLeft=[ SfRL(2:end,ii+1);SfRL(end,ii) ];
    
    [P1adaptive] = ComputePAdaptiveYdir(PedgesBotLeft(:,ii),P1param);
%     [Epx, Epy] = getXYPrabolodid(prevA,prevB,prevC,prevD,prevE,prevF);
    [SaRL(:,ii),SbRL(:,ii),ScRL(:,ii),SdRL(:,ii),SeRL(:,ii),SfRL(:,ii)] = ...
        ComputePropagatedParabolaYdir(SaRL(:,ii),SbRL(:,ii),ScRL(:,ii),SdRL(:,ii),SeRL(:,ii),SfRL(:,ii),...
                                      prevABotLeft,prevBBotLeft,prevCBotLeft,prevDBotLeft,prevEBotLeft,prevFBotLeft,P1adaptive);
    
    
    % Diagonal aggegation
%     [P1adaptive, P2adaptive] = ComputePAdaptiveYdir(PedgesBotLeft(:,ii),prevABotLeft,prevDBotLeft,P1param,P2param);
%     [Epx, Epy] = getXYPraboloid(prevABotLeft,prevBBotLeft,prevCBotLeft,prevDBotLeft,prevEBotLeft,prevFBotLeft);
%     [SaRL(:,ii),SbRL(:,ii),ScRL(:,ii),SdRL(:,ii),SeRL(:,ii),SfRL(:,ii)] = ...
%         ComputePropagatedParabolaYdir(SaRL(:,ii),SbRL(:,ii),ScRL(:,ii),SdRL(:,ii),SeRL(:,ii),SfRL(:,ii),P1adaptive,P2adaptive,Epx,Epy);   
end


% for ii=2:size(SaLR,2)
%     % Bot left
%     prevABotLeft=[ SaRL(2:end,ii-1);SaRL(end,ii) ];
%     prevBBotLeft=[ SbRL(2:end,ii-1);SbRL(end,ii) ];
%     prevCBotLeft=[ ScRL(2:end,ii-1);ScRL(end,ii) ];
%     prevDBotLeft=[ SdRL(2:end,ii-1);SdRL(end,ii) ];
%     prevEBotLeft=[ SeRL(2:end,ii-1);SeRL(end,ii) ];
%     prevFBotLeft=[ SfRL(2:end,ii-1);SfRL(end,ii) ];
%     
%     [P1adaptive] = ComputePAdaptiveYdir(PedgesBotLeft(:,ii),P1param);
% %     [Epx, Epy] = getXYPrabolodid(prevA,prevB,prevC,prevD,prevE,prevF);
%     [SaRL(:,ii),SbRL(:,ii),ScRL(:,ii),SdRL(:,ii),SeRL(:,ii),SfRL(:,ii)] = ...
%         ComputePropagatedParabolaYdir(SaRL(:,ii),SbRL(:,ii),ScRL(:,ii),SdRL(:,ii),SeRL(:,ii),SfRL(:,ii),...
%                                       prevABotLeft,prevBBotLeft,prevCBotLeft,prevDBotLeft,prevEBotLeft,prevFBotLeft,P1adaptive);
%     
%     
%     % Diagonal aggegation
% %     [P1adaptive, P2adaptive] = ComputePAdaptiveYdir(PedgesBotLeft(:,ii),prevABotLeft,prevDBotLeft,P1param,P2param);
% %     [Epx, Epy] = getXYPraboloid(prevABotLeft,prevBBotLeft,prevCBotLeft,prevDBotLeft,prevEBotLeft,prevFBotLeft);
% %     [SaRL(:,ii),SbRL(:,ii),ScRL(:,ii),SdRL(:,ii),SeRL(:,ii),SfRL(:,ii)] = ...
% %         ComputePropagatedParabolaYdir(SaRL(:,ii),SbRL(:,ii),ScRL(:,ii),SdRL(:,ii),SeRL(:,ii),SfRL(:,ii),P1adaptive,P2adaptive,Epx,Epy);   
% end

if dir_flag == 'MD' % Main diagonal
    parab_prop.a = SaLR+SaRL;
    parab_prop.b = SbLR+SbRL;
    parab_prop.c = ScLR+ScRL;
    parab_prop.d = SdLR+SdRL;
    parab_prop.e = SeLR+SeRL;
    parab_prop.f = SfLR+SfRL;
elseif dir_flag == 'OD' % Off-diagional
    parab_prop.a = fliplr(SaLR+SaRL);
    parab_prop.b = fliplr(SbLR+SbRL);
    parab_prop.c = fliplr(ScLR+ScRL);
    parab_prop.d = fliplr(SdLR+SdRL);
    parab_prop.e = fliplr(SeLR+SeRL);
    parab_prop.f = fliplr(SfLR+SfRL);
end



end