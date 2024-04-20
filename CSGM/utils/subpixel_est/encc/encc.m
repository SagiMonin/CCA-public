function [disparity, tau, enccArray]=encc(L, R, window, maxDisp, flag_zm, flag_mode)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENCC Enhanced Normalized Cross Correlation for Stereo Correspondence
%
% This function implements the ENCC stereo correspondence algorithm as it
% is presented in the paper "An Enhanced Correlation-based Method for
% Stereo Corresponce with Subpixel Accuracy" by E.Z. Psarakis and G. D.
% Evangelidis presented in ICCV-2005.
%
% [DISPARITY, TAU, ENCCARRAY]= ENCC(L, R, WINDOW, MAXDISP, FLAG_ZM, FLAG_MODE)
%
% Given two stereo (left-right) images, L and R, it returns in DISPARITY
% the disparity map with subpixel accuracy. TAU is the subpixel correction
% that is already embodied in DISPARITY, but is returned for other purposes
% (i.e. to detect sticky points). ENCCARRAY is the array with the optimum
% ENCC value for each left investigated pixel (block).
% -------------------------------------------------------------------------
% Input arguments:
%     L: left stereo image
%     R: right stereo image
%     WINDOW: a 1x2 vector that defines the size ([rows, columns])
%         of the window (block)
%     MAXDISP: maximum disparity (disparity range: [0, MAXDISP])
%     FLAG_ZM: zero-mean flag.
%         By FLAG_ZM=1, each block reduces to zero-mean (by subtracting
%         the average) before its processing, while FLAG_ZM =0 leaves
%         block as is.
%     FLAG_MODE: implementation-mode flag.
%         When FLAG_MODE = 0 the typical local strategy is followed
%         (each block decides about its center), while FLAG_MODE = 1
%         enables the shiftable-window mode where each point is
%         mathed with respect to the best window that participate into
%         (among WINDOW(1)*WINDOW(2) windows), no matter what its position
%         inside the window.
%  -----------------------------------------------------------------------
% $ Ver: 2.1, 5/6/2012,  released by Georgios D. Evangelidis.
% For any comment, please contact the author
% e-mail: evagelid@ceid.upatras.gr
%
% Note that this software is provided "as is" without any kind of warranty.
% Still, it is provided for research purposes only. In any case, please
% cite the above paper.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% check for multichannel images
if size(L,3)>1
    if ((size(L,3)==2) || (size(L,3)>3))
        error('ENCC: Unknown color format: check the number of channels');
    else
        if strcmp(class(L),'uint8')
        L=rgb2gray(L);
        elseif strcmp(class(L),'double')
                L = L * (255/max(L(:)));
                L=rgb2gray(uint8(L));
        else
            error('ENCC: Unknown type of images (convert images to double or uint8 class)');
        end
        
        
    end
end
if size(R,3)>1
    if ((size(R,3)==2) || (size(R,3)>3))
        error('ENCC: Unknown color format: check the number of channels');
    else
        if strcmp(class(R),'uint8')
        R=rgb2gray(R);
        elseif strcmp(class(R),'double')
                R = R *(255/max(R(:)));
                R=rgb2gray(uint8(R));
        else
            error('ENCC: Unknown type of images (convert images to double or uint8 class)');
        end
        
        
    end
end


L=double(L);
R=double(R);
w=window;


R=double(R);
L=double(L);

szLi=size(L);szRi=size(R);

if (szLi(1)~=szRi(1)) || (szLi(2)~=szRi(2))
    error('ENCC: Images must have the same resolution.');
end


w1=w(1);w2=w(2);
mw1=(w1-1)/2;mw2=(w2-1)/2;

% compute the average of each window
h=fspecial('average',w);

% precompute local characteristics
Mol = imfilter(L, ones(w1,w2), 'same')/sqrt(w1*w2);
Mor = imfilter(R, ones(w1,w2), 'same')/sqrt(w1*w2);

if flag_zm == 0
    Mol = 0*Mol;
    Mor = 0*Mor;
end

VarL = imfilter(L.^2, ones(w1,w2), 'same');
VarR = imfilter(R.^2, ones(w1,w2), 'same');

VaL = sqrt(VarL - Mol.^2);
VaR = sqrt(VarR - Mor.^2);

%ratio of norms between two concecutive windows in the right image
Lambdas = [ones(size(VaR,1),1)  VaR(:,1:end-1)./(VaR(:,2:end)+.000001) ];

% %auto-correlation
R2R = R.*([ones(size(L,1),1) R(:,1:end-1)]);
autoC = imfilter(R2R, ones(w1,w2), 'same')-Mor.*([ones(size(L,1),1) Mor(:,1:end-1)]);
V2V = [ones(size(VaR,1),1)  VaR(:,1:end-1).*(VaR(:,2:end)) ];
rMat = autoC./(V2V+.000001);


%zero padding
LL=zeros(size(L)+[w1-1 w2-1]);
RR=zeros(size(R)+[w1-1 w2-1]);
LL(mw1+1:end-mw1,mw2+1:end-mw2)=L;
RR(mw1+1:end-mw1,mw2+1:end-mw2)=R;
L=LL;
R=RR;
clear RR LL


sizeL=size(L);
disparity=zeros(sizeL);
enccArray=disparity; tau=enccArray;

N=disparity;
OP=disparity;
T=disparity;
N=N-1;
N2=N;

tic


for i=mw1+1:sizeL(1)-mw1 % y-idx
    for j=mw2+1:sizeL(2)-mw2 %x-idx
        %reference window; we look for the conjugate one in the right image
        Ml = L(i-mw1:i+mw1,j-mw2:j+mw2);
        
        cor=-Inf*ones(1,maxDisp+1);
        tValues = zeros(1,maxDisp+1);
        
        if j-mw2-1>maxDisp % Iteration over disparites (if max disp is not limit)
            

            % precompute the cross-products between left and right intensities
            % along the disparity range 
            crossFilter = filter2(Ml, R(i-mw1:i+mw1, j-maxDisp-mw2-1:j+mw2), 'valid');
                        
            for k=j-maxDisp:j
                Mr = R(i-mw1:i+mw1,k-mw2:k+mw2);
                
                % compute the left-adjacent window of Rwin; overlap is (w2-1) columns
                if (k-mw2>1)
                    Mr_m=R(i-mw1:i+mw1,k-mw2-1:k+mw2-1);
               else
                    Mr_m=[zeros(w1,1) Mr(:,1:end-1)];
                end
                
                % optimum solution 
                lambda = Lambdas(i-mw1,k-mw2);
                r = rMat(i-mw1,k-mw2);

                B1= (crossFilter(k-(j-maxDisp)+2)-Mol(i-mw1,j-mw2)*Mor(i-mw1,k-mw2))/(VaL(i-mw1,j-mw2)*VaR(i-mw1,k-mw2)+.00001);
                B2 = (crossFilter(k-(j-maxDisp)+1)-Mol(i-mw1,j-mw2)*Mor(i-mw1,k-mw2-1))/(VaL(i-mw1,j-mw2)*VaR(i-mw1,k-mw2-1)+0.00001);
                
                
                nom = B2-r*B1;
                den = lambda*(r*B2-B1)-nom;
                
                % subpixel correction by maximizing ENCC
                t = nom/den;
                
                % denominator's sign
                a_b=sign(den);
                
                % ENCC value
                encc_value = sqrt((B1^2+B2^2-2*r*B1*B2)/(1-r^2));
                
                if a_b>0
                    t=0;
                    encc_value=B1;
                end
                
                cor(k-j+maxDisp+1)= encc_value;
                tValues(k-j+maxDisp+1)=t;
                
            end
            
            [m,n]=max(cor);

            enccArray(i-mw1,j-mw2)=m;
            
            tau(i-mw1,j-mw2)=tValues(n);
            
            
            clear cor;
            disparity(i-mw1,j-mw2)=maxDisp-n+1;
            
            if (flag_mode) %shiftable-window mode
                N2(i-mw1:i+mw1,j-mw2:j+mw2)=m;
                I=N2>N;
                N(I)=N2(I);
                OP(I)=maxDisp-n+1;
                T(I)=tValues(n);
            end
        else
            for k=mw2+1:j % Iterate over disparities
                Mr = R(i-mw1:i+mw1,k-mw2:k+mw2);
                
                % compute the left-adjacent window of Rwin; the overlap is
                % (w2-1) columns
                if (k-mw2>1)
                    Mr_m = R(i-mw1:i+mw1,k-mw2-1:k+mw2-1);
                else
                    Mr_m=[zeros(w1,1) Mr(:,1:end-1)];
                end
                
                if k-mw2-1<1;
                    lambda = 1;
                    r = 1;
                    B1 = 1;
                    B2 = 1;
                    
                else
                    
                    %optimum solution
                    lambda = Lambdas(i-mw1,k-mw2);
                    r = rMat(i-mw1,k-mw2);
                    B1= (sum(Ml(:).*Mr(:))-Mol(i-mw1,j-mw2)*Mor(i-mw1,k-mw2))/(VaL(i-mw1,j-mw2)*VaR(i-mw1,k-mw2)+.00001);
                    B2= (sum(Ml(:).*Mr_m(:))-Mol(i-mw1,j-mw2)*Mor(i-mw1,k-mw2-1))/(VaL(i-mw1,j-mw2)*VaR(i-mw1,k-mw2-1)+.00001);

                end

                
                nom = B2-r*B1;
                den = lambda*(r*B2-B1)-nom;
                
                % subpixel correction by maximizing ENCC
                t = nom/den;
                
                % denominator's sign
                a_b=sign(den);
                
                % ENCC value
                encc_value = sqrt((B1^2+B2^2-2*r*B1*B2)/(1-r^2));
                
                if a_b>0
                    t=0;
                    encc_value=B1;
                end
                
                cor(k-mw2)=encc_value;
                tValues(k-mw2)=t;
                
            end
            
            [m,n]=max(cor);
            
            enccArray(i-mw1,j-mw2)=m;
            tau(i-mw1,j-mw2)=tValues(n);
            
            
            clear cor
            disparity(i-mw1,j-mw2)=j-mw2-n;
            
            if (flag_mode)%shiftable-window mode
                
                N2(i-mw1:i+mw1,j-mw2:j+mw2)=m;
                I=N2>N;
                N(I)=N2(I);
                OP(I)=j-mw2-n;
                T(I)=tValues(n);
            end
            
        end
    end
end


if ~flag_mode
    disparity=disparity(1:end-w1+1,1:end-w2+1);
    enccArray=enccArray(1:end-w1+1,1:end-w2+1);
    tau=tau(1:end-w1+1,1:end-w2+1);
else
    disparity=OP(mw1+1:end-mw1,mw2+1:end-mw2);
    tau=T(mw1+1:end-mw1,mw2+1:end-mw2);
    enccArray=N(mw1+1:end-mw1,mw2+1:end-mw2);
end

 %final disparity refined by good subpixel corrections
tauMap = abs(tau)<1.5;
disparity(tauMap) = disparity(tauMap)-tau(tauMap);

disparity(disparity<0) = 0;
disparity(disparity>maxDisp) = maxDisp;

toc;

