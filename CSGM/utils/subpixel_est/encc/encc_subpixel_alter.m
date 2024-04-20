function [ tau, enccArray,max_map]=encc_subpixel_alter(L, R, window, maxDisp, flag_zm, flag_mode,dispar_map)

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
disp_map2=zeros(size(L)+[w1-1 w2-1]);
LL(mw1+1:end-mw1,mw2+1:end-mw2)=L;
RR(mw1+1:end-mw1,mw2+1:end-mw2)=R;
disp_map2(mw1+1:end-mw1,mw2+1:end-mw2)=dispar_map;
max_map = zeros(size(disp_map2));
L=LL;
R=RR;
dispar_map = disp_map2;
clear RR LL disp_map2;


sizeL=size(L);
disparity=zeros(sizeL);
enccArray=disparity; tau=enccArray;

N=disparity;
OP=disparity;
T=disparity;
N=N-1;
N2=N;



for i=mw1+1:sizeL(1)-mw1 % y-idx
%     tic;

    for j=mw2+1:sizeL(2)-mw2 %x-idx
        %reference window; we look for the conjugate one in the right image
        Ml = L(i-mw1:i+mw1,j-mw2:j+mw2);
        
        curr_disp = dispar_map(i,j);
        cor=-Inf*ones(1,4);
        tValues = zeros(1,4);
        maxDisp = abs(curr_disp)+1;
        if j-mw2-1>maxDisp % Iteration over disparites (if max disp is not limit)
            

            % precompute the cross-products between left and right intensities
            % along the disparity range 
%             crossFilter = zeros(4,1);
%             tic
            a = R(i-mw1:i+mw1, j-maxDisp-mw2-1:j+mw2-maxDisp+2-3);
            crossFilter(1) = Ml(:)'*a(:);
            a = R(i-mw1:i+mw1, j-maxDisp-mw2:j+mw2-maxDisp+2-2);
            crossFilter(2) = Ml(:)'*a(:);
            a = R(i-mw1:i+mw1, j-maxDisp-mw2+1:j+mw2-maxDisp+2-1);
            crossFilter(3) = Ml(:)'*a(:);
            a = R(i-mw1:i+mw1, j-maxDisp-mw2+2:j+mw2-maxDisp+2-0);
            crossFilter(4) = Ml(:)'*a(:);
            
%             
%             a = R(i-mw1:i+mw1, j-maxDisp-mw2:j+mw2-maxDisp+2-2);
%             crossFilter(1) = Ml(:)'*a(:);
%             a = R(i-mw1:i+mw1, j-maxDisp-mw2+1:j+mw2-maxDisp+2-1);
%             crossFilter(2) = Ml(:)'*a(:);
%             a = R(i-mw1:i+mw1, j-maxDisp-mw2+2:j+mw2-maxDisp+2-0);
%             crossFilter(3) = Ml(:)'*a(:);
%             
            
%             toc
%             tic
%             crossFilter2 = filter2(Ml, R(i-mw1:i+mw1, j-maxDisp-mw2-1:j+mw2-maxDisp+2), 'valid');
%             toc
            
%             tic
%             crossFilter = conv2(R(i-mw1:i+mw1, j-maxDisp-mw2-1:j+mw2-maxDisp+2),fliplr(flipud(Ml)), 'valid');
% 
%             toc
%             crossFilter = filter2(Ml, R(i-mw1:i+mw1, j-maxDisp-mw2-1:j+mw2), 'valid');
                        
%             for k = j+curr_disp:j+curr_disp+1 %k=j-maxDisp:j
            for k = j+curr_disp-1:j+curr_disp+1 %k=j-maxDisp:j
%                 Mr = R(i-mw1:i+mw1,k-mw2:k+mw2);
                
                % compute the left-adjacent window of Rwin; overlap is (w2-1) columns
%                 if (k-mw2>1)
%                     Mr_m=R(i-mw1:i+mw1,k-mw2-1:k+mw2-1);
%                 else
%                     Mr_m=[zeros(w1,1) Mr(:,1:end-1)];
%                 end
                
                % optimum solution 
                lambda = Lambdas(i-mw1,k-mw2);
                r = rMat(i-mw1,k-mw2);

                B1= (crossFilter(k-(j-maxDisp)+2)-Mol(i-mw1,j-mw2)*Mor(i-mw1,k-mw2))/(VaL(i-mw1,j-mw2)*VaR(i-mw1,k-mw2)+.00001);
                B2 = (crossFilter(k-(j-maxDisp)+1)-Mol(i-mw1,j-mw2)*Mor(i-mw1,k-mw2-1))/(VaL(i-mw1,j-mw2)*VaR(i-mw1,k-mw2-1)+0.00001);
%                 B1= (crossFilter(k-(j-maxDisp)+1)-Mol(i-mw1,j-mw2)*Mor(i-mw1,k-mw2))/(VaL(i-mw1,j-mw2)*VaR(i-mw1,k-mw2)+.00001);
%                 B2 = (crossFilter(k-(j-maxDisp))-Mol(i-mw1,j-mw2)*Mor(i-mw1,k-mw2-1))/(VaL(i-mw1,j-mw2)*VaR(i-mw1,k-mw2-1)+0.00001);
                
                
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
                
%                 cor(k-j+2+abs(curr_disp)-1)= encc_value;
%                 tValues(k-j+2+abs(curr_disp)-1)=t;
                cor(k-j+2+abs(curr_disp))= encc_value;
                tValues(k-j+2+abs(curr_disp))=t;
                
            end
%             tic
            [m,n]=max(cor);
%             if n ==1
%                 n = 2;%display('hi')
%             end
%             toc
%             n=1;
            max_map(i,j) = n;
%             enccArray(i-mw1,j-mw2)=m;
            
            tau(i-mw1,j-mw2)=tValues(n);
            
            
            clear cor;
%             disparity(i-mw1,j-mw2)=maxDisp-n+1;
            
%             if (flag_mode) %shiftable-window mode
%                 N2(i-mw1:i+mw1,j-mw2:j+mw2)=m;
%                 I=N2>N;
%                 N(I)=N2(I);
%                 OP(I)=maxDisp-n+1;
%                 T(I)=tValues(n);
%             end
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
            max_map(i,j) = n;
            enccArray(i-mw1,j-mw2)=m;
            tau(i-mw1,j-mw2)=tValues(n);
            
            
%             clear cor
%             disparity(i-mw1,j-mw2)=j-mw2-n;
%             
%             if (flag_mode)%shiftable-window mode
%                 
%                 N2(i-mw1:i+mw1,j-mw2:j+mw2)=m;
%                 I=N2>N;
%                 N(I)=N2(I);
%                 OP(I)=j-mw2-n;
%                 T(I)=tValues(n);
%             end
            
        end
       
    end
%      toc
end


if ~flag_mode
    disparity=disparity(1:end-w1+1,1:end-w2+1);
    enccArray=enccArray(1:end-w1+1,1:end-w2+1);
    tau=tau(1:end-w1+1,1:end-w2+1);
    max_map = max_map(mw1+1:end-mw1,mw2+1:end-mw2);

else
    disparity=OP(mw1+1:end-mw1,mw2+1:end-mw2);
    tau=T(mw1+1:end-mw1,mw2+1:end-mw2);
    enccArray=N(mw1+1:end-mw1,mw2+1:end-mw2);
    max_map = max_map(mw1+1:end-mw1,mw2+1:end-mw2);
end

 %final disparity refined by good subpixel corrections
% tau = abs(tau)<1.5;
tau(tau>1) = 0;
tau(tau<-2) = 0;

% disparity(tauMap) = disparity(tauMap)-tau(tauMap);
% 
% disparity(disparity<0) = 0;
% disparity(disparity>maxDisp) = maxDisp;

toc;

