// (C) Georgios Evangelidis
// e-mail: evagelid@ceid.upatras.gr
// ver. 2.0, 28/5/2012

#include <matrix.h>
#include <mex.h>
#include "common_functions.cpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    if (nrhs!=6){
        mexErrMsgTxt("Function needs six input arguments:\n two image arrays \n window size \n disprity range\n two flags");
    }
    
    if ((mxGetNumberOfDimensions(prhs[0]) >2) || ((mxGetNumberOfDimensions(prhs[1]) >2))) {
        mexErrMsgTxt("Images must be grayscale (single-channel).");
    }
    
    double *dispPr, *Wpr, *LPr, *RPr, *rangePr, *flag2Pr, *enccmatPr;
    double *flagPr, *corPr, *cor2Pr, *accTPr, *accT2Pr, *tauPr;
    int padH, padW, Lw, Lh, Rw, Rh, Wx, Wy;
    mxArray *Lpad, *Rpad, *Lwin, *Rwin, *RwinM, *Disp, *Cor;
    mxArray *Cor2, *Tau, *accT, *accT2, *temp, *temp2, *bin, *Disp2, *Tau2;
    mxArray *ENCCmat2, *ENCCmat;
    
    int range, flag, i, j, k, flag2;
    bool shift_win_mode;//shiftable-window mode
    
    Lh = mxGetM(prhs[0]);
    Lw = mxGetN(prhs[0]);
    
    Rh = mxGetM(prhs[1]);
    Rw = mxGetN(prhs[1]);
    
    if ((Lw != Rw) || (Lh != Rh)) {
        mexErrMsgTxt("Stereo Images must have the same size.");
    }
    
    Wpr = mxGetPr(prhs[2]);
    Wy = (int)Wpr[0];
    Wx = (int)Wpr[1];
    
    rangePr = mxGetPr(prhs[3]);
    range = (int) rangePr[0];//maximum disparity (disparity range: [0, range])
    flagPr = mxGetPr(prhs[4]);
    flag = (int) flagPr[0];//flag=1 implies zero-mean ENCC
    
    if (range>Rw) {
        mexErrMsgTxt("Range cannot exceed the width of images.");
    }
    
    //flag2=1 enables the mode with shiftable windows, while flag2=0 implies
    //the classic mode where each block decides about its center only
    flag2Pr = mxGetPr(prhs[5]);
    flag2 = (int) flag2Pr[0];
    
    if (flag2==1) {shift_win_mode = true;}
    else {shift_win_mode = false;}
    
    if ((flag != 0) && (flag != 1)) {
        mexErrMsgTxt("Zero-mean flag takes values from {0,1}.");
    }
    
    if ((flag2 != 0) && (flag2 != 1)) {
        mexErrMsgTxt("Shiftable-window-mode flag takes values from {0,1}.");
    }
    
    if((!mxIsDouble(prhs[0])) || (!mxIsDouble(prhs[1]))){
        mexErrMsgTxt("Images must be of type double.");
    }
    
    padH = (Wy-1)/2;
    padW = (Wx-1)/2;
    
    Disp = plhs[0] = mxCreateDoubleMatrix(Lh, Lw, mxREAL);
    dispPr = mxGetPr(Disp);
    
    Tau = plhs[1] = mxCreateDoubleMatrix(Lh, Lw, mxREAL);
    tauPr = mxGetPr(Tau);
    
    ENCCmat = plhs[2] = mxCreateDoubleMatrix(Lh, Lw, mxREAL);
    enccmatPr = mxGetPr(ENCCmat);
    
    if (shift_win_mode){
        Disp2 = mxCreateDoubleMatrix(Lh+2*padH, Lw+2*padW, mxREAL);
        Tau2 = mxCreateDoubleMatrix(Lh+2*padH, Lw+2*padW, mxREAL);
        ENCCmat2 = mxCreateDoubleMatrix(Lh+2*padH, Lw+2*padW, mxREAL);
    }
    
    Lpad = mxCreateDoubleMatrix(Lh+2*padH, Lw+2*padW, mxREAL);
    Rpad = mxCreateDoubleMatrix(Lh+2*padH, Lw+2*padW, mxREAL);
    
    zeroPadding( prhs[0], Lpad, padW, padH);
    zeroPadding( prhs[1], Rpad, padW, padH);
    
    LPr =   mxGetPr(Lpad);
    RPr =   mxGetPr(Rpad);
    
    Lwin = mxCreateDoubleMatrix(Wy, Wx, mxREAL);
    Rwin = mxCreateDoubleMatrix(Wy, Wx, mxREAL);
    RwinM = mxCreateDoubleMatrix(Wy, Wx, mxREAL);
    
    Cor = mxCreateDoubleMatrix(range+1, 1, mxREAL);
    corPr = mxGetPr(Cor);
    
    accT = mxCreateDoubleMatrix(range+1, 1, mxREAL);
    accTPr = mxGetPr(accT);
    
    //matrices for shift-mode
    temp = mxCreateDoubleMatrix(Lh+2*padH, Lw+2*padW, mxREAL);
    temp2 = mxCreateDoubleMatrix(Lh+2*padH, Lw+2*padW, mxREAL);
    bin = mxCreateDoubleMatrix(Lh+2*padH, Lw+2*padW, mxREAL);
    
    //initialization to -1.
    setZero(temp);
    subArrS(temp, 1.0);
    setZero(temp2);
    subArrS(temp2, 1.0);
    
    double rho, ar, t, normR, normRM, lambda, rho1, nom, den, encc;
    int index;
    
    for (i=0; i<(Lw); i++){
        for (j=0; j<(Lh); j++){
            
            if (flag==0){
                getPatch(Lpad, Lwin, j, i);
            }
            else {
                getPatchZM(Lpad, Lwin, j, i);
            }
      
            normalizeMatrix(Lwin, Lwin);
            
            if (i>range){
                
                for (k=i-range; k<i+1; k++){
                    
                    if (flag==0) {
                        getPatch(Rpad, Rwin, j, k);
                        getPatch(Rpad, RwinM, j, k-1);
                    }
                    else {
                        getPatchZM(Rpad, Rwin, j, k);
                        getPatchZM(Rpad, RwinM, j, k-1);
                    }
                    
                    normR = normOfArr(Rwin);
                    normRM = normOfArr(RwinM);
                    
                    lambda = normRM/normR;
                    normalizeMatrix(Rwin, Rwin);
                    normalizeMatrix(RwinM, RwinM);
                    
                    rho = dot2d(Lwin, Rwin);
                    ar = dot2d(Rwin, RwinM);
                    rho1 = dot2d(Lwin, RwinM);
                    
                    nom = rho1-(ar*rho);
                    den = lambda*((ar*rho1)-rho)-nom;
                    
                   // subpixel correction by maximizing ENCC
                    t = nom/den;
                    
                    // ENCC value
                    encc = sqrt((rho*rho+rho1*rho1-2*ar*rho*rho1)/(1-ar*ar));
                   
                    if (den>0){
                        // positive sign for t's denominator => ENCC cannot be
                        // maximized => the best value we can achieve is rho.
                        t=0.0;
                        encc=rho;
                    }
                    
                    corPr[i-k] = encc;
                    accTPr[i-k] = t;
                    
                }
                
                index = maxIndexOfArr(Cor);
                
                
                if (!shift_win_mode){
                    dispPr[i*Lh+j]=index;
                    tauPr[i*Lh+j]=accTPr[index];
                    enccmatPr[i*Lh+j] = corPr[index];
                    
                }
                else {
                    
                    fillROI(temp2, j, i, Wx, Wy, corPr[index]);
                    compArr(temp2, temp, bin); //I need something faster here
                    equalArr(temp2, bin, temp);
                    fillArr(Disp2, bin, (double) index);
                    fillArr(Tau2, bin, accTPr[index]);
                    fillArr(ENCCmat2, bin, corPr[index]);
                }
                
            }
            else {
                
                Cor2 = mxCreateDoubleMatrix(i+1, 1, mxREAL);
                cor2Pr = mxGetPr(Cor2);
                cor2Pr[i]=-1;
                accT2 = mxCreateDoubleMatrix(i+1, 1, mxREAL);
                accT2Pr = mxGetPr(accT2);
                accT2Pr[i]=0;
                
                for (k=1; k< (i+1); k++){
                    
                    if (flag==0) {
                        getPatch(Rpad, Rwin, j, k);
                        getPatch(Rpad, RwinM, j, k-1);
                    }
                    else {
                        getPatchZM(Rpad, Rwin, j, k);
                        getPatchZM(Rpad, RwinM, j, k-1);
                    }
                    
                    normR = normOfArr(Rwin);
                    normRM = normOfArr(RwinM);
                    
                    lambda = normRM/normR;
                    
                    normalizeMatrix(Rwin, Rwin);
                    normalizeMatrix(RwinM, RwinM);
                    
                    rho = dot2d(Lwin, Rwin);
                    ar = dot2d(Rwin, RwinM);
                    rho1 = dot2d(Lwin, RwinM);
                    
                    nom = rho1-(ar*rho);
                    den = lambda*((ar*rho1)-rho)-nom;
                    
                    // subpixel correction by maximizing ENCC
                    t = nom/den;
                    
                    // ENCC value
                    encc = sqrt((rho*rho+rho1*rho1-2*ar*rho*rho1)/(1-ar*ar));
                    
                    if (den>0){
                        //positive sign for t's denominator => ENCC cannot be
                        //maximized => the best value you can achieve is rho.
                        t=0.0;
                        encc=rho;
                    }
                    
                    cor2Pr[i-k] = encc;
                    accT2Pr[i-k] = t;
                }
                    
                index = maxIndexOfArr(Cor2);


                if (!shift_win_mode){
                    dispPr[i*Lh+j] = index;
                    tauPr[i*Lh+j] = accT2Pr[index];
                    enccmatPr[i*Lh+j] = cor2Pr[index];
                    
                }
                else {
                    
                    fillROI(temp2, j, i, Wx, Wy, cor2Pr[index]);
                    compArr(temp2, temp, bin);//I need something faster here
                    equalArr(temp2, bin, temp);
                    fillArr(Disp2, bin, (double) index);
                    fillArr(Tau2, bin, accT2Pr[index]);
                    fillArr(ENCCmat2, bin, cor2Pr[index]);
                    
                }
                
                mxDestroyArray(Cor2);
                mxDestroyArray(accT2);
                
            }
        }
        
    }
    
  
    if (shift_win_mode){
        getPatch(Disp2, Disp, padH, padW);
        getPatch(Tau2, Tau, padH, padW);
        getPatch(ENCCmat2, ENCCmat, padH, padW);
        mxDestroyArray(Disp2);
        mxDestroyArray(ENCCmat2);
        mxDestroyArray(Tau2);
    }
    
    // refine disparities with all tau's (not good for any point)
    //subArr(Disp, Tau, Disp);
    
    
    
    // refine disparities with subpixel estimates when abs(tau)<1.5;
    for (i=0; i<(Lw); i++){
        for (j=0; j<(Lh); j++){
            if (fabs(tauPr[i*Lh+j])<1.5){
                dispPr[i*Lh+j] -= tauPr[i*Lh+j];
            }
            
            if (dispPr[i*Lh+j]<0){ //ignore negative disparities
                dispPr[i*Lh+j] = 0;
            }
            
            if (dispPr[i*Lh+j]>range){ //check for overflow
                dispPr[i*Lh+j] = range;
            }
            
        }
    }
    
    
    mxDestroyArray(Lpad);
    mxDestroyArray(Rpad);
    mxDestroyArray(Lwin);
    mxDestroyArray(Rwin);
    mxDestroyArray(RwinM);
    mxDestroyArray(Cor);
    mxDestroyArray(accT);
    mxDestroyArray(temp);
    mxDestroyArray(temp2);
    mxDestroyArray(bin);
    
}


