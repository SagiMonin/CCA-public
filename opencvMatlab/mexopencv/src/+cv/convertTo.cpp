/**
 * @file convertTo.cpp
 * @brief mex interface for cv::Mat::convertTo
 * @ingroup core
 * @author Amro
 * @date 2016
 */
#include "mexopencv.hpp"
using namespace std;
using namespace cv;

/**
 * Main entry called from Matlab
 * @param nlhs number of left-hand-side arguments
 * @param plhs pointers to mxArrays in the left-hand-side
 * @param nrhs number of right-hand-side arguments
 * @param prhs pointers to mxArrays in the right-hand-side
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Check the number of arguments
    nargchk(nrhs>=1 && (nrhs%2)==1 && nlhs<=1);

    // Argument vector
    vector<MxArray> rhs(prhs, prhs+nrhs);

    // Option processing
    int rtype = -1;
    double alpha = 1.0;
    double beta = 0.0;
    for (int i=1; i<nrhs; i+=2) {
        string key(rhs[i].toString());
        if (key == "RType")
            rtype = (rhs[i+1].isChar()) ?
                ClassNameMap[rhs[i+1].toString()] : rhs[i+1].toInt();
        else if (key == "Alpha")
            alpha = rhs[i+1].toDouble();
        else if (key == "Beta")
            beta = rhs[i+1].toDouble();
        else
            mexErrMsgIdAndTxt("mexopencv:error",
                "Unrecognized option %s", key.c_str());
    }

    // Process
    Mat src(rhs[0].toMat()), dst;
    src.convertTo(dst, rtype, alpha, beta);
    plhs[0] = MxArray(dst);
}
