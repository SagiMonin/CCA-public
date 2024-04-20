function [dispLeft,dispRight]=leftRightConsistency(dispLeft,dispRight,lrThresh)

[xs,ys]=meshgrid(1:size(dispLeft,2),1:size(dispLeft,1));
rightInLeft=interp2(xs,ys,dispRight,xs+dispLeft,ys,'linear',nan);
leftInRight=interp2(xs,ys,dispLeft, xs+dispRight,ys,'linear',nan);

dispLeft(isnan(rightInLeft) | abs(rightInLeft+dispLeft)>lrThresh)=nan;
dispRight(isnan(leftInRight) | abs(leftInRight+dispRight)>lrThresh)=nan;

