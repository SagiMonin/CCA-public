function [aNew,bNew,cNew]=ComputePropagatedParabola(a,b,c,P1,Ep)

% Propagation
aNew=a- P1;
bNew=b+  2.*P1.*Ep;
cNew=c - P1.*Ep.^2 ;
end