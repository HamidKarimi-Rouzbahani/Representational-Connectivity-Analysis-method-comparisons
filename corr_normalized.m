function [pXYz_normalized]=corr_normalized(X,Y,maxXY)

pXYn=corr(X,Y,'row','complete')./maxXY;

pXYz_normalized=pXYn;

end
