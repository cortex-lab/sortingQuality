


function plotAllMeasures(cgs, uQ, cR, isiV)

figure; 

subplot(3,1,1);
plotByCG(cgs, uQ, cR)

subplot(3,1,2);
plotByCG(cgs, uQ, isiV)

subplot(3,1,3);
plotByCG(cgs, cR, isiV)

function plotByCG(cgs, X, Y)
plot(X(cgs==0), Y(cgs==0), 'k.')
hold on; plot(X(cgs==1), T(cgs==1), 'ro')
hold on; plot(X(cgs==2), Y(cgs==2), 'go')
legend({'noise', 'mua', 'good'})