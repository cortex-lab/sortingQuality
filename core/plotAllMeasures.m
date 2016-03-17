


function plotAllMeasures(cgs, uQ, cR, isiV)

figure; 

subplot(3,1,1);
plotByCG(cgs, uQ, cR)
xlabel('iso distance'); ylabel('contamination from mahal');

subplot(3,1,2);
plotByCG(cgs, uQ, isiV)
xlabel('iso distance'); ylabel('contamination from ISI');

subplot(3,1,3);
plotByCG(cgs, cR, isiV)
xlabel('contamination from mahal'); ylabel('contamination from ISI');

function plotByCG(cgs, X, Y)
plot(X(cgs==0), Y(cgs==0), 'k.')
hold on; plot(X(cgs==1), Y(cgs==1), 'ro')
hold on; plot(X(cgs==2), Y(cgs==2), 'go')
hold on; plot(X(cgs==3), Y(cgs==3), 'm.')
[n,~] = hist(cgs,0:3);
legStrings = {'noise', 'mua', 'good', 'unsorted'};
legend(legStrings(n>0));