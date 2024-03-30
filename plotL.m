function plotL(W,limit)
% изобрашение ассимптотической ЛАЧХ
[L0,~,lgw1,lgws] = freqasymp(W, limit, 1);
figure
h = plot(lgw1,L0); grid
title('ЛАЧХ')
xlabel('lg \omega'), ylabel('dB')
ax = axis; kf = length(lgws);
line(ones(2, kf)*diag(lgws),diag([ax(3) ax(4)])*ones(2,kf),...
'Color', [1 0 0], 'LineStyle', '--')
set(h,'LineWidth', 2.5)