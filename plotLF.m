function plotLF(W,limit)
% изобрашение ассимптотической ЛАЧХ и ЛФЧХ
[A,F,lgw1,lgws] = freqasymp(W, limit, 1);
figure
    subplot(2,1,1)
        h(1) = plot(lgw1,A); grid
        xlabel('lg \omega'), ylabel('dB')
        ax = axis; kf = length(lgws);
        line(ones(2, kf)*diag(lgws),diag([ax(3) ax(4)])*ones(2,kf),...
        'Color', [1 0 0], 'LineStyle', '--')
    subplot(2,1,2)
        h(2) = plot(lgw1,F); grid
        xlabel('lg \omega')
        ylabel('phase (deg)')
        ax = axis;
        line(ones(2, kf)*diag(lgws),diag([ax(3) ax(4)])*ones(2,kf),...
        'Color', [1 0 0], 'LineStyle', '--')
        set(h,'LineWidth', 2.5)