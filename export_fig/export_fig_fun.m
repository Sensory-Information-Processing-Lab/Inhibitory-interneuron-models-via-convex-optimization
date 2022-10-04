function export_fig_fun(fig1, fig2, h)
% Testing export_fig
% Mengchen Zhu
figure(h)

subplot(2,1,1)
imagesc(fig1);
axis image
axis off

subplot(2,1,2)
imagesc(fig2);
axis image
axis off

