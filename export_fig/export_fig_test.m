% Mengchen Zhu.
clc
clear all
close all

fig1 = rand(20,20);
fig2 = rand(20,20);

h = figure;
export_fig_fun(fig1, fig2, h);
set(h, 'Color', 'w')

export_fig(h, './export_fig_subplot.eps');

close(h)