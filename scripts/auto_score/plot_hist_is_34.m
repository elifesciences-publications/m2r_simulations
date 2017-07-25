function plot_hist_is_34(hist_4bin_is, hist_3bin_is, caption, range_bins)
if ~exist('range_bins','var') || isempty(range_bins); range_bins = 0:0.05:1; end;
range_bins = range_bins + 0.025;

y_lim_4 = max(hist_4bin_is(:)) + 5 - mod(max(hist_4bin_is(:)), 5);
y_lim_3 = max(hist_3bin_is(:)) + 5 - mod(max(hist_4bin_is(:)), 5);

figure(); set_print_page(gcf,0);
subplot(4,2,1);
h = bar(range_bins,hist_4bin_is(4,:), 1,'stacked');
set(h,'FaceColor',[0.62 0.78 0.02]);
set(gca,'fontsize',16); xlim([0 1.05]);
if length(range_bins)<=50; ylim([0 y_lim_4]); end;
xlabel('bpp(i,j)'); ylabel('#');
legend('Confirmed', 'location','best');
title([caption, ': 4bin in silico'],'fontsize',20,'fontweight','bold');
subplot(4,2,3);
h = bar(range_bins,hist_4bin_is(1,:), 1,'stacked');
set(h,'FaceColor',[1 0.32 0.17]);
set(gca,'fontsize',16); xlim([0 1.05]);
if length(range_bins)<=50; ylim([0 y_lim_4]); end;
legend('Not Confirmed', 'location','best');
xlabel('bpp(i,j)'); ylabel('#');
subplot(4,2,5);
h = bar(range_bins,hist_4bin_is(2,:), 1,'stacked');
set(h,'FaceColor',[0.73 0.72 0.77]);
set(gca,'fontsize',16); xlim([0 1.05]);
if length(range_bins)<=50; ylim([0 y_lim_4]); end;
legend('Uncertain', 'location','best');
xlabel('bpp(i,j)'); ylabel('#');
subplot(4,2,7);
h = bar(range_bins,hist_4bin_is(3,:)', 1,'stacked');
set(h,'FaceColor',[1 0.76 0.03]);
set(gca,'fontsize',16); xlim([0 1.05]);
if length(range_bins)<=50; ylim([0 y_lim_4]); end;
legend('Partial', 'location','best');
xlabel('bpp(i,j)'); ylabel('#');
subplot(4,2,2);
h = bar(range_bins,hist_3bin_is(3,:), 1,'stacked');
set(h,'FaceColor',[0.62 0.78 0.02]);
set(gca,'fontsize',16); xlim([0 1.05]);
if length(range_bins)<=50; ylim([0 y_lim_3]); end;
xlabel('bpp(i,j)'); ylabel('#');
legend('Confirmed', 'location','best');
title([caption, ': 3bin in silico'],'fontsize',20,'fontweight','bold');
subplot(4,2,4);
h = bar(range_bins,hist_3bin_is(1,:), 1,'stacked');
set(h,'FaceColor',[1 0.32 0.17]);
set(gca,'fontsize',16); xlim([0 1.05]);
if length(range_bins)<=50; ylim([0 y_lim_3]); end;
legend('Not Confirmed', 'location','best');
xlabel('bpp(i,j)'); ylabel('#');
subplot(4,2,6);
h = bar(range_bins,hist_3bin_is(2,:), 1,'stacked');
set(h,'FaceColor',[0.73 0.72 0.77]);
set(gca,'fontsize',16); xlim([0 1.05]);
if length(range_bins)<=50; ylim([0 y_lim_3]); end;
legend('Uncertain', 'location','best');
xlabel('bpp(i,j)'); ylabel('#');
