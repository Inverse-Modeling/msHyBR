close all
clear all

load('Ex2_results')
title_ex = 'Ex2_1';

c1 = [0 0.4470 0.7410];
c2 = [0.8500 0.3250 0.0980];
c3 = [0.4660 0.6740 0.1880];
c4 = [0.4940 0.1840 0.5560];
c5 = [0.3010 0.7450 0.9330];
c6 = [0.6350 0.0780 0.1840];
c7 = [0.9290 0.6940 0.1250];

lw=4;
fs=22;

%% Plot true solution
% plot mean
visual_no_Africa_no_sea(mean_true(land), land, lon, lat, 1); 
title('\bf true mean ({\boldmath$X \beta$})','interpreter','latex','FontSize',fs)
clim([min(s_true(land)) max(s_true(land))])
cc = colorbar; set(cc,'FontSize',fs); 
set(gca,'FontSize',fs)
saveas(gcf,['Figures/',title_ex,'_truemean'],'epsc')

% plot stochastic component
visual_no_Africa_no_sea(stochastic_true(land), land, lon, lat, 1); 
title('\bf true stochastic component ({\boldmath$\zeta$})','interpreter','latex','FontSize',fs)
set(gca,'FontSize',fs)
clim([min(s_true(land)) max(s_true(land))])
set(gca,'FontSize',fs)
cc = colorbar; set(cc,'FontSize',fs); 
saveas(gcf,['Figures/',title_ex,'_truestochastic'],'epsc')

% plot full true solution
visual_no_Africa_no_sea(s_true(land), land, lon, lat, 1); 
title('\bf{true emissions ({\boldmath$s$})}','interpreter','latex','FontSize',fs)
clim([min(s_true(land)) max(s_true(land))])
set(gca,'FontSize',fs)
cc = colorbar; set(cc,'FontSize',fs); 
saveas(gcf,['Figures/',title_ex,'_true'],'epsc')

%% Plot zonation model basis representation
grid_to_plot_ind = zeros(170,70);
grid_to_plot_ind(1:10:170,:) = 1;
grid_to_plot_ind(:,1:7:70) = 1;
grid_to_plot = ones(size(grid_to_plot_ind));
grid_to_plot(grid_to_plot_ind==1)=10;
visual_no_Africa_no_sea(grid_to_plot(land)/10, land, lon, lat, 1); 
title('\bf{zonation model}','interpreter','latex','FontSize',fs)
set(gca,'FontSize',fs)
saveas(gcf,['Figures/',title_ex,'_basis'],'epsc')

%% Plot error norms
figure
hold on; box off
semilogy(output_mshybr.Enrm,'-','Color', c2,'LineWidth',lw);
semilogy(output_genGKmean.Enrm,'-.','Color', c3,'LineWidth',lw); 
axis([1 50 0 1])
legend('{\bf msHyBR}','{\bf two-step}','interpreter','latex','FontSize',fs)
set(gca,'FontSize',fs)
xlabel('{\bf iterations}','interpreter','latex','FontSize',fs) 
ylabel('{\bf error norm}','interpreter','latex','FontSize',fs) 
saveas(gcf,['Figures/',title_ex,'_errors'],'epsc')

%% Plot betas
beta_twostep_full=zeros(size(beta_true));
beta_twostep_full(info_selected.indices_selected)=beta_twostep;
figure
plot(beta_true,'--','Color', [0 0 0],'LineWidth',3)
hold on
plot(beta_mshybr,'-','Color', c2,'LineWidth',lw)
plot(beta_twostep_full,'-.','Color', c3,'LineWidth',lw); 
axis([0 size(beta_true,1) -0.4 1.6])
set(gca,'FontSize',fs)
legend('{\bf true}','{\bf msHyBR}','\bf two-step','interpreter','latex','FontSize',fs)
xlabel('{\bf index}','interpreter','latex','FontSize',fs) 
ylabel('{\bf beta values {\boldmath${\beta}$} }','interpreter','latex','FontSize',fs) 
saveas(gcf,['Figures/',title_ex,'_estimated_betas'],'epsc')


%% Plot reconstructions using msHyBR
Xused_new_ind = zeros(size(beta_true));
Xused_new_ind(abs(beta_mshybr)>0.20)=1;
Xselected_mshybr = X*Xused_new_ind;
mean_mshybr = X*beta_mshybr;
stoc_mshybr = x_mshybr-X*beta_mshybr;

% plot subset selection
visual_no_Africa_no_sea(Xselected_mshybr(land), land, lon, lat, 1); 
title('\bf{basis vectors representation}','interpreter','latex','FontSize',fs)
set(gca,'XTick',[], 'YTick', []); %cc = colorbar; set(cc,'FontSize',fs); 
saveas(gcf,['Figures/',title_ex,'1'],'epsc')

% plot full true solution
visual_no_Africa_no_sea(x_mshybr(land), land, lon, lat, 1); 
title('\bf Reconstructions of {\boldmath${s}$}','interpreter','latex','FontSize',fs)
set(gca,'XTick',[], 'YTick', []); cc = colorbar; set(cc,'FontSize',fs); clim([min(s_true(land)) max(s_true(land))])
saveas(gcf,['Figures/',title_ex,'2'],'epsc')

% plot stochastic component
visual_no_Africa_no_sea(mean_mshybr(land), land, lon, lat, 1); 
title('\bf mean ({\boldmath$X \beta$})','interpreter','latex','FontSize',fs)
set(gca,'XTick',[], 'YTick', []); cc = colorbar; set(cc,'FontSize',fs); clim([min(s_true(land)) max(s_true(land))])
saveas(gcf,['Figures/',title_ex,'3'],'epsc')

% plot mean 
visual_no_Africa_no_sea(stoc_mshybr(land), land, lon, lat, 1); 
title('\bf stochastic component ({\boldmath$\zeta$})','interpreter','latex','FontSize',fs)
set(gca,'XTick',[], 'YTick', []);  cc = colorbar; set(cc,'FontSize',fs); 
clim([min(s_true(land)) max(s_true(land))])
saveas(gcf,['Figures/',title_ex,'4'],'epsc')


%% Plot reconstructions using a two-step process (forward BIC + genHyBRmean)
Xused_twostep = X(:,info_selected.indices_selected)*ones(size(info_selected.indices_selected,2),1);
mean_twostep = X(:,info_selected.indices_selected)*beta_twostep;
stoc_twostep = x_twostep-X(:,info_selected.indices_selected)*beta_twostep;

% plot subset selection
visual_no_Africa_no_sea(Xused_twostep(land), land, lon, lat, 1); 
title('\bf{basis vectors representation}','interpreter','latex','FontSize',fs)
set(gca,'XTick',[], 'YTick', []);
saveas(gcf,['Figures/',title_ex,'5'],'epsc')

% plot subset selection
visual_no_Africa_no_sea(x_twostep(land), land, lon, lat, 1); 
title('\bf Reconstructions of {\boldmath${s}$}','interpreter','latex','FontSize',fs)
set(gca,'XTick',[], 'YTick', []); cc = colorbar; set(cc,'FontSize',fs); clim([min(s_true(land)) max(s_true(land))])
saveas(gcf,['Figures/',title_ex,'6'],'epsc')

% plot mean
visual_no_Africa_no_sea(mean_twostep(land), land, lon, lat, 1); 
set(gca,'XTick',[], 'YTick', []); cc = colorbar; set(cc,'FontSize',fs); clim([min(s_true(land)) max(s_true(land))])
title('\bf mean ({\boldmath$X \beta$})','interpreter','latex','FontSize',fs)
saveas(gcf,['Figures/',title_ex,'7'],'epsc')

% plot stochastic component
visual_no_Africa_no_sea(stoc_twostep(land), land, lon, lat, 1); 
title('\bf stochastic component ({\boldmath$\zeta$})','interpreter','latex','FontSize',fs)
set(gca,'XTick',[], 'YTick', []); cc = colorbar; set(cc,'FontSize',fs); clim([min(s_true(land)) max(s_true(land))])
saveas(gcf,['Figures/',title_ex,'8'],'epsc')

