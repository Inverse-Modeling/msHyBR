close all

exp_name = 'deblurring_example';

c1 = [0 0.4470 0.7410];
c2 = [0.8500 0.3250 0.0980];
c3 = [0.4660 0.6740 0.1880];
c4 = [0.4940 0.1840 0.5560];
c5 = [0.3010 0.7450 0.9330];
c6 = [0.6350 0.0780 0.1840];
c7 = [0.9290 0.6940 0.1250];

lw=4;
fs=20;

%%
% Plot basis
figure
plot(domain,X,'LineWidth',lw)
set(findall(gcf,'-property','FontSize'),'FontSize',fs)
title('\bf{basis vectors}','interpreter','latex','FontSize',fs)
saveas(gcf,['Figures/', exp_name,'_basis'],'epsc')

% Plot measurements 
figure
hold on
plot(domain, b,'-','LineWidth',lw)
plot(domain, bn,'-.','LineWidth',lw)
xl=legend('blurred', 'blurred and noisy');
set(findall(gcf,'-property','FontSize'),'FontSize',fs)
xl.Position(1) = xl.Position(1) - 0.01;
title('\bf{measurements}','interpreter','latex','FontSize',fs)
saveas(gcf,['Figures/',exp_name,'_measurements'],'epsc')

%% Plot reconstructions 
% Plot reconstructions - iterative methods
figure
plot(domain,x_true,'-','LineWidth',lw,'Color','k')
hold on
plot(domain,x_ss, '--','LineWidth',lw,'Color',c2)
axis([0 1 -0.15 2.5])
xl=legend('true', 'ssHyBR');
set(findall(gcf,'-property','FontSize'),'FontSize',fs)
xl.Position(1) = xl.Position(1) - 0.04;
saveas(gcf,['Figures/',exp_name,'_reconstructions'],'epsc')

% Plot reconstructions - two step methods - forward
figure
plot(domain,x_true,'-','LineWidth',lw,'Color','k')
hold on
plot(domain,x_ss, '--','LineWidth',lw,'Color',c2)
plot(domain,s_forward_VT, '-','LineWidth',lw,'Color',c1)
plot(domain,s_forward_BIC, '-.','LineWidth',lw,'Color',c3)
plot(domain,s_forward_AIC, ':','LineWidth',lw,'Color', c7)
title('solutions of forward selection methods')
axis([0 1 -0.15 2.5])
xl=legend('true', 'ssHyBR', 'forward VT','forward AIC','forward BIC');
set(findall(gcf,'-property','FontSize'),'FontSize',fs)
xl.Position(1) = xl.Position(1) - 0.06;
saveas(gcf,['Figures/',exp_name,'_reconstructions2'],'epsc')

% Plot reconstructions - two step methods - exhaustive
figure
plot(domain,x_true,'-','LineWidth',lw,'Color','k')
hold on
plot(domain,x_ss, '--','LineWidth',lw,'Color',c2)
plot(domain,s_exhaustive_BIC, '-.','LineWidth',lw,'Color',c3)
plot(domain,s_exhaustive_AIC, ':','LineWidth',lw,'Color', c7)
axis([0 1 -0.15 2.5])
title('solutions of exhaustive selection methods')
xl=legend('true', 'ssHyBR', 'exhaustive AIC','exhaustive BIC');
set(findall(gcf,'-property','FontSize'),'FontSize',fs)
xl.Position(1) = xl.Position(1) - 0.06;
saveas(gcf,['Figures/',exp_name,'_reconstructions3'],'epsc')

%% Plot betas 

% Plot beta iterative methods
figure
plot(beta_true,'-','LineWidth',lw,'Color','k')
hold on
plot(beta_ss,'--','LineWidth',lw,'Color',c2);
xl= legend('true', 'ssHyBR');
title('beta for forward methods')
axis([1 7 -0.4 1.45])
set(findall(gcf,'-property','FontSize'),'FontSize',fs)
xl.Position(1) = xl.Position(1) - 0.44;
xl.Position(2) = xl.Position(2) - 0.04;
saveas(gcf,['Figures/',exp_name,'_betas'],'epsc')

% Plot beta reconstructions - two step methods - forward
figure
plot(beta_true,'-','LineWidth',lw,'Color','k')
hold on
plot(beta_ss,'--','LineWidth',lw,'Color',c2)
plot(beta_forward_VT_for_comparison,'-','LineWidth',lw, 'Color',c1)
plot(beta_forward_BIC_for_comparison,'-.','LineWidth',lw, 'Color',c3)
plot(beta_forward_AIC_for_comparison,':','LineWidth',lw, 'Color',c7)
xl=legend('true', 'ssHyBR', 'forward VT','forward AIC','forward BIC');
title('beta for forward methods')
axis([1 7 -0.4 1.45])
set(findall(gcf,'-property','FontSize'),'FontSize',fs)
xl.Position(1) = xl.Position(1) - 0.43;
saveas(gcf,['Figures/',exp_name,'_betas2'],'epsc')

% Plot beta reconstructions - two step methods - forward
figure
plot(beta_true,'-','LineWidth',lw,'Color','k')
hold on
plot(beta_ss,'--','LineWidth',lw,'Color',c2)
plot(beta_exhaustive_BIC_for_comparison,'-.','LineWidth',lw,'Color',c3)
plot(beta_exhaustive_AIC_for_comparison,':','LineWidth',lw,'Color',c7)
title('beta for exhaustive methods')
axis([1 7 -0.4 1.45])
xl=legend('true', 'ssHyBR','exhaustive AIC','exhaustive BIC');
set(findall(gcf,'-property','FontSize'),'FontSize',fs)
xl.Position(1) = xl.Position(1) - 0.40;
xl.Position(2) = xl.Position(2) + 0.01;
saveas(gcf,['Figures/',exp_name,'_betas3'],'epsc')


%% Plot norms 
figure
hold on
plot(output_ss.Enrm, '--','LineWidth',lw,'Color',c2)
legend('\bf ssHyBR','interpreter','latex','FontSize',fs);
xlabel('{\bf iterations}','interpreter','latex','FontSize',fs) 
ylabel('{\bf error norm}','interpreter','latex','FontSize',fs) 
set(findall(gcf,'-property','FontSize'),'FontSize',fs)
axis([1 50 0.05 0.2])
saveas(gcf,['Figures/',exp_name,'_Enrm'],'epsc')

%% Compute tables for the evaluation of binary classifiers

thr_beta = 0.2;
beta_ind_msHyBR=zeros(7,1);
beta_ind_msHyBR(abs(beta_ss)>thr_beta ) = 1;

beta_ind_exh_BIC= zeros(7,1);
beta_ind_exh_BIC(info_exhaustive_BIC.indices_selected)=1;
beta_ind_exh_AIC= zeros(7,1);
beta_ind_exh_AIC(info_exhaustive_AIC.indices_selected)=1;

beta_ind_for_BIC= zeros(7,1);
beta_ind_for_BIC(info_forward_BIC.indices_selected)=1;
beta_ind_for_AIC= zeros(7,1);
beta_ind_for_AIC(info_forward_AIC.indices_selected)=1;
beta_ind_for_VT= zeros(7,1);
beta_ind_for_VT(info_forward_VT.indices_selected)=1;

%%
name = {'msHyBR','exh_BIC','exh_AIC','for_BIC','for_AIC','for_VT'};

fprintf('\n <strong> Binary classifier metrics </strong> \n \n');
fprintf('   TP    FP   TN   FN   F1     name \n');
for i = 1:6
   [TPR,TNR,PPV, F1, TP, TN, FP, FN]=binary_class_eval( eval(['beta_ind_',name{i}]), beta_true);
   fprintf(['%5d%5d%5d%5d    %4.2f   ',name{i},'\n'], TP, FP, TN, FN, F1)
end
