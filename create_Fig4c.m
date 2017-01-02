% function create_Fig4c
% 
% Reproduces Fig 4c of the paper "Recent is more: a negative time-order effect 
% in non-symbolic numerical judgment" (Van den Berg et al. 2017, JEP:HPP)
%   
% For questions, bug reports, etc, please email nronaldvdberg@gmail.com

function create_Fig4c

% get list of included subjects
subjlist = get_included_subjects(3);

% compute mean and std dev of over estimation for each of the three
% presented numerosities
uN = [8 11 14];
for ii=1:numel(subjlist)
    data = read_data(3,subjlist(ii),1);
    for jj=1:numel(uN)
        o0 = 100*(data.N1_hat(data.N1==uN(jj) & data.single)./uN(jj)-1);
        o1 = 100*(data.N1_hat(data.N1==uN(jj) & ~data.single)./uN(jj)-1);
        o2 = 100*(data.N2_hat(data.N2==uN(jj) & ~data.single)./uN(jj)-1);
        mean_o0(ii,jj) = mean(o0);
        mean_o1(ii,jj) = mean(o1);
        mean_o2(ii,jj) = mean(o2);
        mean_std0(ii,jj) = std(data.N1_hat(data.N1==uN(jj) & data.single));
        mean_std1(ii,jj) = std(data.N1_hat(data.N1==uN(jj) & ~data.single));
        mean_std2(ii,jj) = std(data.N2_hat(data.N2==uN(jj) & ~data.single));
    end
end

% plot
figure
set(gcf,'Position',get(gcf,'Position').*[0.1 0.1 1.2 .7]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[0.1 0.1 1.2 .7]);
subplot(1,2,1);
hold on
Y_mean0 = mean(mean_o0);
Y_std0 = std(mean_o0)/sqrt(numel(subjlist));
Y_mean1 = mean(mean_o1);
Y_std1 = std(mean_o1)/sqrt(numel(subjlist));
Y_mean2 = mean(mean_o2);
Y_std2 = std(mean_o2)/sqrt(numel(subjlist));
% plot([6 16],[6 16],'k')
errorbar(uN-0,Y_mean0,Y_std0,'k--','color',[.7 .7 .7]);
errorbar(uN-0.2,Y_mean1,Y_std1,'ks-','markerfacecolor','k');
errorbar(uN+0.2,Y_mean2,Y_std2,'ro-','markerfacecolor','r');
ylim([0 35])
set(gca,'Xtick',uN);
ylabel('mean(overestimation)');

subplot(1,2,2);
hold on
Y_mean0 = mean(mean_std0);
Y_std0 = std(mean_std0)/sqrt(numel(subjlist));
Y_mean1 = mean(mean_std1);
Y_std1 = std(mean_std1)/sqrt(numel(subjlist));
Y_mean2 = mean(mean_std2);
Y_std2 = std(mean_std2)/sqrt(numel(subjlist));
errorbar(uN-0,Y_mean0,Y_std0,'k--','color',[.7 .7 .7]);
errorbar(uN-0.2,Y_mean1,Y_std1,'ks-','markerfacecolor','k');
errorbar(uN+0.2,Y_mean2,Y_std2,'ro-','markerfacecolor','r');
% ylim([10 30])
set(gca,'Xtick',uN);
ylabel('std(overestimation)');
ylim([0 4.5])
