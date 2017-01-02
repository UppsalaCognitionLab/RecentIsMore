function create_Fig2c

condidx=1:3;
subjlist = get_included_subjects(2);

figure
set(gcf,'Position',get(gcf,'Position').*[0.1 0.1 2 .8]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[0.1 0.1 2 .8]);
ISI = [50 300 2000];
X = linspace(1.7,2.71,25);
for ii=1:numel(condidx)
    subplot(1,3,ii);
    hold on
    ylim([-30 15]);
    xlim([X(1)-.02 X(end)+.02]);
    set(gca,'Xtick',[1.70 2.0127 2.3937 2.7058],'Xticklabel',[11 15 22 30]);
    xlabel('Stimulus magnitude');
    ylabel('TOE (in percentage)');
    title(sprintf('ISI=%dms',ISI(ii)));
    for jj=1:numel(subjlist)
        fitinfo = fit_model_exp2(subjlist(jj),condidx(ii),0);
        Y(jj,:) = fitinfo.fitpars(3) + fitinfo.fitpars(4)*X;
        Y(jj,:) = -100*(exp(Y(jj,:))-1); % bias in terms of percentage overestimation
        plot(X,Y(jj,:),'k-','color',[.75 .75 .75]);
        beta(jj) = fitinfo.fitpars(4);
        drawnow;
    end
    plot(X,mean(Y),'k-','Linewidth',2);
end