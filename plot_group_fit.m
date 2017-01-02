function plot_group_fit(expnr)

subjlist = get_included_subjects(expnr);

% create group plot
shade_blue = [.8 .8 1];
shade_red   = [1 .8 .8];
figure
set(gcf,'Position',get(gcf,'Position').*[0.1 0.1 2.4 .8]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[0.1 0.1 2.4 .8]);
ncondidx_vec = [2 3];
for condidx=1:ncondidx_vec(expnr)
    clear all_Y_emp1 all_Y_emp2 all_Y_fit1 all_Y_fit2;
    % get plot for each subject
    cnt=0;
    for ii=1:numel(subjlist)
        subjidx = subjlist(ii);
        if expnr==1
            [fitinfo, plotinfo] = fit_model_exp1(subjidx,condidx,0);
        elseif expnr==2
            [fitinfo, plotinfo] = fit_model_exp2(subjidx,condidx,0);
        end
        cnt=cnt+1;
        
        all_Y_emp1(cnt,:) = plotinfo.Y_emp(1,:);
        all_Y_emp2(cnt,:) = plotinfo.Y_emp(2,:);
        
        all_Y_fit1(cnt,:) = plotinfo.Y_fit(1,:);
        all_Y_fit2(cnt,:) = plotinfo.Y_fit(2,:);
    end
    
    uRatio_emp = plotinfo.X_emp';
    uRatio_fit = plotinfo.X_fit';
    
    nsubj = size(all_Y_emp1,1);
    
    % compute means and std    
    Y1_emp_mean = nanmean(all_Y_emp1);
    Y2_emp_mean = nanmean(all_Y_emp2);
    Y1_emp_eb = nanstd(all_Y_emp1)/sqrt(nsubj);
    Y2_emp_eb = nanstd(all_Y_emp2)/sqrt(nsubj);
    Y1_fit_mean = nanmean(all_Y_fit1);
    Y2_fit_mean = nanmean(all_Y_fit2);
    Y1_fit_eb = nanstd(all_Y_fit1)/sqrt(nsubj);
    Y2_fit_eb = nanstd(all_Y_fit2)/sqrt(nsubj);
    
    % plot
    subplot(1,3,condidx)
    hold on
    XX=[log(uRatio_fit) log(uRatio_fit(end:-1:1))];
    YY=[Y1_fit_mean+Y1_fit_eb Y1_fit_mean(end:-1:1)-Y1_fit_eb(end:-1:1)];
    fill(XX,YY,shade_red,'linestyle','none');
    plot(log(uRatio_fit),Y1_fit_mean,'r','color',[1 0 0]);
    YY=[Y2_fit_mean+Y2_fit_eb Y2_fit_mean(end:-1:1)-Y2_fit_eb(end:-1:1)];
    fill(XX,YY,shade_blue,'linestyle','none');
    plot(log(uRatio_fit),Y2_fit_mean,'b','color',[0 0 1]);
    q(1)=errorbar(log(uRatio_emp),Y1_emp_mean,Y1_emp_eb,'ro','markerfacecolor','r');
    q(2)=errorbar(log(uRatio_emp),Y2_emp_mean,Y2_emp_eb,'bo','markerfacecolor','b');

    xlabel('Blue/Yellow ratio (log)');
    if expnr==1
        xlim([-.75 .75]);
    elseif expnr==2
        xlim([-.4 .4])
    end
    ylim([0 1]);

    grid on
    if expnr==1 && condidx==2
        ylabel('Proportion "same" responses');
    else
        ylabel('Proportion "blue" responses');
    end
    drawnow;
end
