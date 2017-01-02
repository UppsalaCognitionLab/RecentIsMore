% function [fitinfo, plotinfo] = fit_model_exp1(subjidx,condidx,makeplot)
% 
% Fits an ideal observer model to data from experiment 1. See paper for
% details about the experiment.
%
% INPUT
%  subjidx : subject index (integer in  range 1-30)
%  condidx : condition index (1=larger/smaller task, 2=same/different task)
%  makeplot : 0=no plot, 1=produce plot of maximum-likelihood fit 
%
% OUTPUT 
%  fitinfo  : structure with parameter estimates and AIC of the fit
%  plotinfo : data used for making subject-averaged group plot 
%
% This file is part of the code published with the paper "Recent is more: 
% a negative time-order effect in non-symbolic numerical judgment" by 
% R. van den Berg, M. Lindskog, L. Poom, and A. Winman (JEP:HPP, 2017).
%
% For questions, bug reports, etc, please email nronaldvdberg@gmail.com

function [fitinfo, plotinfo] = fit_model_exp1(subjidx,condidx,makeplot)

if ~exist('makeplot','var')
    makeplot=1;
end

% get data
data = read_data(1,subjidx,condidx);

% fit model
cnames = {'larger/smaller','same/different'};
fprintf('Fitting data of subject %d (experiment 1, %s task)...',subjidx,cnames{condidx}); tic
[fitpars, mLLH] = fminsearch(@(pars) -LLH_fun(pars, data, condidx), get_initpars(data,condidx));
fprintf(' (took %2.1f seconds)\n',toc);
LLH = -mLLH;
AIC = -2*LLH+2*numel(fitpars);

% plot
[X_emp, X_fit, Y_emp, Y_fit]=plot_fit(fitpars,data,condidx,makeplot);
plotinfo.X_emp = X_emp;
plotinfo.X_fit = X_fit;
plotinfo.Y_emp = Y_emp;
plotinfo.Y_fit = Y_fit;
fitinfo.LLH = LLH;
fitinfo.fitpars = fitpars;
fitinfo.AIC = AIC;

%--------------------- HELPER FUNCTIONS -------------------------%

function LLH = LLH_fun(pars, data, condidx)
% condidx: 1=LS, 2=SD
sigma=pars(1);
bias_blue = pars(2);
alpha=pars(3);
beta=pars(4);
if condidx==2
    crit=pars(5);
end
if sigma<0
    LLH=-Inf;
    return
end

bias_2nd = alpha + beta*data.Nbar;
if condidx==1
    % larger/smaller task
    pblue1 = 1-normcdf(0,log(data.Nb)+bias_blue-bias_2nd-log(data.Ny),sqrt(2*sigma^2)); % probability of "blue" response when "blue" presented first
    pblue2 = 1-normcdf(0,log(data.Nb)+bias_blue+bias_2nd-log(data.Ny),sqrt(2*sigma^2)); % probaiblity of "blue" response when "blue" presented last
    p_resp(data.C_hat==1 & data.order==1) = pblue1(data.C_hat==1 & data.order==1); % "blue" response, "blue" presented first
    p_resp(data.C_hat==1 & data.order==2) = pblue2(data.C_hat==1 & data.order==2); % "blue" response, "blue" presented last
    p_resp(data.C_hat==2 & data.order==1) = 1-pblue1(data.C_hat==2 & data.order==1); % "yellow" response, "blue" presented first
    p_resp(data.C_hat==2 & data.order==2) = 1-pblue2(data.C_hat==2 & data.order==2); % "yellow" response, "blue" presented last
elseif condidx==2
    % same/different task
    psame1 = normcdf(crit,log(data.Nb)+bias_blue-bias_2nd-log(data.Ny),sqrt(2*sigma^2))-normcdf(-crit,log(data.Nb)+bias_blue-bias_2nd-log(data.Ny),sqrt(2*sigma^2)); % probability of "same" response when "blue" presented first
    psame2 = normcdf(crit,log(data.Nb)+bias_blue+bias_2nd-log(data.Ny),sqrt(2*sigma^2))-normcdf(-crit,log(data.Nb)+bias_blue+bias_2nd-log(data.Ny),sqrt(2*sigma^2)); % probability of "same" response when "blue" presented first
    p_resp(data.C_hat==1 & data.order==1) = psame1(data.C_hat==1 & data.order==1); % "same" response, "blue" presented first
    p_resp(data.C_hat==1 & data.order==2) = psame2(data.C_hat==1 & data.order==2); % "same" response, "blue" presented last
    p_resp(data.C_hat==0 & data.order==1) = 1-psame1(data.C_hat==0 & data.order==1); % "diff" response, "blue" presented first
    p_resp(data.C_hat==0 & data.order==2) = 1-psame2(data.C_hat==0 & data.order==2); % "diff" response, "blue" presented last
end
LLH=sum(log(max(p_resp,1e-3)));

function initpars = get_initpars(data,condidx)
nsteps=8;
sigma_vec=linspace(.01,.5,nsteps);
bias_blue_vec=linspace(-.2,.2,nsteps);
alpha_vec = linspace(-1,1,nsteps);
beta_vec= linspace(-.5,.5,nsteps);
if condidx==1
    crit_vec=0;
else
    crit_vec=linspace(0,.5,nsteps);
end
maxLLH=-Inf;
for ii=1:numel(sigma_vec)
    for jj=1:numel(bias_blue_vec)
        for kk=1:numel(alpha_vec)
            for ll=1:numel(beta_vec)
                for nn=1:numel(crit_vec)
                    pars = [sigma_vec(ii) bias_blue_vec(jj) alpha_vec(kk) beta_vec(ll) crit_vec(nn)];
                    LLH=LLH_fun(pars, data, condidx);
                    if LLH>maxLLH
                        maxLLH=LLH;
                        initpars=pars;
                    end
                end
            end
        end
    end
end
if condidx==1
    initpars=initpars(1:4);
end

function [X_emp, X_fit, Y_emp, Y_fit] = plot_fit(pars,data,condidx,makeplot)
% condidx: 1=LS, 2=SD
uRatio = unique(data.ratio);
eb_blue = [.8 .8 .8];
eb_red = [1 .8 .8];
sigma=pars(1);
bias_blue = pars(2);
alpha=pars(3);
beta=pars(4);
X_emp = uRatio;
X_fit = uRatio;
if makeplot
    figure
end
if condidx==1
    % empirical
    for ii=1:numel(uRatio)
        Y_emp(1,ii) = mean(data.C_hat(data.ratio==uRatio(ii) & data.order==1)==1); % probability of "blue" response when blue presented first
        Y_emp(2,ii) = mean(data.C_hat(data.ratio==uRatio(ii) & data.order==2)==1); % probability of "blue" response when blue presented first
        Y_emp_eb(1,ii) = 0;
        Y_emp_eb(2,ii) = 0;
    end
    % fit 
    bias_2nd = alpha + beta*data.Nbar;
    pblue1 = 1-normcdf(0,log(data.Nb)+bias_blue-bias_2nd-log(data.Ny),sqrt(2*sigma^2)); % probability of "blue" response when "blue" presented first
    pblue2 = 1-normcdf(0,log(data.Nb)+bias_blue+bias_2nd-log(data.Ny),sqrt(2*sigma^2)); % probaiblity of "blue" response when "blue" presented last
    p_blue_resp(data.order==1) = pblue1(data.order==1); % "blue" response, "blue" presented first
    p_blue_resp(data.order==2) = pblue2(data.order==2); % "blue" response, "blue" presented last
    for ii=1:numel(uRatio)
        Y_fit(1,ii) = mean(p_blue_resp(data.ratio==uRatio(ii) & data.order==1)); % probability of "blue" response when blue presented first
        Y_fit(2,ii) = mean(p_blue_resp(data.ratio==uRatio(ii) & data.order==2)); % probability of "blue" response when blue presented first
        Y_fit_eb(1,ii) = 0;
        Y_fit_eb(2,ii) = 0;
    end
    if makeplot
        errorbar(log(uRatio),Y_emp(1,:),Y_emp_eb(1,:),'ro','markerfacecolor','r','color',eb_red);
        hold on
        errorbar(log(uRatio),Y_emp(2,:),Y_emp_eb(2,:),'bo','markerfacecolor','b','color',eb_blue);
        plot(log(uRatio),Y_emp(1,:),'ro','markerfacecolor','r');
        plot(log(uRatio),Y_emp(2,:),'bo','markerfacecolor','b');
        plot(log(uRatio),Y_fit(1,:),'r-');
        plot(log(uRatio),Y_fit(2,:),'b-');
        ylabel('Proportion "blue" responses');
    end
elseif condidx==2
    crit=pars(5);
    % empirical
    for ii=1:numel(uRatio)
        Y_emp(1,ii) = mean(data.C_hat(data.ratio==uRatio(ii) & data.order==1)==1); % probability of "same" response when blue presented first
        Y_emp(2,ii) = mean(data.C_hat(data.ratio==uRatio(ii) & data.order==2)==1); % probability of "same" response when blue presented first
        Y_emp_eb(1,ii) = 0;
        Y_emp_eb(2,ii) = 0;
    end
    % fit
    bias_2nd = alpha + beta*data.Nbar;
    psame1 = normcdf(crit,log(data.Nb)+bias_blue-bias_2nd-log(data.Ny),sqrt(2*sigma^2))-normcdf(-crit,log(data.Nb)+bias_blue-bias_2nd-log(data.Ny),sqrt(2*sigma^2)); % probability of "same" response when "blue" presented first
    psame2 = normcdf(crit,log(data.Nb)+bias_blue+bias_2nd-log(data.Ny),sqrt(2*sigma^2))-normcdf(-crit,log(data.Nb)+bias_blue+bias_2nd-log(data.Ny),sqrt(2*sigma^2)); % probability of "same" response when "blue" presented first
    p_same_resp(data.order==1) = psame1(data.order==1); % "same" response, "blue" presented first
    p_same_resp(data.order==2) = psame2(data.order==2); % "same" response, "blue" presented last    
    for ii=1:numel(uRatio)
        Y_fit(1,ii) = mean(p_same_resp(data.ratio==uRatio(ii) & data.order==1)); % probability of "same" response when blue presented first
        Y_fit(2,ii) = mean(p_same_resp(data.ratio==uRatio(ii) & data.order==2)); % probability of "same" response when blue presented first
        Y_fit_eb(1,ii) = 0;
        Y_fit_eb(2,ii) = 0;
    end
    if makeplot
        errorbar(log(uRatio),Y_emp(1,:),Y_emp_eb(1,:),'ro','markerfacecolor','r','color',eb_red);
        hold on
        errorbar(log(uRatio),Y_emp(2,:),Y_emp_eb(2,:),'bo','markerfacecolor','b','color',eb_blue);
        plot(log(uRatio),Y_emp(1,:),'ro','markerfacecolor','r');
        plot(log(uRatio),Y_emp(2,:),'bo','markerfacecolor','b');
        plot(log(uRatio),Y_fit(1,:),'r-');
        plot(log(uRatio),Y_fit(2,:),'b-');
        ylabel('Proportion "same" responses');
    end
end
if makeplot
    l=legend({'blue first','blue last'});
    set(l,'location','northwest','fontsize',8);
    xlim([-.75 .75]);
    ylim([0 1]);
    xlabel('Blue/Yellow ratio (log)');
    grid on
end

