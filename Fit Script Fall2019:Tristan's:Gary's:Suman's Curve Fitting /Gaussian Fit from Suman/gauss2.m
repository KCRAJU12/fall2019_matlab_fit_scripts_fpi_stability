function [fitresult, gof] = gauss2(xdata, ydata, initialguess, T)
[xdata, ydata] = prepareCurveData( xdata, ydata );

% Set up fittype and options.
ft = fittype( 'a1*exp(-0.5*((x- b1)/w1 )^2)+ a2*exp(-0.5*((x- b2)/w2 )^2)','independent','x','coefficients',{'a1','b1','w1','a2','b2','w2' });
%ft = fittype( 'gauss2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Upper = [2.6  15.2  0.09   12   15.8   0.1];
opts.Lower = [0.1  0    0    0.1   0.1     0];
opts.Normalize = 'off'; %% this was causing me the problem always off
opts.StartPoint = initialguess ; %[1  0.2  0.3  4.9  0.7  0.4];
% opts.Weights = weight;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fit model to data.
[fitresult, gof ] = fit( xdata, ydata, ft, opts );
va = abs (coeffvalues(fitresult) );
uncer = diff(confint(fitresult)/2);%uncertanity is given by subtracting max and min and dividing by 2

%%
figure(12)
%figure( 'Name', ' Double Gaussian' );
plot( fitresult,'b.');
 hold on
 plot( xdata, ydata,'b.', 'markers',12) ;
 hold on
errorbar(xdata, ydata, sqrt(ydata)./sum(ydata), -sqrt(ydata)./sum(ydata), 'b.', 'markers', 12)
legend off
 set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.02 .02], ...
                                 'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off','XGrid', 'off', ...
                                 'XColor', [0 0 0], 'YColor', [0 0 0] )
%title([ 'N_{1}=' num2str( round(va(1),3) ) ',' ' ' '\mu_{1}=' num2str(round(va(2),3)) ',' ' ' '\sigma_{1}=' num2str(round(va(3),3)) ',' ' ' 'N_{2}=' num2str(round(va(4),3)) ',' ' '  '\mu_{2}=' num2str(round(va(5),3)) ',' ' ' '\sigma_{2}=' num2str(round(va(6),3))   ])
 xlabel('\bf Current (\muA)','interpreter', 'tex')
 %xlabel(['$(\mathrm{\overline{u}})(m/s)$'],'interpreter','latex')
 ylabel('\bf P','Fontsize',14)
 set(gca, 'fontsize', 14)
 %grid on 
 %$(\mathrm{\overline{u}})(m/s)$'
 %%
fit_str1= ['\bf N_{1}=  ' , num2str( round(va(1),3),'%+.3f'),'\pm',num2str(uncer(1),'%.3f')] ;
fit_str2= ['\bf \mu_{1}=  ' , num2str( round(va(2),3),'%+.3f'),'\pm',num2str(uncer(2),'%.3f')] ; 
fit_str3= ['\bf \sigma_{1}=  ' , num2str( round(va(3),3),'%+.3f'),'\pm',num2str(uncer(3),'%.3f')] ;
fit_str4= ['\bf N_{2}=  ' , num2str( round(va(4),3),'%+.3f'),'\pm',num2str(uncer(4),'%.3f')] ;
fit_str5= ['\bf \mu_{2}=  ' , num2str( round(va(5),3),'%+.3f'),'\pm',num2str(uncer(5),'%.3f')] ;
fit_str6= ['\bf \sigma_{2}=  ' , num2str( round(va(6),3),'%+.3f'),'\pm',num2str(uncer(6),'%.3f') ];
fit_str7= ['\bf T=   ', num2str( round(T,3)),' ' 'K' ];
hleg = legend({'foo', 'bar'}, 'Location', 'best');
axes('position', get(hleg, 'position'));
text(0, 0.2,fit_str1,'color',[0 0 1]);
text(0,-0.2,fit_str2,'color',[0 0 1]);
text(0,-0.6,fit_str3,'color',[0 0 1]);
text(0,-1.0,fit_str4,'color',[0 0 1]);
text(0,-1.4,fit_str5,'color',[0 0 1]);
text(0,-1.8,fit_str6,'color',[0 0 1]);
text(0,-2.6,fit_str7,'color',[0 0 1]);
axis off
delete(hleg);

%% chi square test
yfitup = feval(fitresult, xdata); % gives the fitted y value
for kk=1:numel(ydata)
chisquaree(kk)=((ydata(kk) - yfitup(kk)))^2/yfitup(kk);
end
chisquaree(chisquaree == inf)= 0;
chi2up = sum(chisquaree);
%%
nn = length(xdata);                % number of data points
m = 6;                                 % number of fitting parameters a1, a2, ...
reduce_chi2up = chi2up/(nn-m);
disp(reduce_chi2up);
fit_st=['\bf{\bf \chi^2}_{reduced} =  ', num2str( reduce_chi2up,'%+.3f')];
text(0,-2.2,fit_st,'color','b');
%% Goodness-of-Fit
tgof = '   *** Acceptable Fit   ***';
if reduce_chi2up > 1, tgof = '   ??? Fit may not be acceptable ???'; end
if reduce_chi2up <= 1, tgof = '   ? Fit may be too good ?'; end
disp(tgof);


 
 
 
 
%[fitresult, gof ] = fit( xdata', ydata', ft, opts );
%% threshold voltage and finding root of quadratic equation
%va = coeffvalues(fitresult);
%aa =  ( va(6)^2- va(3)^2  );
%bb =  2*( -va(2) * va(6)^2 + va(5) * va(3)^2 ) ;
%cc = - va(5)^2 * va(3)^2 + va(2)^2 * va(6)^2 - 2* va(3)^2 * va(6)^2 * log(va(1)./va(4));
%p = [aa bb cc];
%r = roots(p) ;
%vthh = r(1) ;
%vth = r(2) ;


%%
% set(0,'DefaultFigureVisible','on')% turn  on 
% %%
% % Plot fit with data.
% figure( 'Name', ' Double Gaussian' );
% plot( fitresult,'r.');
% hold on
% plot( xData, yData,'r.', 'markers',12) ;
% %legend( h, 'ydata vs. xdata', 'untitled fit 1', 'Location', 'NorthEast' );
% 
% % Label axes
% xlabel('\bf V_{out}','Fontsize',12);
% ylabel('\bf PDF','Fontsize',12);
% get(gca, 'XTick'); 
% set(gca, 'Fontsize',14);
% grid on
% plotbrowser('off')

% %% another method for fitting matlab
% fun = @(x,xdata)x(1)*exp (( xdata - x(2))./(sqrt(2)*x(3))).^2 + x(4)*exp (( xdata - x(5))./(sqrt(2)*x(6))).^2;
% x0 = [2,0.5, 0.02, 2, 0.5, 0.02];
% options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective');
% lb = [3,-0.1, 0, 2, 0.1, 0 ];
% ub = [6, 0.8, 1, 8, 0.9, 1];
% x = lsqcurvefit(fun,x0,xdata,ydata); %, lb,ub,options) ;
% %x = lsqnonlin(fun,x0,xdata,ydata,lb,ub,options);

