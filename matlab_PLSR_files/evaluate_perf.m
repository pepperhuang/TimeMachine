function rsquared = evaluate_perf(actual, predicted, debug)

if nargin < 3,
    debug=1;
end;


D1=GetDifferenceAngle(actual,repmat(GetMeanAngle(actual),1,length(actual)));
TSS = sum((D1).^2);

D2=GetDifferenceAngle(actual,predicted);
RSS_PLS=sum((D2).^2);

rsquared = 1 - RSS_PLS/TSS;
rsquared

if debug,
    set(gca,'fontsize',16);
    plot (actual, predicted,'+');
    xlabel('Observed Response');
    ylabel('Predicted Response');
    hold on;
    axis_= axis;
    plot([-50 400],[-50 400],'k--');
    title(sprintf('R-squared: %1.3f',rsquared));
end;



