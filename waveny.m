dt = 1;
load 'NY.dat';
s0 = NY(40:127);
t = [0:length(s0)-1];
variance = std(s0)^2;
s = (s0 - mean(s0))/sqrt(variance);

dj = 0.01;
scales = [0:dj:4.5];
scales = 2.^scales;
n = length(s);
xlim = [0,(n-1)];
[wave,period,coi]=wavelet2(s,1, scales);
power = abs(wave).^2;

global_ws = (sum(power')/n);

Cdelta = 0.776;
scale_avg = (scales')*(ones(1,n));
scale_avg = power ./ scale_avg;
scale_avg = dj*dt/Cdelta*sum(scale_avg);

%Significance test
nsample = 100;
sig1 = zeros(length(scales), nsample);
sig2 = zeros(nsample, length(s));
for i=1:nsample
    st = datasample(s, length(s));
    [wt,pt,ct] = wavelet2(st,dt,scales);
    powert = (abs(wt)).^2;    
    sig1(:,i)=(sum(powert')/n);
    temp = (scales')*(ones(1,n));
    temp = powert ./ temp;
    temp = dj*dt/Cdelta*sum(temp);
    sig2(i,:) = temp;
end
sig1_l = quantile(sig1, 0.05, 2);
sig1_u = quantile(sig1, 0.95, 2);
sig2_l = quantile(sig2,0.05,1);
sig2_u = quantile(sig2,0.95,1);

subplot('position',[0.1 0.75 0.65 0.2])
%subplot(2,1,1);
plot(t,s0);
set(gca,'XLim',xlim(:))
xlabel('Time (day)')
ylabel('Count')
xline(21,'--');%NY 21, CA 18, WA 22
title('NY COVID-19 Daily New Cases')

Yticks = 2.^(fix(log2(min(period))):fix(log2(max(period))));
subplot('position',[0.1 0.37 0.65 0.28])
%subplot(2,1,2);
imagesc(t,log2(period),power);
%contour(t,log2(period),power); 
colorbar;
caxis([0 3.0]);
xlabel('Time (day)')
ylabel('Log_2(Period (day))')
title('b) Wavelet Power Spectrum')
%set(gca,'XLim',xlim(:))
%set(gca,'YLim',log2([min(period),max(period)]), ...
%	'YDir','reverse', ...
%	'YTick',log2(Yticks(:)), ...
%	'YTickLabel',Yticks)
hold on
plot(t,log2(coi),'w--')

%[signif,fft_theor] = wave_signif(1.0,dt,scales,0,0.72,-1,-1,'MORLET');
sig95 = (sig1_u)*(ones(1,n));
sig95 = power ./ sig95;
contour(t,log2(period),sig95,[-99,1],'k');

sig5 = (sig1_l)*(ones(1,n));
sig5 = -power./sig5;
contour(t,log2(period),sig5,[-99,-1],'k--');
hold off

subplot('position',[0.77 0.37 0.2 0.28])
plot(global_ws,log2(period));
hold on
%plot(true_ws, log2(period),'r--');
plot(sig1_l, log2(period),'--');
plot(sig1_u, log2(period),'--');
%global_signif = wave_signif(1.0,dt,scales,0,0.72,-1,20,'MORLET');
%plot(global_signif,log2(period),'--')
hold off
xlabel('Power')
title('c) Global Wavelet Spectrum')
legend({'From Data', '5% Bootstrap', '95% Bootstrap'});
%set(gca, 'XLim', [0, 1.0]);
set(gca, 'YLim',log2([min(period),max(period)]), 'YDir','reverse');
%set(gca,'YLim',log2([min(period),max(period)]), ...
%	'YDir','reverse', ...
%	'YTick',log2(Yticks(:)), ...
%	'YTickLabel','');


subplot('position',[0.1 0.07 0.65 0.2])
plot(t,scale_avg)
hold on
plot(t, sig2_l, '--');
plot(t,sig2_u,'--');
legend({'From Data', '5% Bootstrap', '95% Bootstrap'},'Location','best');
set(gca,'XLim',xlim(:))
xlabel('Time (s)')
ylabel('Avg variance')
title('d) Scale-average Time Series')
hold off


