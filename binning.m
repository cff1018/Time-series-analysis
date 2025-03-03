clear;
clf;

load 'GMSSUSC.txt'
a=GMSSUSC(:,1);
b=GMSSUSC(:,2);

[c,ic,ia]=unique(a);
sums=accumarray(ia,b);
wgts=accumarray(ia,ones(size(b)));
avgs=sums./wgts;

%data_file=load('haf.txt')
day4a=c;
day4b=avgs;


dt_day4 = detrend(day4b,'linear',1);
trend = day4b - dt_day4;

opol = 4;
[p,s,mu] = polyfit(day4a,day4b,opol);
[f_y,delta] = polyval(p,day4a,s,mu);
dt_ecgnl = day4b - f_y;


%plot(day4a,dt_ecgnl,'-',day4a,dt_ecgnl + 2 * delta,':',day4a,dt_ecgnl - 2 * delta,':'), grid

%plot(day4a,f_y,'g');set(gca,'XDir','reverse');

subplot(3,1,1);
plot(a,b,'.');set(gca,'XDir','reverse');

subplot(3,1,2);
plot(c,avgs,'-r',day4a,f_y,'g');set(gca,'XDir','reverse');

subplot(3,1,3);
plot(day4a,dt_ecgnl,'-');set(gca,'XDir','reverse');

out=[c,avgs];
out2=[c,dt_ecgnl];


save('GMSSUSC_binned.txt','out','-ascii');
save('GMSSUSC_binned_detrended_degree4.txt','out2','-ascii');