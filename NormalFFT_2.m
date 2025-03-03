clear;
clf;
tmp=load('subtracted119-99.txt'); %input data file
%tmp=load('DAY_pdt.txt'); %input data file
c=tmp(:,2); %name the color vector c
h=tmp(:,1); %name the height vector h

set_fmax=.03;
figure(4);
plot(h,c);

%interpolation to evenly spaced ages
a=linspace(min(h),max(h),2*length(h));
cc=interp1(h,c,a);


%subtract the mean then plot
cc=cc-mean(cc);
figure(1);
plot(a,cc);
%zoom on;

%calculate FFT and Power after padding
npt=2^14;
H=fft(cc,npt);
P=H.*conj(H);
%calculate frequency array
Nyquist=abs(0.5/(a(2)-a(1)));
f=linspace(0,Nyquist,npt/2);
%normalize to unit mean power
P=P/mean(P(1:npt/2));
%plot
fmax=set_fmax;
num=round(npt/2*fmax/Nyquist);
figure(2);
plot(f(1:num),P(1:num));
xlabel('frequency')
ylabel('spectral power')
title('DAY2')

mark1xx

t=h;

cm=h;

t=(cm*1); % conversion to time in kyr assuming that meter level 0 is
% 36000 kyr and using an inverse sed rate of 41kyr/46cm


%cm=cm*100; % centimeters

%convert meters to age in kyr
%m=(100*cm/100); % inverse sed rate -- kyr/m
%t=m;
%convert meters to relative age in kyr
%t=(35700 + (-94.34*cm/100));




figure(1); plot(t,c); legend('DAY data','95% confidence');

%interpolation to evenly spaced ages
t2=linspace(min(t),max(t),2*length(t));
c1=interp1(t,c,t2); %this is now the carbonate

t=t2;
c=detrend(c1);

figure(2); plot(t2,c); legend('DAY time series','95% confidence');

%calculate frequency array
npt=2^14;
Nyquist=0.5/abs((t2(2)-t2(1)));
f=linspace(0,Nyquist,npt/2);
fmax=set_fmax
num=round(npt/2*fmax/Nyquist);
freq=f(1:num);


%calculate FFT and Power after padding
H=fft(c,npt);
P=H.*conj(H);

%normalize to unit mean power
P=P/mean(P(1:npt/2));

%calculate the noise with Monte Carlo
t0=t; % just so I don't mess up the t vector
t3=linspace(min(t0),max(t0),2*length(t0))';
%calculate frequency array
Nyquist3=0.5/abs(t3(2)-t3(1)); 
npt=2^14;
f=linspace(0,Nyquist3,npt/2);
num3=round(npt/2*fmax/Nyquist3);
freq3=f(1:num3);
%
steps=1000; %number of MC simulations 
Power=zeros(num3,steps); %creates an empty array for the spectral power
w=0; % a counter

for i=1:steps
	r=rand([(length(t0)) 1]); % generates a string of random numbers with a certain spacing
	hh=interp1(t0,r,t3); % interpolates between random numbers
    % this above step is called for because with the carbonate data, I did
    % a linear interpolation and increased the number of points by a factor
    % of 2, just to give better spectral resolution -- so I do the same
    % thing with the random number strings
	hh=hh-mean(hh); % subtracts the mean
	w=w+1; %advances the counter
	%calculate FFT and Power after padding
	H3=fft(hh,npt)'; % puts the result into a column vector
	P3=H3.*conj(H3); % gets rid of the imaginary part
	%normalize to unit mean power
	P3=P3/mean(P3(1:npt/2));
	Power(:,w)=P3(1:num3)'; % loads the Power column into the larger array
end
PQ=mean(Power'); % calculates the mean of each column of the transpose of Power

PR=std(Power');  % calculates the standard deviation of each column of the transpose of Power
PQ3=2*PQ; % THREE TIMES THE MEAN  -- this is what theoretically gives you the 95% confidence level
PQ3b=3*PQ; % THREE TIMES THE MEAN  -- this is what theoretically gives you the 95% confidence level
PQ3c=4*PQ; % THREE TIMES THE MEAN  -- this is what theoretically gives you the 95% confidence level, 4 is 99%, 2 is 90%?
%plot
figure(3);plot(freq,P(1:num),freq3,PQ3,freq3,PQ3b,freq3,PQ3c)
legend('GMS lab FFT2','90% confidence','95% confidence','99% confidence')
xlabel('Frequency in cycles per cm')
ylabel('Spectral power')
mark1xx

% bandpass filter to remove noise
%first define the range of frequencies to keep, between f1 and f2
f1=1/830; %1600 1/longer period
f2=1/350;%350 1/shorter period
f3=1/10;%350
f4=1/70;%150
n=length(a);
fc=fft(cc);
fc2=fc;
fre=linspace(0,2*Nyquist,n)';
fre2=linspace(0,2*Nyquist,n)';
[m,k1]=min(abs(fre-f1));
[m,k2]=min(abs(fre-f2));
k3=n-k2+2;
k4=n-k1+2;
% zero all but the desired band
fc(1:k1)=zeros(k1,1);
fc(k2:k3)=zeros(k3-k2+1,1);
fc(k4:n)=zeros(k1-1,1);
%take inverse FFT; remove the imaginary part
newc1=real(ifft(fc));

[m,k1]=min(abs(fre2-f3));
[m,k2]=min(abs(fre2-f4));
k3=n-k2+2;
k4=n-k1+2;
% zero all but the desired band
fc2(1:k1)=zeros(k1,1);
fc2(k2:k3)=zeros(k3-k2+1,1);
fc2(k4:n)=zeros(k1-1,1);

%take inverse FFT; remove the imaginary part


newc2=real(ifft(fc2));
newc3=newc1+newc2;
figure(5);
plot(a,cc,a,newc1,a,newc2,a,newc3);set(gca,'XDir','reverse');
title([' DAY = ',num2str(1/f1),' to ',num2str(1/f2),' and ',num2str(1/f3),' to ',num2str(1/f4),' units?'],'FontSize',12)

%calculate difference between data and 

dd=cc-newc1;
ee=cc-newc2;
ddd=[a',dd'];
eee=[a',ee'];
fff=[a',newc1'];
save('subtracted1.txt','ddd','-ascii');
save('subtracted2.txt','eee','-ascii');
save('filter.txt','fff','-ascii');

figure(6);
subplot(2,2,1)       % add first plot in 2 x 2 grid
plot(a,cc,a,newc1);set(gca,'XDir','reverse');           % line plot
title('Data with 1 filter')

subplot(2,2,3)       % add second plot in 2 x 2 grid
plot(ddd(:,1),ddd(:,2));set(gca,'XDir','reverse');        % scatter plot
title('Residual 1')

subplot(2,2,2)       % add third plot in 2 x 2 grid
plot(a,cc,a,newc2);set(gca,'XDir','reverse');           % stem plot
title('Data with 2 filter')

subplot(2,2,4)       % add fourth plot in 2 x 2 grid
%yyaxis left          % plot against left y-axis  
plot(eee(:,1),eee(:,2));set(gca,'XDir','reverse');%xlim([0 3500])              
%yyaxis right         % plot against right y-axis
%plot(x,y2)
title('Residual 2')

fff=[a',newc1'];
save('filter.txt','fff','-ascii');