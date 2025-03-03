% this routine loads the data, and then performs a fft with a sliding time window
clear;
load filter.txt; %this is the upper portion of the section
t=filter(:,1); %meters "here:age in Ma"
%cm=cm*100; % centimeters
c=filter(:,2); %carbonate

%convert meters to age in kyr
%m=(100*cm/100); % inverse sed rate -- kyr/m


%t=m;

%convert meters to relative age in kyr
%t=(35700 + (-94.34*cm/100));

t2=t;

%interpolation to evenly spaced ages
t2=linspace(min(t),max(t),2*length(t));
c1=interp1(t,c,t2); %this is now the carbonate


%t=t2;
c=detrend(c1);

%calculate frequency array
npt=2^14;
Nyquist=0.5/abs((t(2)-t(1)));
f=linspace(0,Nyquist,npt/2);
fmax=.03;
num=round(npt/2*fmax/Nyquist);
freq=f(1:num);

% bandpass filter to remove noise
%first define the range of frequencies to keep, between f1 and f2
f1=1/22; % 1/longer period (148.4cm=140kyr)
f2=1/13;% 1/shorter period (90.1cm=85kyr)
f3=1/22;
f4=1/13;
n=length(t);
fc=fft(c);
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

%multiply by gaussian window
fc(k1:k2)=fc(k1:k2).*gausswin(length(fc(k1:k2)))';

bpwinx=linspace(f1,f2,length(fc(k1:k2)));
bpwiny=gausswin(length(fc(k1:k2)));

[m,k1]=min(abs(fre2-f3));
[m,k2]=min(abs(fre2-f4));
k3=n-k2+2;
k4=n-k1+2;
% zero all but the desired band
fc2(1:k1)=zeros(k1,1);
fc2(k2:k3)=zeros(k3-k2+1,1);
fc2(k4:n)=zeros(k1-1,1);

%take inverse FFT; remove the imaginary part
newc1=real(ifft(fc)); % gaussian windowed bp
newc2=real(ifft(fc2));
%figure(1);
%plot(t,newc1,t,newc2);
%title('BP comparison');
%legend('gaussian','rectangular')

% hilbert transforms
h1=imag(hilbert(newc1));% gaussian windowed bp
h2=imag(hilbert(newc2)); % imaginary part of the hilbert transform

t=t2;
 
h1=abs(h1);h2=abs(h2); % absolute values of the hilbert
r=findpeaks(h1);r2=findpeaks(h2); %find the x-axis values corresponding to peaks values of the previous vector
tr=t(r); tr2=t(r2); %finds the ages that correspond to the peaks
y=h1(r);y2=h2(r2); %find the y-axis values of the peaks




%interpolation to evenly spaced ages
a=linspace(min(tr),max(tr),4*length(tr));a2=linspace(min(tr2),max(tr2),4*length(tr2));
h=interp1(tr,y,a); h2=interp1(tr2,y2,a2);
%figure(2);plot(a,h,a2,h2);
%title('envelopes of band-passes');
%legend('gaussian bp','rect. bp');
hp=h; %saves the undetrended envelope for later plotting
h=detrend(h);h2=detrend(h2);

%calculate frequency array
npt=2^14;
Nyquist=0.5/abs((a(2)-a(1)));Nyquist2=0.5/abs((a2(2)-a2(1)));
f=linspace(0,Nyquist,npt/2);f2=linspace(0,Nyquist2,npt/2);
fmax=.01;fmax2=.01;
num=round(npt/2*fmax/Nyquist);num2=round(npt/2*fmax2/Nyquist2);
freq=f(1:num);freq2=f2(1:num2);

%calculate FFT and Power after padding
npt=2^14;
H=fft(h,npt);H2=fft(h2,npt);
P=H.*conj(H);P2=H2.*conj(H2);

%normalize to unit mean power
P=P/mean(P(1:npt/2));P2=P2/mean(P2(1:npt/2));

%plot
figure(3);subplot(3,1,2);plot(freq,P(1:num))
title(['FFT of Hilbert Transform GaussianBP 71-55'],'FontSize',14)
xlabel('frequency in cycles per kyr')
ylabel('spectral power')
mark1xx

F = freq';
Pow = P(1:num);
Pow = Pow';
FFTenvMAS = [F,Pow];

eq=newc1(1)-hp(1);

figure(3); subplot(3,1,3); plot(t,newc1,a,hp+eq);
xlabel('Stratigraphic height (cm)')
ylabel('BP')

%plot
%figure(4);plot(freq2,P2(1:num2))
%title(['FFT of Hilbert Transform RectangularBP'],'FontSize',14)
%xlabel('frequency in cycles per kyr')
%ylabel('spectral power')
%mark1x

%calculate FFT and Power of the whole data set after padding
npt=2^14;
H=fft(c,npt);
P=H.*conj(H);

%normalize to unit mean power
P=P/mean(P(1:npt/2));

bpwiny=bpwiny*(.5*max(P));

Nyquist=0.5/abs((t(2)-t(1)));
f=linspace(0,Nyquist,npt/2);
fmax=.6;
num=round(npt/2*fmax/Nyquist);
freq=f(1:num);



%plot
figure(3);subplot(3,1,1);plot(freq,P(1:num),bpwinx,bpwiny)
title(['Massignano 15-23m CaCo3'],'FontSize',14)
xlabel('frequency in cycles per kyr')
ylabel('spectral power')
mark1xx

%lk=ecc(3001:7001); % this creates a subset of the La2004 ecc that goes from 93-97Ma
%lk=lk/max(lk); %normalizes the ecc (not really necessary)
%newa=linspace(min(a),max(a),1916); %creates a new time vector for the Prec envelope (our section is 1.916Ma in duration)
%newh=interp1(a,h,newa); %resamples h to create more points, one for every 1kyr
%newh=newh/max(newh); %normalizes the prec envelope( not really necessary)
%[crf,lags]=xcorr(lk,newh);  % cross-correlation between laskar eccentricity and envelope of precession
%lags=(lags+93000)/1000; % converts lags to millions of years
%figure(5); plot(lags,crf); title('cross-correlation of precession envelope and lask ecc');