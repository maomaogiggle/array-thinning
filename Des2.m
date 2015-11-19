clc;
clear;

N = 41; %N antenna elements
lambda = 1;
%position of antenna array
d=zeros(1,N);
d(1) = 10;
d(2) = 9.6065;
d(3) = 8.8098;
d(4) = 8.2995;
d(5) = 7.8973;
d(6) = 7.3497;
d(7) = 6.8494;
d(8) = 6.5302;
d(9) = 5.6299;
d(10)=5.3749;
d(11)=5;
d(12)=4.6065;
d(13)=3.8098;
d(14)=3.2995;
d(15)=2.8973;
d(16)=2.3497;
d(17)=1.8494;
d(18)=1.5302;
d(19)=0.6299;
d(20)=0.3749;
d(21)=0;
d(41) = -10;
d(40) = -9.6065;
d(39) = -8.8098;
d(38) = -8.2995;
d(37) = -7.8973;
d(36) = -7.3497;
d(35) = -6.8494;
d(34) = -6.5302;
d(33) = -5.6299;
d(32)=-5.3749;
d(31)=-5;
d(30)=-4.6065;
d(29)=-3.8098;
d(28)=-3.2995;
d(27)=-2.8973;
d(26)=-2.3497;
d(25)=-1.8494;
d(24)=-1.5302;
d(23)=-0.6299;
d(22)=-0.3749;

%number of mainlobe samples
theta_m = 85:1:95;
L_m = length(theta_m);
u_m = cos(theta_m*pi/180);
A_m = exp(1i*(2*pi/lambda)*u_m'*d);
%null
theta_n = [120:0.5:130];
L_n = length(theta_n);
u_n = cos(theta_n*pi/180);
A_n = exp(1i*(2*pi/lambda)*u_n'*d);
%sidelobe
theta_s = [0:1:80,100:1:115,135:1:180];
L_s = length(theta_s);
u_s = cos(theta_s*pi/180);
A_s = exp(1i*(2*pi/lambda)*u_s'*d);
L = L_m + L_n + L_s;
%transformation matrix
theta = [theta_m,theta_s,theta_n];
u = cos(theta*pi/180);
A = exp(1i*(2*pi/lambda)*u'*d);
%desired beampattern
%amplitude
fm = [ones(L_m,1);10^(-50/20)*ones(L_s,1);10^(-80/20)*ones(L_n,1)];
%phase
fp = ones(L,1);
f = fm.*fp;
%emphasis vector
w = ones(L,1);
W = diag(w);
mu = 0.4;
%ripple tolerance
gamma = 0;
epsilon = 0.05;
km = 40;
kp = 4000;
kn = 1e6;
x0 = inv(A'*A)*A'*f;

%change weight every iteration
error = 1;
tic;
for i = 1:2
    xp = inv(A'*W*A)*(A'*W)*f;
    xp = xp./abs(xp);
    iawa = inv(real(diag(xp)'*A'*W*A*diag(xp)));
    awf = real(diag(xp)'*A'*W*f)-(1/mu)*ones(N,1);
    xa = max(iawa*awf,gamma);
    xc = xa.*xp;
    %update the weight corresponding to mainlobe
    fi = abs(A*xc)./max(abs(A*xc));
    em = fi(1:L_m)-f(1:L_m);
    for i = 1:L_m
        if (em(i) < -epsilon)
            w(i) = w(i) + km * abs(em(i));
        end
    end
    %update the weight corresponding to sidelobe
    es = fi(L_m+1:L_m+L_s)-f(L_m+1:L_m+L_s);
    w(L_m+1:L_m+L_s) = max(kp*es,0);
    %update the weight corresponding to null
    en = fi(L_m+L_s+1:end)-f(L_m+L_s+1:end);
    w(L_m+L_s+1:end) = max(kn*en,0);    
    W = diag(w);
    error = norm(x0-xc);
    x0 = xc;
end

t = toc;
%calculate the number of selected antennas
lam = 0;
for i = 1:N
    if (abs(xc(i))<0.001)
        xc(i)=0;
    else
        lam = lam + 1;
    end
end

%define the transformation matrix for plotting
theta_plot = [0:0.1:180];
u_plot = cos(theta_plot*pi/180);
A_plot = exp(1i*(2*pi/lambda)*u_plot'*d);
figure(2);hold on;
y = 20*log10(abs(A_plot*xc)/max(abs(A_plot*xc)));
plot(theta_plot,y);
hold on;
fd  =  10^(-40/20)*ones(length(theta_plot),1);
for i = 1:length(theta_plot)
    if ((theta_plot(i)>=85)&&(theta_plot(i)<=95))
        fd(i) = 1;
    elseif((theta_plot(i)>=120)&&(theta_plot(i)<=130))
        fd(i) = 10^(-60/20);
    end
end
plot(theta_plot,20*log10(fd),'r');



