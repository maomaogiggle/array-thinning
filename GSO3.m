%Fast Iterative Shrinkage Method
%unsymmetric linear array

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
theta_m = [85:1:95,120:1:130];
L_m = length(theta_m);
u_m = cos(theta_m*pi/180);
A_m = exp(1i*(2*pi/lambda)*u_m'*d);
%sidelobe
theta_s = [0:1:80,100:1:115,113.8,136.8,135:1:180];
L_s = length(theta_s);
u_s = cos(theta_s*pi/180);
A_s = exp(1i*(2*pi/lambda)*u_s'*d);
L = L_m + L_s;
%transformation matrix
theta = [theta_m,theta_s];
u = cos(theta*pi/180);
A = exp(1i*(2*pi/lambda)*u'*d);
%desired beampattern
%amplitude
fm = [ones(L_m,1);10^(-40/20)*ones(L_s,1)];
%phase
fp = ones(L,1);
f = fm.*fp;
%emphasis vector
w = fm.^(-1);

%transform from complex domain to real domain
Ar = [real(A),-imag(A);imag(A),real(A)];
wr = [w;w];
Mr = diag(wr);
x0 = zeros(2*N,1);
xc = inv(A'*A)*A'*f;
%tau
tau = 2*max(eig(Ar'*Mr*Ar));
%mu
mu = 0.5;
%ripple tolerance
epsilon = 0.05;
km = 40;
kp = 4000;

%change weight every iteration
error = 1;
y1 = x0;
t1 = 1;
tic;
for k = 1:5
    while (error > 1e-4)
        %calculating the phase of the desired beampattern
        fp = A*xc;
        fp = fp./abs(fp);
        fr = [fm.*real(fp);fm.*imag(fp)];
        %shrinkage
        g = y1-(1/tau)*Ar'*Mr*(Ar*y1-fr);
        x = zeros(2*N,1);
        for i = 1:N
            r = [g(i);g(i+N)] ;
            z = max(norm(r)-1/(tau*mu),0)*(r/norm(r));
            x(i) = z(1);
            x(i+N) = z(2);
        end
        t = (1+sqrt(1+4*t1*t1))/2;
        y = x + (t1-1)/t*(x-x0);
        error = norm(x0-x); 
        x0 = x;
        y1 = y;
        t1 = t;       
        xc = x0(1:N)+1i*x0(N+1:2*N);  
    end
    %update the weight corresponding to mainlobe
    fi = abs(A*xc)./max(abs(A*xc));
    em = fi(1:L_m)-fm(1:L_m);
    for i = 1:L_m
        if (em(i) < -epsilon)
            w(i) = w(i) + km * abs(em(i));
        end
    end
    %update the weight corresponding to sidelobe
    es = fi(L_m+1:end)-fm(L_m+1:end);
    w(L_m+1:end) = max(kp*es,0);  
    wr = [w;w];
    Mr = diag(wr);
    %tau
    tau = 2*max(eig(Ar'*Mr*Ar));
    error = 1;
end

t = toc;
%calculate the number of selected antennas
lam = 0;
for i = 1:N
    if (abs(xc(i))<0.005)
        xc(i)=0;
    else
        lam = lam + 1;
    end
end

%define the transformation matrix for plotting
theta_plot = [0:0.1:180];
u_plot = cos(theta_plot*pi/180);
A_plot = exp(1i*(2*pi/lambda)*u_plot'*d);
figure(1);hold on;
y = 20*log10(abs(A_plot*xc)/max(abs(A_plot*xc)));
plot(theta_plot,y);
hold on;
fd  =  10^(-30/20)*ones(length(theta_plot),1);
for i = 1:length(theta_plot)
    if ((theta_plot(i)>=85)&&(theta_plot(i)<=95))
        fd(i) = 1;
    elseif((theta_plot(i)>=120)&&(theta_plot(i)<=130))
        fd(i) = 1;
    end
end
plot(theta_plot,20*log10(fd),'r');

