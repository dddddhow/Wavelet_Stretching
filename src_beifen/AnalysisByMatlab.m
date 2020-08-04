clc,clear;
%% 
nt=1000; %时间采样点数（个）
dt=0.002; %时间采样间隔（秒）

nw=81; %子波长度（点）

v=1000;%速度（米秒）

h=100; %反射界面深度（米）
t0=2.0*h/1.0/v;%反射界面旅行时（秒）
nh=t0/dt;%反射界面点数（个）

dx=5;%横向采样间隔（米）
dz=5;%纵向采样间隔（米）
nx=1001;%横向采样（个）
nz=nt*dt*v/dz/2;%纵向采样（个）


%% CMP道集
ix=-(nx-1)/2*dx:dx:(nx-1)/2*dx;
seis_xt_1=sqrt(t0^2+4.0*ix.*ix/v/v);
seis_xt_2=sqrt((t0+nw*dt)^2+4.0*ix.*ix/v/v);

figure(1)
subplot(121)
plot(ix,seis_xt_1,"linewidth",2),hold on;
plot(ix,seis_xt_2,'r',"linewidth",2),legend('Upper Bound','Lower Bound'),
xlabel("Offset / m"),ylabel("Time / s"),
ylim([0 nt*dt]),
set(gca,'linewidth',2,'fontsize',30,'fontname','Times New Roman');
axis ij

subplot(122)
plot(ix,seis_xt_1+nw*dt-seis_xt_2,"linewidth",2),
legend('Subtracting Lower Bound from Upper Bound'),
set(gca,'linewidth',2,'fontsize',30,'fontname','Times New Roman'),
xlabel("Offset / m"),ylabel("Time / s"),
ylim([0 nt*dt]),
axis ij
% 


%% T-Theta道集

n_theta=nx;
seis_thetat_1=zeros(1,n_theta);
seis_thetat_2=zeros(1,n_theta);
theta_max=atan(((nx-1)/2*dx)/(h/2));
d_theta=theta_max/n_theta;

theta_max_here_1=atan(((nx-1)/2*dx)/(h));
theta_max_here_2=atan(((nx-1)/2*dx)/(h+nw*dt*v/2));
theta_max_here=theta_max_here_1;

if theta_max_here_2 < theta_max_here_1
    theta_max_here=theta_max_here_2;
end

for ix=1:nx
    
    theta=d_theta*(ix-(nx-1)/2);
    
    if abs(theta)<theta_max_here_1
        seis_thetat_1(ix)=sqrt(t0^2+4.0*h*h*sin(theta).*sin(theta)./cos(theta)./cos(theta)/v/v);
    end
    
    if abs(theta)<theta_max_here_2
        seis_thetat_2(ix)=sqrt((t0+nw*dt)^2+4.0*(h+nw*dt*v/2)*(h+nw*dt*v/2)*sin(theta).*sin(theta)./cos(theta)./cos(theta)/v/v);
    end   
end

figure(2)
subplot(121)
plot(((1:n_theta)*d_theta-(n_theta-1)/2*d_theta)/pi*180,seis_thetat_1,"linewidth",2),hold on;
plot(((1:n_theta)*d_theta-(n_theta-1)/2*d_theta)/pi*180,seis_thetat_2,"linewidth",2),hold on;
legend('Upper Bound','Lower Bound'),
xlabel("Angle / degree"),ylabel("Time / s"),
ylim([0 nt*dt]),
set(gca,'linewidth',2,'fontsize',30,'fontname','Times New Roman');
axis ij

subplot(122)
plot(((1:n_theta)*d_theta-(n_theta-1)/2*d_theta)/pi*180,seis_thetat_1+nw*dt-seis_thetat_2,"linewidth",2),hold on;
legend('Subtracting Lower Bound from Upper Bound'),
xlabel("Angle / degree"),ylabel("Time / s"),
ylim([-50*dt 50*dt]),
set(gca,'linewidth',2,'fontsize',30,'fontname','Times New Roman');
axis ij


