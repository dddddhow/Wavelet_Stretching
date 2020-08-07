clc,clear;
dt=0.001;
nx=201;
nt=2000;
dx=5;
v=2000;

for ix=1:nx
    for it=1:nt
        t0=it*dt;
        x=(ix-nx/2)*dx;
        t_xt(it,ix)=sqrt(t0^2+4*x^2/v/v);
    end
end

figure(1)
for it=1:100:nt
subplot(121),plot(1:nx,t_xt(it,:)),hold on,
end
ylim([0 nt*dt]),
title("X-T");
axis ij
for ix=1:10:100
subplot(122),plot(1:nt,t_xt(:,ix)),hold on
end
title("T0-T");
axis ij

figure(2)
surfc(t_xt)