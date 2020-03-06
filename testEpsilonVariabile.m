clear all;
close all;

%All u's and v's needs to be interpreted as average value over the
%volumes/cj's

L=6;
m = 200; %space steps
x = linspace(0,L,m)';
T = 0.1; %final time
dx = x(3)-x(2);

global epsilon;
epsilon = 0.1*ones(length(x)-1,1);
% epsilon = 1*(x<=3) + 0.1 * (x>3);
% epsilon = epsilon(1:end-1);

if min(epsilon)<=1e-6
    dt = min(epsilon)*dx;
else
%     dt = dx^2/5;
    dt = dx^2;
end

time = [0:dt:T];
n = length(time);

delta=.5;
% u0 = sin(2*pi*x).*exp(-x.^2/(2*delta)); %initial datum
% v0=zeros(m,1);

v0 = zeros(m,1);
u0 = sin(4*pi*(x-3)).*exp(-4*(x-3).^2);

% F = @(u) (x(1:length(u)).^2-6*x(1:length(u))) .* u.^2;
% F = @(u) u.^2;
F = @(u) u;
u = u0;
v = v0;

c = epsilon.^2 ./ (epsilon.^2+dt);

figure('units','normalized','outerposition',[0 0 1 1]) %run immediately full screen

i = 2;

for t = time(1:end-1)
    y = F(u);
    vold = v;
    uhalfP = 0.5 * (u(1:end-1)+[u(2:end-1);u(1)]) - epsilon/2 .* ([v(2:end-1);v(1)]-v(1:end-1)); %Vector with N-1, u_{j+1/2}
    uhalfM = 0.5 * ([u(end-1);u(1:end-2)]+u(1:end-1)) - epsilon/2 .* (v(1:end-1)-[v(end-1);v(1:end-2)]); %Vector with N-1, u_{j-1/2}

    v(1:end-1) = v(1:end-1) .* c - 1/dx * (1-c) .* (0.5 * ([u(2:end-1);u(1)]-[u(end-1);u(1:end-2)])- ...
                 epsilon/2 .* ([v(2:end-1);v(1)]-2*v(1:end-1)+[v(end-1);v(1:end-2)]) ) + (1-c) .* y(1:end-1);
    u(1:end-1) = u(1:end-1) - dt * (c/dx .* ( 0.5*([vold(2:end-1);vold(1)]-[vold(end-1);vold(1:end-2)])  - 1./(2*epsilon).*([u(2:end-1);u(1)]-2*u(1:end-1)+[u(end-1);u(1:end-2)])  ) - ....
                (1-c) .* ( ([u(2:end-1);u(1)]-2*u(1:end-1)+[u(end-1);u(1:end-2)])/dx^2 - (1/dx) *(F(uhalfP) - F(uhalfM) ) ) );
            
    u(end) = u(1);
    v(end) = v(1);
        
    s = "Current time t = "+t;
    
    plot(x,u,'b-',x,v,'k-')
    title(s);
    legend('u','v');
    pause(0.00001);
end

