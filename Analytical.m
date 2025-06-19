% Analytical solution of 1-D ADE with Laplace Transform by Bright Takyi
%Advecttion-Dispersion-Reaction Equation; s_t = -c*s_x + u*s_xx-lambs
% Initial and Boundary Conditions; 
% s(x,0)=Q0, s(0,t)=Q1, dQ/dx=0
%************************************************************************************
 
 
%Solution to Advection-Dispersion Equation with Laplace Transform
%********************************************************************
 
clear
 close all
 clc
L = 10000; % lenght of river [m]
%L = 1000; % lenght of river [m]
%oneday = 24*3600; % converting one-day into seconds
%ndays = 31;
%T = ndays*oneday; % simulation time [sec]
T = 2*3600; % simulation time [hrs]
%T = 172800; % simulation time (20 days)[s]
%oneday = 24*3600; % converting one-day into seconds
dt = 1;
dx = 20;
n = T/dt;
m = L/dx;
%tmax = 9;
%m = 500; % number of subintervals of length of river
%n = 15; % number of subintervals of simulation time
%n = 360; %1728; % number of subintervals of simulation time
%dx = L/m; % space step-size
%dt = T/n; % time step-size
%oneday = 24*3600;
x = linspace(0,L,m+1); % space discretization
x1 = x(2:end);
%t = linspace(0,T,n+1); % time discretization
t = linspace(0,T,n+1); % time discretization
%rx = dt/(2.*dx);
%rxx = dt/(dx^2);
%halflife = 3.8*24*3600;%halflife of Ra-222
%halflife = 3.8*24*3600;%halflife of Ra-222
%halflife = 3600; % halflife of Pa-230
%k = log(2)/halflife; % decay rate of Ra-222
k = 0.00002; % decay rate of Ra-222
%e = 1-k*dt;
v = 3;%.*ones(m+1,1); % wave celerity [m/s]                                        
%D1 = xlsread('C:\Users\G\Desktop\Diffusivity.xls','A4:B23');
%D=repmat(D1',m+1,1);% diffusion coefficient [m^3/s]
%D = 0.003;%.*ones(m+1,1); % diffusion coefficient [m^3/s]
D = 0.05;%.*ones(m+1,1); % diffusion coefficient [m^3/s]
% Setting up Q matrix, initial and boundary conditions
s = zeros(m+1,n+1);  
sa = 0; % initial concentration of radon-220 [mol/m^3]
sb = 100;% boundary condition at X = 0
%sb = 1;% boundary condition at X = 0
s(:,1) = sa;
s(1,:) = sb;
%sa = s(1);
s(2:end,1) = sa;
%Q(:,1) = 850; % initial condition
%Q(1,:) = 850; % boundary condition
%Qa = Q(1);
 
%for i = 1:size(t,2)
%    h = 1./(2.*D.*t(i).*sqrt(v.^2.*x.^2./(4.*D.^2) + k.*x.^2./D));
 %   p = sqrt(v.^2.*x.^2./(4.*D.^2) + k.*x.^2./D);
  %  r = (v.*x./2.*D);
   % s(:,i) = e.*(exp(-k.*t(i))) + exp(r-p).*erfc(x.^2./h)-exp(r-p).*e.*(exp(-k.*t(i))).*erfc(x.^2./h);
    
%end
 
p = ((v.^2.*x1.^2./(4.*D.^2))+(k.*x1.^2)./D);
q = x1.^2./D;
n = v.*x1./(2.*D);
%h = p.*x1./(sqrt(D));
h = q./(2.*(sqrt(p)));
%h = x.^2./(D.*(sqrt(p)));
for i = 2:size(t,2)
    
    A = sb.*exp(n-(p.^(1./2))).*erfc(h./(2.*t(i)));
    B = -sa.*exp(n-(p.^(1./2))).*exp(-k.*t(i)).*erfc(h./(2.*t(i)));
    C = sa.*exp(-k.*t(i));
    
    
    s(2:end,i) = A' + B' + C';
end
s = s./sb;
spA = s(51,:);
save spA
%figure
%plot(t(1:7201),spA(1:7201))
%xlabel('Time [s]')
%ylabel('Radionuclide Concentration s/so [Bq]')
%title('Graph of Single Radionuclide Concentration profile')
%s(1,:) = sb;
%plot(s,t)
%xlabel('time [sec]')
%ylabel('concentration [Bq]')
%title('A Graph of concentration vs time')
%figure
%plot(t(1:7201),s(1,1:7201),'LineWidth',3);
%plot(t(1:7201),s(51,1:7201),'LineWidth',3);
%hold on
%plot(t(1:7201),s(151,1:7201),'LineWidth',3);
%plot(t(1:7201),s(301,1:7201),'LineWidth',3);
%plot(t(1:7201),s(501,1:7201),'LineWidth',3);
%xlabel('Time [s]')
%ylabel('Relative Concentration [s/so]')
%legend('x = 1000 m','x = 3000 m', 'x = 6000 m','x = 10000 m','location','SouthEast')
 
 
figure
%plot(t(1:7201),s(1,1:7201),'LineWidth',3);
plot(t(1:7201),s(151,1:7201),'LineWidth',3);
hold on
plot(t(1:7201),s(251,1:7201),'LineWidth',3);
plot(t(1:7201),s(301,1:7201),'LineWidth',3);
plot(t(1:7201),s(401,1:7201),'LineWidth',3);
xlabel('t/tmax [s]')
ylabel('Relative Concentration [s/so]')
legend('x = 50 m','x = 75 m', 'x = 100 m','x = 10000 m','location','SouthEast')


%figure
%plot(x,s(:,601),'LineWidth',3);
%hold on 
%plot(x,s(:,1801),'LineWidth',3);
%plot(x,s(:,3001),'LineWidth',3);
%plot(x,s(:,7201),'LineWidth',3);
%xlabel('Distance from initial section [m]')
%ylabel('Relative Concentration [s/so]')
%legend('After 600 s','After 1800 s', 'After 3000 s', 'After 7200 s','SouthWest')

%figure
%plot(x,s(:,1301),'LineWidth',3);
%hold on 
%plot(x,s(:,1301),'LineWidth',3);
%plot(x,s(:,3001),'LineWidth',3);
%plot(x,s(:,4001),'LineWidth',3);
%xlabel('Distance from initial section [m]')
%ylabel('Relative Concentration [s/so]')
%legend('After 5 days','After 31 days', 'After 31 days', 'After 7200 s','SouthWest')



