% The following program solves the  Poisson equation describing
% the two-dimensional steady state groundwater head profile 
% throughout th island which is approximated as a rectangle with 
% x=600m and y=1000m 

% Define variables
xlength = 600;  % length in x- direction (m) 
ylength = 1000;  % length in y-direction (m)
del = 50;       % delx = dely = 50m
d=del^2; %del^2=2500m^2
xnodes=xlength/del -1;
ynodes=ylength/del -1;
totalnodes=xnodes*ynodes;
% Let the transmissivity be 200 m2/day
T=200;
% Let the initial aquifer recharge be 40cm/year
R=40;
S=R*0.01/365;
V=d*S/T;

% Allocate space for matrices
h=zeros(totalnodes,1);
rhs=(-V-10/d)*ones(totalnodes,1);

% Fill in rhs matrix for right BC
for i=1:ynodes
        rhs(i*xnodes)=(-V-10/d); 
end

% the four corners of the island and the top and bottom rows 
% will have the following boundary conditions

rhs(1)=-V-20/d;
rhs(2:10)=-V-10/d;
rhs(11)=-V-20/d;
rhs(199)=-V-20/d;
rhs(200:208)=-V-10/d;
rhs(209)=-V-20/d;

  
% Fill in coefficient matrix (as sparse matrix)
A=sparse(totalnodes);
M_diag=sparse(1:totalnodes,1:totalnodes,-4,totalnodes,totalnodes);
L_diag_1=sparse(2:totalnodes,1:totalnodes-1,1,totalnodes,totalnodes);
L_diag_2=sparse(xnodes+1:totalnodes,1:totalnodes-xnodes,1,totalnodes,totalnodes);
U_diag_1=L_diag_1';
U_diag_2=L_diag_2';
A=M_diag+L_diag_1+L_diag_2+U_diag_1+U_diag_2;

% Modify coefficient matrix for BCs
for i=1:ynodes-1
    A(i*xnodes,i*xnodes+1)=0;   % right BC
    A(i*xnodes+1,i*xnodes)=0;   % 
end

% Solve for h
h_old=A\rhs;


Sy=0.8;

% Define the End time
tmax=550;

% Define a hh matrix to store
hh=sparse(totalnodes,tmax);

nt=551;
dt=tmax/(nt-1);
r=(T*dt)/(Sy*d);

% Fill in coefficient matrix (as sparse matrix)
G=sparse(totalnodes);
M_diag=sparse(1:totalnodes,1:totalnodes,1+2*r,totalnodes,totalnodes);
L_diag_1=sparse(2:totalnodes,1:totalnodes-1,-r/2,totalnodes,totalnodes);
U_diag_1=sparse(1:totalnodes-1,2:totalnodes,-r/2,totalnodes,totalnodes);
L_diag_2=L_diag_1;
U_diag_2=U_diag_1;
G=M_diag+L_diag_1+U_diag_1+L_diag_2+U_diag_2;

% Modify coefficient matrix for BCs
for i=1:ynodes-1
    G(i*xnodes,i*xnodes+1)=0;   % right BC
    G(i*xnodes+1,i*xnodes)=0;   % 
end

% Fill in coefficient matrix (as sparse matrix)
H=sparse(totalnodes);
M_diag=sparse(1:totalnodes,1:totalnodes,1-2*r,totalnodes,totalnodes);
L_diag_1=sparse(2:totalnodes,1:totalnodes-1,r/2,totalnodes,totalnodes);
U_diag_1=sparse(1:totalnodes-1,2:totalnodes,r/2,totalnodes,totalnodes);
L_diag_2=L_diag_1;
U_diag_2=U_diag_1;
H=M_diag+L_diag_1+U_diag_1+L_diag_2+U_diag_2;

% Modify coefficient matrix for BCs
for i=1:ynodes-1
    H(i*xnodes,i*xnodes+1)=0;   % right BC
    H(i*xnodes+1,i*xnodes)=0;   % 
end

n=dt;

for n=1:dt:tmax
    h_new=G\(H*h_old);
    hh(:,n)=h_new;
    h_old=h_new;
end

% Specify the time at which you want the groundwater profile
k=input('Please enter the time at which you want the groundwater profile (specify an integer value between 1 and 550 days)''\n')

h=hh(:,k)+10;

% Reshape h into grid and plot
hh_new=reshape(h,xnodes,ynodes)';
xplot=del:del:xlength-del;
yplot=del:del:ylength-del;
mesh(xplot,yplot,hh_new)
xlabel('x distance, m')
ylabel('y distance, m')
zlabel('Groundwater head')
title({'Groundwater head at time after'; k; 'days'})