function [X]=sim_VARMA_2_1(N)
%
% simulate VARMA(2,1) process directly 
%
% coefficients are those found in Eichler (2008) JMVA, 99, 968-1009
%
% N    number of sample vectors to generate
%
%randn('state',2) % state must be defined in calling program
p=3;
Phi_1=[1.3 0.2 0.0; 0 -1.4 0.1; 0 0 0.7];
Phi_2=[-0.9 -0.1 0; 0 -0.9 0.1; 0 0 -0.9];
Theta_1=[0.5 0 0; 0 0.5 0; 0 0 -0.5];

nsim=N+1000;% use series much longer than needed, throw away first part
z(1:p,1)=zeros;
z(1:p,2)=zeros;
epsold(1:p,1)=zeros;
for j=3:nsim
    epsnew=randn(p,1);
    z(1:p,j)=Phi_1*z(1:p,j-1)+Phi_2*z(1:p,j-2)+epsnew+Theta_1*epsold;
    epsold=epsnew;
end
X(1:p,1:N)=z(1:p,1001:N+1000);
