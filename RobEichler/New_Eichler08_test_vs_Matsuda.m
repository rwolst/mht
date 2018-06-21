% this compares Eichler's 2008 results (testing partial noncorrelation)
% versus Matsuda


%
% get zero mean X
%
clear all
close all

% uses VARMA(2,1) model from Eichler 2008

deltat=1;
p=3; 


percent=0.2;
%percent=0;
nsims=10000;
% specify partial correlation indices to examine
a=1; b=3;

Eichler=0; 

if Eichler==1 %%%%%%%%%% Eichler test statistic %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% generate VAR(1) series of length N
%
% to choose M run VARMA_2_1_test_freq_wftd_Eichler08.m

%N=4000; M=55;

%N=2000; M=100;
N=2000; M=30;
%N=200; M=4;
%N=2000; M=20;
%N=2000; M=80;
Q(1:nsims)=zeros;
randn('seed', 1);
count=0;
for j=1:nsims

[X]=sim_VARMA_2_1(N);

[p,N]=size(X);

[Sd C_h] = direct_spectrum(X, percent, deltat);
%
%get weights
%
[w C_w2 C_w4]=quadratic_weights(M);
%
% now weighted direct spectrum
%
S  = freq_wgtd_dir_spec(Sd, w);
%
% now compute Eichler's statistic Q_T
%
QQ=sumcoh_ab(S);
%S_T=S_T*2*pi;
%
% form Q_T
%
% Note B_T=2M/N so M_T=1/B_T=N/(2M)
M_T=N./(2.*M);
Q(j)= N.*QQ(a,b)-M_T.*(C_h.*C_w2);
Q(j)=Q(j)./(C_h.*sqrt(M_T.*2.*C_w4));

if abs(Q(j)) > 1.96
    count=count+1;
end

end
'rejection rate'
rejection_rate=100*count/nsims

figure(1)
fig=1;
qqplot_stdnorm(Q, fig)

exportfig(gcf,'Eichler_statN2000M30.eps','FontMode','fixed',...
            'FontSize',7,'Height',2.5,'Width',2.5);



else %%%%%%%%%%%%%%% Matsuda test statistic %%%%%%%%%%%%%%%%%%%%%%

%N=500; m=40;
% to choose m run VARMA_2_1_test_wgtd_period_Eichler08.m

N=2000; m=40;
T(1:nsims)=zeros;
randn('seed', 1);
count=0;
for j=1:nsims

[X]=sim_VARMA_2_1(N);
%
% now carry out Matsuda test 
%
% define weights
for k=-m/2:m/2
    w(k+ (m/2)+1)=cos(pi*k/m);
end

S = WtdPdgm(X, w, deltat);

%
% run CompCuDu to get constants C and D
%
C=0.6169;
D=0.4464;
%
% set up matrices such that 0 indicates no connection
%
E1(1:p,1:p)=ones; %  saturated
E2(1:p,1:p)=ones; 
E2(a,b)=0;
E2(b,a)=0;
T(j) = real( TestStat( X,E1,E2,S,m,C,D ));
if abs(T(j)) > 1.96
    count=count+1;
end

end

'rejection rate'
rejection_rate=100*count/nsims

fig=1;
qqplot_stdnorm(T, fig)

exportfig(gcf,'Matsuda_statN2000m40.eps','FontMode','fixed',...
            'FontSize',7,'Height',2.5,'Width',2.5);

end

