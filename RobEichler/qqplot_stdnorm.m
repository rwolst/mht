function qqplot_stdnorm(x, f)
%
% produces a qqplot for standard normal testing
%
n=length(x);
p(1:n)=zeros;
x=sort(x);
p=(1:n)./(n+1); % use positions j/(n+1)
q=norminv(p,0,1); % quantiles of N(0,1)
u=max(abs(x));
v=max(abs(q));
w=max(u,v);
w=ceil(w);

figure(f)

plot(q,x,'+')
axis([-w w -w w])
hold on
plot([-w w],[-w w],'-.')
xlabel('standard normal quantiles')
ylabel('sample quantiles')