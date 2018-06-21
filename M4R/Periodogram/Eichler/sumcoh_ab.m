function QQ=sumcoh_ab(S)
%
% computes the statistic Q_T on p981 of Eichler(2008)
% Spectral matrix is structured as
% S(1:p,1:p,1:(n/2)+1)
%
[p,p,nreq]=size(S);
QQ(p,p)=zeros;
gamsq(1:nreq)=zeros;

% Precalculate the inverse of S
Sinv = zeros(size(S));
for j=1:nreq
    Sinv(:,:,j)=inv(squeeze(S(:,:,j)));
end

for a=1:p-1
    for b=a+1:p
        Q=0;
        for j=1:nreq
            gamsq(j)=real(abs(Sinv(a,b,j)).^2./( Sinv(a,a,j).*Sinv(b,b,j)));
        end            
%        gamsq(1:nreq)=abs(S(a,b,1:nreq)).^2./( S(a,a,1:nreq).*S(b,b,1:nreq));
        Q=sum(gamsq);% covers 0 to Nyquist, now do rest [negative freqs]
        Q=Q+sum(gamsq(2:nreq-1));
        Q=Q./(2.*(nreq-1));
         QQ(a,b)=Q;
    end
   
end