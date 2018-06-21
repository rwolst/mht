function y=BartlettPriestley_p4(x)
%
% this is k^4 for BartlettPriestley and is zero after \pm 2
%
if x==0
    y=1;
else
    y=(3./((pi.*x).^2)).*( (sin(pi.*x)./(pi.*x))-cos(pi.*x));
    y=y.^4;
end
