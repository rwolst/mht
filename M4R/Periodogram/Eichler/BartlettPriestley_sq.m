function y=BartlettPriestley_sq(x)
%
% this is k^2 for BartlettPriestley and is zero after \pm 4
%
if x==0
    y=1;
else
    y=(3./((pi.*x).^2)).*( (sin(pi.*x)./(pi.*x))-cos(pi.*x));
    y=y.^2;
end
