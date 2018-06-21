r = 100;
alpha = 0.1;
count = 0;

for i = 1:(r-1)
    for j = (i+1):r
        count = count + 1;
        T(count,1) = i;
        T(count,2) = j;
        T(count,3) = randn();
    end
end

P = sortrows(T,3);

exit = 0;
l = 0;
L = length(T(:,1));
split = L + 1;

while (exit == 0 && l < L)
    l = l + 1;
    CritLevel = norminv(1 - (alpha/(L-l+1)),0,1);
    %%CritLevel = norminv((1 - alpha)^(1/(L-l+1)),0,1);
    if P(L+1-l,3) < CritLevel
        exit = 1;
        split = l;
    end
end

E = ones(r,r);

if split ~= L + 1
    for i = 1:(L-split+1)
        E(P(i,1),P(i,2)) = 0;
        E(P(i,2),P(i,1)) = 0;
    end
end
