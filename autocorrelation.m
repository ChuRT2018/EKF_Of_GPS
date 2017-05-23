n = length(NEU);
Rneu = ones(3,n);

Mn = mean(NEU(1,:));
Me = mean(NEU(2,:));
Mu = mean(NEU(3,:));
Dn = var(NEU(1,:));
De = var(NEU(2,:));
Du = var(NEU(3,:));



for m = 1 : n
    SUMn = 0;
    SUMe = 0;
    SUMu = 0;
    for k = m+1 : n
        SUMn = SUMn + (NEU(1,k) - Mn) * (NEU(1,k-m) - Mn);
        SUMe = SUMe + (NEU(2,k) - Me) * (NEU(2,k-m) - Me);
        SUMu = SUMu + (NEU(3,k) - Mu) * (NEU(3,k-m) - Mu);
    end
    Rneu(1,m) = SUMn / (n - m) / Dn;
    Rneu(2,m) = SUMe / (n - m) / De;
    Rneu(3,m) = SUMu / (n - m) / Du;
end
