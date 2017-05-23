function [ B,L,H ] = XYZ2BLH( XC,YC,ZC )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    a = 6378137;
    b = 6356752.3142;
    ep2 = 0.00673949674227;
    e2 = 0.00669437999013;

    B = zeros(length(XC),1);
    L = B;
    H = B;

    parfor i = 1 : length(XC)

        Sxy = sqrt(XC(i) * XC(i) + YC(i) * YC(i));

        u = atan(a * ZC(i) * (1 + b * ep2 / sqrt(XC(i) * XC(i) + YC(i) * YC(i) + ZC(i) * ZC(i))) / (b * Sxy));

        Su = sin(u);
        Cu = cos(u);

        B(i) = atan((ZC(i) + ep2 * b * Su * Su * Su) / (Sxy - e2 * a * Cu * Cu * Cu));

        Sb = sin(B(i));

        L(i) = atan(YC(i) / XC(i));
        H(i) = Sxy * cos(B(i)) + ZC(i) * Sb - a * sqrt(1 - e2 * Sb * Sb);

        B(i) = B(i) * 180 / pi;
        L(i) = L(i) * 180 / pi;
    end

end

