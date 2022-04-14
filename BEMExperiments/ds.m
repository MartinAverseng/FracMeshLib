function [mGamma,M] = ds()


O = [0 0 0];
I = [1 0 0];
P = [0 1 0];
K = [0 0 1];
V = [O; I; P; K; -I; -P; -K];

T = [
    1 2 3;
    1 3 5;
    1 5 6;
    1 6 2;
    1 2 4;
    1 3 4;
    1 5 4;
    1 6 4;
    1 2 7;
    1 3 7;
    1 5 7;
    ];

mGamma = msh(V,T);

M = intrinsicInflation(mGamma);

end

