

clear all; close all;
load('mGamma');

M = intrinsicInflation(mGamma);
J = localToGlobal(M,0);
disp('successs');
