function [As,Bs,Cs]=AnalyticalSolutionGeneralComplex(kD,m,n,A0,B0)

% This gives the equ concentration of mA+nB<->C where kD=kd/ka
%labeled As, Bs and Cs


p = @(C) ((A0-m*C)^m)*((B0-n*C)^n)-kD*C;

Cs = fzero(p,min(A0/m,B0/n)); 

As=A0-m*Cs;

Bs=B0-n*Cs; 
