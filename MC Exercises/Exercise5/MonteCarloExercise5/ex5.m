clc;clear;clf;
file = load ( 'r_vs_nstep.txt'); 
n = file(:,3); 
r2 = file(:,1); 
error  = file(:,2);

errorbar(n,r2,error,'.');