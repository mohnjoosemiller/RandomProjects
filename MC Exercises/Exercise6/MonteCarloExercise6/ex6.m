clc;clear;clf;
file = load( 'out.txt'); 
T = file(:,3); 
U = file(:,6); 
dU = file(:,7); 
errorbar ( T, U, dU, '.--')