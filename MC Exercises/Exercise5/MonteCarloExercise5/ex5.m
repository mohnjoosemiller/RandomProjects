clc;clear;clf;
file = load ( 'r_vs_nstep.txt'); 
n = file(:,3); 
r2 = file(:,1); 
error  = file(:,2);

errorbar(n,r2,error,'.','DisplayName', 'Theta = 0.9');
p = polyfit(n,r2,1);
f = @(x) polyval(p,x);
hold on 
plot(n,f(n), '--r','DisplayName', 'LinearFit');
legend('show','location', 'NorthWest');
text(200, 2.7, ['D = ' num2str(p(1),4)], 'FontWeight' ,'Bold');