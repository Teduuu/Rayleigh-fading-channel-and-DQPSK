clc;
M = csvread('dbpsk_data.csv');
t = M(:,1);
rn = M(:,2);
phi = M(:,3);
hi_t = M(:,4);
hq_t = M(:,5);
%[a,xi] = ksdensity(rn);
%figure;
%plot(xi,a);
figure(1);
semilogy(t,rn);
title('amplitude of r_n');
xlabel('t (s)');
ylabel('r_n');
xlim([0 0.15]);


figure(2);


p2 = histogram(rn,'Normalization','pdf','BinWidth',0.005);
p2.EdgeColor = "#0072BD";
title('distribution of r_n');
ylabel('probability density');
xlim([-1 3]);
hold on
x = [0:0.01:3];
p = raylpdf(x,0.707);
plot(x,p,'r');
hold off

figure(3);
p3 = histogram(phi,'Normalization','pdf','BinWidth',0.005);
p3.EdgeColor = "#0072BD";
title('distribution of \phi _n');
ylabel('probability density');
xlim([-4 4]);
ylim([0 0.4]);


figure(4);
p4 = histogram(hi_t,'Normalization','pdf','BinWidth',0.005);
p4.FaceColor = "#0072BD";
p4.EdgeColor = "#0072BD";
title('distribution of h_I_n');
ylabel('probability density');
xlim([-3 3]);
hold on
x2 = [-3:0.01:3];
p2 = normpdf(x2,0,0.707);
plot(x2,p2,'r');
hold off

figure(5);
p5 = histogram(hq_t,'Normalization','pdf','BinWidth',0.005);
p5.EdgeColor = "#0072BD";
title('distribution of h_Q_n');
ylabel('probability density');
xlim([-3 3]);
hold on
x3 = [-3:0.01:3];
p3 = normpdf(x3,0,0.707);
plot(x3,p3,'r');
hold off
