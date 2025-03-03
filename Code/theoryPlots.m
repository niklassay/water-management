clc;
clear;
f = @(t) sin(t*pi);
x1 = linspace(0,1, 101);
x2 = linspace(0,1, 11);
figure
hold on
bar(x2,f(x2), 1, "w")
plot(x2,f(x2), "o", "color", "k")
plot(x1,f(x1), "color", "#0072BD")
ylabel('f(t)');
xlabel('t');
set(gca,'FontSize',8)
set(gcf,'Units','Centimeters','position',[0,0,16,6]);
saveas(gcf,'quadratureA.eps', 'epsc');
sum(f(linspace(0,1,101))/101)