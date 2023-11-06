H = [0.2 0.1 0.05 0.025 0.0125 0.01];
L2 = [0.2074 0.0371 0.0077 0.0019 4.7143e-04 3.0223e-04];
H1 = [0.1135 0.0328 0.0090 0.0026 8.7547e-04 6.2295e-04];
Y1 = log(L2);
Y2 = log(H1);
X = log(1./H);
p1 = polyfit(X,Y1,1)
p2 = polyfit(X,Y2,1)
y1 = p1(1)*X+p1(2);
y2 = p2(1)*X+p2(2);
figure(1)
hold on
plot(X,Y1,'r-*')
plot(X,Y2,'b-o')
plot(X,y1,'r--')
plot(X,y2,'b--')
xlabel('log(1/h)');
ylabel("logarithme de l'erreur");
legend({'norme L2','norme H1'},'Location','southwest')