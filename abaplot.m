clear all;
clc;
yy=[0 4.61548 11.9563 22.3358 43.7058 88.888 186.263 398.98 440.571];
xx=[1 1.593255 1.846082 2.003516 2.161809 2.31815 2.470736 2.61858 2.637441];
plot(xx,yy,'-*','linewidth',1);
xlabel('Primary Stretch','FontSize',18,'Color','r')
ylabel('Principal Stress','FontSize',18, 'Color','r')
hold on;
lam=[1:0.1:2.7];
sig=zeros(1,length(lam))
for i=1:length(lam)
    sig(i)=1.5*(lam(i)^2*(1/lam(i)))
end
plot(lam,sig,'linewidth',1.5);
legend('Abaqus results','Analytical results');
hold on;