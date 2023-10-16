clear

% basic domain

% circles
xc = [0.2 0.5 0.8 0.2 0.8 0.2 0.5 0.8];
yc = [0.2 0.2 0.2 0.5 0.5 0.8 0.8 0.8];
N = length(xc);

radius = 0.13;
theta = linspace(0,2*pi,200);

for n=1:N
    plot(xc(n) + radius*cos(theta), yc(n) + radius*sin(theta),'-k','LineWidth',2);
    hold on
    text(xc(n)-0.05,yc(n)+0.01,strcat('# ',num2str(n)),'FontSize',15)
end

% central square
xsquare = [0.4 0.6 0.6 0.4 0.4];
ysquare = [0.4 0.4 0.6 0.6 0.4];
plot(xsquare,ysquare,'-k','LineWidth',2);
text(0.485,0.5,'F','FontSize',15)

axis square

saveas(gcf,'cookies_domain','png')

