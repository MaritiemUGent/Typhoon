fid = fopen('calc_equi_07072020.txt','rt');
C = textscan(fid, '%f%f%f%f%f%f');
fclose(fid);
A = cell2mat(C)
[s1 s2] = size(A);
A(:,2) = A(:,2).*180/pi;
for i=1:s1

clf 
figure(1)
% plot3(A(1:i,1),A(1:i,2),A(1:i,3));
% hold on
% plot3(A(1:i,1),A(1:i,2),A(1:i,4));
% xLimits = get(gca,'XLim');  yLimits = get(gca,'YLim');  zLimits = get(gca,'ZLim');
% [X,Y]=meshgrid(xLimits(1):0.001:xLimits(2),yLimits(1):0.001:yLimits(2));
% Z=Y.*0;                  
% s=surf(X,Y,Z,'EdgeColor','none','FaceAlpha',0.8,'FaceColor','b');
clf
figure(1)                  
subplot(1,2,1)
plot(A(1:i,1),A(1:i,3));
ylabel('Z-force [N]');
xlabel('Relative elavation [m]')
%xlim([-0.4 -0.2])
hold on
scatter(A(i,1),A(i,3));
% subplot(2,2,2)
% plot(A(1:i,2),A(1:i,3));
% hold on
% scatter(A(i,2),A(i,3));
% subplot(2,2,3)
% plot(A(1:i,1),A(1:i,4));
% hold on
% scatter(A(i,1),A(i,4));
subplot(1,2,2)

plot(A(1:i,2),A(1:i,4));
ylabel('Moment around y-axis [Nm]');
xlabel('Relative pitch [deg.]')
%xlim([-0.5 1.5])
hold on
scatter(A(i,2),A(i,4));

set(gcf,'color','w');

pause



end    