clear
clc
PathStr='./data/';
textFiles=dir([PathStr '*.dat']);
%xloc(1:16)=[362,540,718,895,1073,1251,1429,1607,1784,1962,2140,2318,2496,2673,2851,3029];
xloc(1:3)=[5,210,216];


fid=dlmread ("./post/nodeData/node_10.dat");

ver=dlmread ("./jb_data/vertical.csv");

hor=dlmread ("./jb_data/horizontal.csv");

[n,m]=size(fid);

for i=1:n
  disp(i)=sqrt((fid(i,2))^2+(fid(i,3))^2+(fid(i,4))^2);
end
disp;
force = 200:200:2000;
%force = 100:100:2000;
force = force/2000;
figure
plot(abs(fid(:,4)),force,'LineWidth',3)
xlabel('Displacement')
ylabel('Point Load (P/P_{max})')
hold on
plot(abs(fid(:,3)),force,'LineWidth',3)
hold on
plot(hor(:,1),hor(:,2),'LineWidth',3)
hold on
plot(ver(:,1),ver(:,2),'LineWidth',3)
legend('vertical','horizontal','JB horizontal','JB vertical','location','se')