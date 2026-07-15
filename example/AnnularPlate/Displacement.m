clear
clc
PathStr='./data/';
textFiles=dir([PathStr '*.dat']);
%xloc(1:16)=[362,540,718,895,1073,1251,1429,1607,1784,1962,2140,2318,2496,2673,2851,3029];
xloc(1:3)=[5,210,216];


fidA=dlmread ("./post/nodeData/node_6.dat");
fidB=dlmread ("./post/nodeData/node_0.dat");
verA=dlmread ("./JB_Data/vertical_A.csv");

verB=dlmread ("./JB_Data/vertical_B.csv");
MITC3verA=dlmread ("./JB_Data/MITC3_vertical_A.csv");

[n,m]=size(fidA);

for i=1:n
 % disp(i)=sqrt((9.96533-fid(i,2))^2+(0.0168766+fid(i,3))^2+(fid(i,4))^2)-10;
end
force = 0.032:0.032:3.2;
force = force/3.2;
%force = force/4;

figure
plot(abs(fidA(:,4)),force,'LineWidth',3)
xlabel('Displacement')
ylabel('Point Load (P/P_{max})')
hold on
plot(abs(fidB(:,4)),force,'LineWidth',3)
hold on
plot(verA(:,1),verA(:,2),'LineWidth',3)
hold on
plot(verB(:,1),verB(:,2),'LineWidth',3)
hold on
plot(MITC3verA(:,1),MITC3verA(:,2),'LineWidth',3)
legend('vertical A','vertical B','JB vertical A','JB vertical B','MITC3 vertical A','location','se')
%xlim([0 6])
%ylim([0 1])