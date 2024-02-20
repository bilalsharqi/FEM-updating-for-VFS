clc;
close all;
set(0, 'DefaultAxesFontName', 'Times');
clear;

Errors = [11.37,0.47,4.52,0.12,0.16;
                    11.70,0.04,0.66,0.25,0.03;
                    10.22,0.09,1.75,0.88,0.27;
                    12.10,1.91,3.66,0.89,1.31;
                    7.12,1.04,5.99,1.04,1.34];

lab_gvt_ini = Errors(:,1);
lab_gvt_fin = Errors(:,2);
wt_gvt_ini = Errors(:,3);
wt_gvt_fin = Errors(:,4);
lab_gvt_fin_split = Errors(:,5);

figure()
x = [1 2 3 4 5];
bar(x,Errors(:,1),0.4);
hold on
bar(x,Errors(:,2),0.4);

bar(x,Errors(:,5),0.4);
% bar(Errors(:,1),Errors(:,4))
xlabel('Mode number')
% xticks([1:15]);
ylabel('Error [%]')
set(gca,'FontSize',20)
ylim([0,13])
legend('Initial FEM','Updated FEM - constant properties',...
    'Updated FEM - variable properties','Location','NorthWest')
legend boxoff
set(gcf,'Position',[100 100 1200 700])
box off
hold off