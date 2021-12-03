clc;
clear;
close all;

% data_lin_v_nl = [
% 1	0.677	0.683	-0.83	0.676	0.23
% 2	2.873	2.873	0.03	2.876	-0.10
% 3	4.252	4.270	-0.42	4.243	0.21
% 4	11.707	11.734	-0.24	11.688	0.16
% 5	17.203	17.197	0.04	17.196	0.04
% 6	18.044	18.035	0.05	18.048	-0.02
% 7	22.179	22.218	-0.17	22.155	0.11
% 8	34.814	34.906	-0.26	34.800	0.04
% 9	47.895	47.871	0.05	47.928	-0.07
% 10	48.229	48.577	-0.72	48.254	-0.05
% 11	52.229	52.236	-0.01	52.181	0.09
% 12	54.887	58.010	-5.69	55.070	-0.33
% 13	60.567	62.543	-3.26	60.642	-0.12
% 14	63.920	63.999	-0.12	63.879	0.06
% 15	67.155	69.970	-4.19	67.252	-0.14];
% 
% data_SD_wo_w = [0.14	0.15	0.16	0.14	0.04	0.10	0.14	0.06
%                 0.03	0.07	0.04	0.04	0.31	3.01	2.94	0.29
%                 0.05	0.00	0.06	0.01	0.34	0.34	0.40	4.98
%                 0.07	0.06	0.02	0.07	0.34	0.35	0.35	5.18
%                 0.05	0.02	0.03	0.07	0.33	0.44	0.58	0.30
%                 0.01	0.10	0.15	0.03	0.12	0.41	1.40	0.23
%                 0.01	0.02	0.04	0.01	4.63	0.18	0.23	0.10
%                 0.05	0.05	0.05	0.05	0.29	0.31	0.31	0.29
%                 0.14	0.06	0.07	0.01	2.10	0.27	0.28	0.36
%                 0.06	0.17	0.09	0.06	0.26	0.84	0.68	0.26
%                 0.03	0.03	0.03	0.03	0.15	0.15	0.15	0.15
%                 0.05	0.05	0.04	0.20	0.27	0.28	0.28	2.10
%                 0.04	0.36	0.04	0.05	0.12	2.56	0.28	0.27
%                 0.05	0.05	0.26	0.05	0.25	0.27	1.78	0.25
%                 0.17	0.03	0.02	0.24	3.56	0.19	0.19	1.78];
% 
% str = ["Mode" "Reference" "single"  "Percent diff(ref_v_single)"  "multi" "Percent diff(ref_v_multi)"];
% %% plot frequencies for linear vs nonlinear optimization
% figure()
% hold on
% plot(data_lin_v_nl(1:end,1),data_lin_v_nl(:,2),'k.','MarkerSize',20);
% plot(data_lin_v_nl(1:end,1),data_lin_v_nl(:,3),'b.','MarkerSize',20);
% plot(data_lin_v_nl(1:end,1),data_lin_v_nl(:,5),'r.','MarkerSize',20);
% ylabel('Frequency [Hz]','FontSize',14); xlabel('Mode number','FontSize',14)
% legend('Reference','Typical FEM updating methodolgy','Proposed FEM updating methodology')
% ax = gca;
% ax.FontSize = 24;
% xticks(linspace(1,15,15));
% 
% %% plot errors for including static deflection in the objective function
% figure()
% hold on
% plot(linspace(1,15,15),data_SD_wo_w(:,4),'k.','MarkerSize',20);
% plot(linspace(1,15,15),data_SD_wo_w(:,8),'r.','MarkerSize',20);
% ylabel('Error [%]','FontSize',14); xlabel('Mode number','FontSize',14)
% legend('Baseline method','Static deflection in objective function')
% ax = gca;
% ax.FontSize = 24;
% xticks(linspace(1,15,15));

%% static displacement plot for comparing single vs multi

grids = load('test/grids.txt');
connectivity = load('test/connect.txt');

% load reference case data
ref_data = load('test/reference_case_complex_load.mat');
% load single case data
single_data = load('test/single_case_complex_load.mat');
% load multi-case data
multi_data = load('test/multi_case_complex_load.mat');

%subcase number
m = 1;

figure();
% add reference FEM
h1 = trisurf(connectivity,squeeze(ref_data.static_displacement(m,1,:))+grids(:,1),...
    squeeze(ref_data.static_displacement(m,2,:))+grids(:,2),...
    squeeze(ref_data.static_displacement(m,3,:))+grids(:,3),'FaceAlpha',0.5,'FaceColor','k','FaceAlpha',1.0,'EdgeAlpha',0.0);
axis equal tight
xlabel('Span [m]','FontSize',15)
ylabel('Chord [m]','FontSize',15)
zlabel('Vertical displacement','FontSize',15)
hold on

% add single case FEM
h2 = trisurf(connectivity,squeeze(single_data.static_displacement(m,1,:))+grids(:,1),...
    squeeze(single_data.static_displacement(m,2,:))+grids(:,2),...
    squeeze(single_data.static_displacement(m,3,:))+grids(:,3),'FaceColor','b','FaceAlpha',1.0,'EdgeAlpha',0.0);

% add multi case FEM 
h3 = trisurf(connectivity,squeeze(multi_data.static_displacement(m,1,:))+grids(:,1),...
    squeeze(multi_data.static_displacement(m,2,:))+grids(:,2),...
    squeeze(multi_data.static_displacement(m,3,:))+grids(:,3),'FaceColor','r','FaceAlpha',0.5,'EdgeAlpha',0.0);
axis equal tight
xlabel('Span [m]','FontSize',15)
ylabel('Chord [m]','FontSize',15)
zlabel('Vertical displacement','FontSize',15)
legend([h1, h2, h3], {'Reference FEM', 'Single FEM', 'Multi FEM'},'FontSize',15);
box off
hold off