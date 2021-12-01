clc;
clear;
close all;

figure()
load('test/grids.txt')
load('test/connect.txt')
trisurf(connect,grids(:,1),grids(:,2),grids(:,3));
axis equal tight

%subcase number
m = 3;
% mode number
n = 5;
% scaling factor
sf = 20;

%% mode shape plot
figure();
% add reference data
load('test/ref_mode_shapes_data_wingbox.mat');
h1 = trisurf(connect,sf*squeeze(mode_shapes(m,n,1,:))+grids(:,1),...
    sf*squeeze(mode_shapes(m,n,2,:))+grids(:,2),...
    sf*squeeze(mode_shapes(m,n,3,:))+grids(:,3),'FaceColor','k','FaceAlpha',0.5);
axis equal tight
box off
hold on

% add initial mistuned data
load('test/initial_mistuned_mode_shapes_data_wingbox.mat');
h2 = trisurf(connect,-sf*squeeze(mode_shapes(m,n,1,:))+grids(:,1),...
    -sf*squeeze(mode_shapes(m,n,2,:))+grids(:,2),...
    -sf*squeeze(mode_shapes(m,n,3,:))+grids(:,3),'FaceColor','b');

% add converged FEM data
load('test/final_mode_shapes_data_wingbox.mat');
h3 = trisurf(connect,sf*squeeze(mode_shapes(m,n,1,:))+grids(:,1),...
    sf*squeeze(mode_shapes(m,n,2,:))+grids(:,2),...
    sf*squeeze(mode_shapes(m,n,3,:))+grids(:,3),'FaceColor','r');
xlabel('Span [m]','FontSize',15)
ylabel('Chord [m]','FontSize',15)
zlabel('Amplitude','FontSize',15)
legend([h1, h2, h3], {'Reference FEM', 'Initial mistuned FEM', 'Calibrated FEM'},'FontSize',15);
hold off

%% static displacement plot for comparing 
figure();
% add reference FEM
load('test/ref_mode_shapes_data_wingbox.mat');
h1 = trisurf(connect,squeeze(static_displacement(m,1,:))+grids(:,1),...
    squeeze(static_displacement(m,2,:))+grids(:,2),...
    squeeze(static_displacement(m,3,:))+grids(:,3),'FaceAlpha',0.5,'FaceColor','k');
axis equal tight
xlabel('Span [m]','FontSize',15)
ylabel('Chord [m]','FontSize',15)
zlabel('Vertical displacement','FontSize',15)
hold on

% add initial mistuned FEM
load('test/initial_mistuned_mode_shapes_data_wingbox.mat');
h2 = trisurf(connect,squeeze(static_displacement(m,1,:))+grids(:,1),...
    squeeze(static_displacement(m,2,:))+grids(:,2),...
    squeeze(static_displacement(m,3,:))+grids(:,3),'FaceColor','b');

% add final mistuned FEM
load('test/final_mode_shapes_data_wingbox.mat');
h3 = trisurf(connect,squeeze(static_displacement(m,1,:))+grids(:,1),...
    squeeze(static_displacement(m,2,:))+grids(:,2),...
    squeeze(static_displacement(m,3,:))+grids(:,3),'FaceColor','r');
xlabel('Span [m]','FontSize',15)
ylabel('Chord [m]','FontSize',15)
zlabel('Vertical displacement','FontSize',15)
legend([h1, h2, h3], {'Reference FEM', 'Initial mistuned FEM', 'Calibrated FEM'},'FontSize',15);
box off
hold off

%% static displacement plot for the different load cases
figure();
% add reference FEM
load('test/ref_mode_shapes_data_wingbox.mat');

% subcase 1
h1 = trisurf(connect,squeeze(static_displacement(1,1,:))+grids(:,1),...
    squeeze(static_displacement(1,2,:))+grids(:,2),...
    squeeze(static_displacement(1,3,:))+grids(:,3),'FaceAlpha',0.5,'FaceColor','k');
axis equal tight
xlabel('Span [m]','FontSize',15)
ylabel('Chord [m]','FontSize',15)
zlabel('Vertical displacement','FontSize',15)
hold on

% subcase 2
h2 = trisurf(connect,squeeze(static_displacement(2,1,:))+grids(:,1),...
    squeeze(static_displacement(2,2,:))+grids(:,2),...
    squeeze(static_displacement(2,3,:))+grids(:,3),'FaceColor','b');

% subcase 3
load('test/final_mode_shapes_data_wingbox.mat');
h3 = trisurf(connect,squeeze(static_displacement(3,1,:))+grids(:,1),...
    squeeze(static_displacement(3,2,:))+grids(:,2),...
    squeeze(static_displacement(3,3,:))+grids(:,3),'FaceColor','r');
xlabel('Span [m]','FontSize',15)
ylabel('Chord [m]','FontSize',15)
zlabel('Vertical displacement','FontSize',15)
legend([h1, h2, h3], {'No Load', '-1 g', '2.5 g'},'FontSize',15);
box off
hold off