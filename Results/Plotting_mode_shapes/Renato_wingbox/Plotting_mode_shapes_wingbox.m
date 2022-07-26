clc;
clear;
close all;

grids = load('test/grids.txt');
connectivity = load('test/connect.txt');

% Test plot to check connectivity and grid points
% figure()
% trisurf(connectivity,grids(:,1),grids(:,2),grids(:,3));
% axis equal tight

% load reference data
ref_data = load('test/ref_mode_shapes_data_wingbox.mat');
% load initial mistuned data
ini_data = load('test/initial_mistuned_mode_shapes_data_wingbox.mat');
% load converged FEM data
fin_data = load('test/final_mode_shapes_data_wingbox.mat');

ini_ref_mode_shapes = ref_data.mode_shapes - ini_data.mode_shapes;
fin_ref_mode_shapes = ref_data.mode_shapes - fin_data.mode_shapes;

%subcase number
m = 1;
% mode number
n = 1;
% scaling factor
sf = 20;

ini_comparison_data = sqrt((sf*squeeze(ini_ref_mode_shapes(m,n,1,:))).^2 + ...
    (squeeze(sf*ini_ref_mode_shapes(m,n,2,:))).^2 + ...
    (squeeze(sf*ini_ref_mode_shapes(m,n,3,:))).^2);
ini_max = max(ini_comparison_data);
ini_comparison_data = ini_comparison_data./ini_max;

fin_comparison_data = sqrt((squeeze(sf*fin_ref_mode_shapes(m,n,1,:))).^2 + ...
    (squeeze(sf*fin_ref_mode_shapes(m,n,2,:))).^2 + ...
    (squeeze(sf*fin_ref_mode_shapes(m,n,3,:))).^2);
fin_comparison_data = fin_comparison_data./ini_max;


%% comparison plot for initial vs final FEM mode shapes

figure();
get(gca,'fontname')  % shows you what you are using.
set(gca,'fontname','times')  % Set it to times
subplot(1,2,1);
h(1) = trisurf(connectivity,sf*squeeze(ref_data.mode_shapes(m,n,1,:))+grids(:,1),...
    sf*squeeze(ref_data.mode_shapes(m,n,2,:))+grids(:,2),...
    sf*squeeze(ref_data.mode_shapes(m,n,3,:))+grids(:,3),ini_comparison_data,'FaceAlpha',0.75,'EdgeAlpha',0.5);
axis equal tight
% title('Initial vs reference')
box off
xlabel('Span [m]','FontSize',15)
ylabel('Chord [m]','FontSize',15)
zlabel('Amplitude','FontSize',15)
ax = gca;
ax.FontSize = 20;
set(gca,'FontName','Times','FontSize',28);
% axis off

subplot(1,2,2);
% figure();
h(2) = trisurf(connectivity,sf*squeeze(ref_data.mode_shapes(m,n,1,:))+grids(:,1),...
    sf*squeeze(ref_data.mode_shapes(m,n,2,:))+grids(:,2),...
    sf*squeeze(ref_data.mode_shapes(m,n,3,:))+grids(:,3),fin_comparison_data,'FaceAlpha',0.75,'EdgeAlpha',0.5);
axis equal tight
% title('Final vs reference')
% xlabel('Span [m]','FontSize',15)
% ylabel('Chord [m]','FontSize',15)
% zlabel('Amplitude','FontSize',15)
% h = colorbar;
c = colorbar;
c.FontSize = 20;
caxis([0 1.0])
% legend([h1, h2, h3], {'Reference FEM', 'Initial mistuned FEM', 'Calibrated FEM'},'FontSize',15);
hold off
% axis off
ax = gca;
ax.FontSize = 20;
get(gca,'fontname')  % shows you what you are using.
set(gca,'fontname','times')  % Set it to times
set(gca,'FontName','Times','FontSize',28);

%% mode shape plot

%subcase number
m = 1;
% mode number
n = 3;
% scaling factor
sf = 20;

figure();
h1 = trisurf(connectivity,sf*squeeze(ref_data.mode_shapes(m,n,1,:))+grids(:,1),...
    sf*squeeze(ref_data.mode_shapes(m,n,2,:))+grids(:,2),...
    sf*squeeze(ref_data.mode_shapes(m,n,3,:))+grids(:,3),'FaceColor','g','FaceAlpha',1.0);
axis equal tight
box off
hold on
h2 = trisurf(connectivity,-sf*squeeze(-ini_data.mode_shapes(m,n,1,:))+grids(:,1),...
    -sf*squeeze(-ini_data.mode_shapes(m,n,2,:))+grids(:,2),...
    -sf*squeeze(-ini_data.mode_shapes(m,n,3,:))+grids(:,3),'FaceColor','b');
h3 = trisurf(connectivity,sf*squeeze(fin_data.mode_shapes(m,n,1,:))+grids(:,1),...
    sf*squeeze(fin_data.mode_shapes(m,n,2,:))+grids(:,2),...
    sf*squeeze(fin_data.mode_shapes(m,n,3,:))+grids(:,3),'FaceColor','r');
xlabel('Span [m]','FontSize',24)
ylabel('Chord [m]','FontSize',24)
zlabel('Amplitude','FontSize',24)
legend([h1, h2, h3], {'Reference FEM', 'Initial mistuned FEM', 'Calibrated FEM'},'FontSize',24);
ax = gca;
ax.FontSize = 24; 
legend boxoff
box off
hold off

%% static displacement plot for comparing 
figure();
% add reference FEM
h1 = trisurf(connectivity,squeeze(ref_data.static_displacement(m,1,:))+grids(:,1),...
    squeeze(ref_data.static_displacement(m,2,:))+grids(:,2),...
    (squeeze(ref_data.static_displacement(m,3,:))+grids(:,3)),'FaceAlpha',0.9,'FaceColor','k','EdgeAlpha',0.1);
% axis equal tight
% xlabel('Span [m]','FontSize',15)
% ylabel('Chord [m]','FontSize',15)
% zlabel('Vertical displacement','FontSize',15)
hold on

% add initial mistuned FEM
h2 = trisurf(connectivity,squeeze(ini_data.static_displacement(m,1,:))+grids(:,1),...
    squeeze(ini_data.static_displacement(m,2,:))+grids(:,2),...
    (squeeze(ini_data.static_displacement(m,3,:))+grids(:,3)),'FaceColor','b','FaceAlpha',1.0,'EdgeAlpha',0.1);

% add final mistuned FEM
h3 = trisurf(connectivity,squeeze(fin_data.static_displacement(m,1,:))+grids(:,1),...
    squeeze(fin_data.static_displacement(m,2,:))+grids(:,2),...
    (squeeze(fin_data.static_displacement(m,3,:))+grids(:,3)),'FaceColor','r','FaceAlpha',1.0,'EdgeAlpha',0.1);
xlabel('Span [m]','FontSize',20)
ylabel('Chord [m]','FontSize',20)
zlabel('Vertical displacement [% of span]','FontSize',20)
legend([h1, h2, h3], {'Reference FEM', 'Initial mistuned FEM', 'Calibrated FEM'},'FontSize',20);
legend boxoff
box off
hold off
axis equal
zticklabels({'0' '5' '10' '15' '20'});
zlim([0 3.]);
grid off
ax = gca;
ax.FontSize = 20; 
% 
%% static displacement plot for the different load cases
figure();
% add reference FEM
% subcase 1
h1 = trisurf(connectivity,squeeze(ref_data.static_displacement(1,1,:))+grids(:,1),...
    squeeze(ref_data.static_displacement(1,2,:))+grids(:,2),...
    squeeze(ref_data.static_displacement(1,3,:))+grids(:,3),'FaceAlpha',0.5,'FaceColor','k','EdgeAlpha',0.25);
axis equal tight
xlabel('Span [m]','FontSize',15)
ylabel('Chord [m]','FontSize',15)
zlabel('Vertical displacement','FontSize',15)
hold on

% subcase 2
h2 = trisurf(connectivity,squeeze(ref_data.static_displacement(2,1,:))+grids(:,1),...
    squeeze(ref_data.static_displacement(2,2,:))+grids(:,2),...
    squeeze(ref_data.static_displacement(2,3,:))+grids(:,3),'FaceColor','b','EdgeAlpha',0.25);

% subcase 3
h3 = trisurf(connectivity,squeeze(ref_data.static_displacement(3,1,:))+grids(:,1),...
    squeeze(ref_data.static_displacement(3,2,:))+grids(:,2),...
    squeeze(ref_data.static_displacement(3,3,:))+grids(:,3),'FaceColor','[0.8500 0.3250 0.0980]','FaceAlpha',1.0,'EdgeAlpha',1.0);
xlabel('Span [m]','FontSize',15)
ylabel('Chord [m]','FontSize',15)
zlabel('Vertical displacement','FontSize',15)
axis equal tight
axis off
% legend([h1, h2, h3], {'No Load', '-1 g', '2.5 g'},'FontSize',15);
box off
hold off