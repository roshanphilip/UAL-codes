clc
clear all
close all
format compact

% This script can be used to generate crossplots for various 1D and 2D problems 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: Name your input files as Data n_e x where x ranges from 1 to n
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0,'DefaultAxesFontSize', 16, 'DefaultAxesLineWidth', 1, ...
'DefaultLineLineWidth', 1.5, 'DefaultAxesFontName', 'Latin Modern Math', ...
'DefaultTextFontName', 'Latin Modern Math', ...
'DefaultTextFontSize', 16, ...
'defaultAxesTickLabelInterpreter', 'latex', ...
'defaultLegendInterpreter', 'latex', ...
'defaultTextInterpreter', 'latex', ...
'defaultColorbarTickLabelInterpreter', 'latex', ...
'defaultPolaraxesTickLabelInterpreter', 'latex', ...
'defaultTextarrowshapeInterpreter', 'latex', ...
'defaultTextboxshapeInterpreter', 'latex', ...
'DefaultLegendBox','on', 'DefaultLegendFontSize', 16, ...
'DefaultAxesBoxStyle', 'back', 'DefaultAxesBox', 'off')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. Choose correct file path 
cd 'C:\Users\rs7625\Desktop\Plotting function\FD plot'

% 2. Load data files
load(['Data n_e 1.mat']);
x1=plot_storage((1:end),1);
y1=plot_storage((1:end),2);
plot_x1_length=size(x1,1);
% time_1=runtime
load(['Data n_e 2.mat']);
x2=plot_storage((1:end),1);
y2=plot_storage((1:end),2);
plot_x2_length=size(x2,1);
% time_2=runtime
load(['Data n_e 3.mat']);
x3=plot_storage((1:end),1);
y3=plot_storage((1:end),2);
plot_x3_length=size(x3,1);
% time_3=runtime
load(['Data n_e 4.mat']);
x4=plot_storage((1:end),1);
y4=plot_storage((1:end),2);
plot_x4_length=size(x4,1);
% time_4=runtime

% 3. Figure name
file_name='UAL, 1D, Local, Coarse, $\phi$=0.80, DL=4mm, tolerance study';
save_file_name='UAL, 1D, Local, Coarse, phi=0.80, DL=4mm, tolerance study';
f = figure("Position",[300 100 720 500]); % The first 2 arguments define the position, the second 2 arguments define the size
t = tiledlayout(1,1,'TileSpacing','compact','Padding','compact'); 
title0 = title(file_name);
set(title0,'interpreter','latex','fontsize', 16,'Color','k');
hold on; box on; grid on;

% 4. Axis label
xlabel ('Displacement (mm)','FontSize', 16)
ylabel ('Force (N)','FontSize', 16)

% 5. Axis limits
x_end = 0.016;
y_end = 3;
xlim([0 x_end])
ylim([0 y_end])
xticks(0:2e-3:x_end)
yticks(0:0.5:y_end)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =============================== Plots ==================================%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A. Plotting data with markers
plot(x1,y1,'--ro','LineWidth',1,'Color','red','LineStyle','-','Marker','square','MarkerSize',9,'MarkerIndices',[1:round(0.01*plot_x1_length):plot_x1_length]);
plot(x2,y2,'--ro','LineWidth',1.5,'Color','black','LineStyle','-.','Marker','o','MarkerSize',9,'MarkerIndices',[1:round(0.01*plot_x2_length):plot_x2_length]);
plot(x3,y3,'--ro','LineWidth',2,'Color','blue','LineStyle','--','Marker','v','MarkerSize',9,'MarkerIndices',[1:round(0.01*plot_x3_length):plot_x3_length]);
plot(x4,y4,'--ro','LineWidth',2.5,'Color',[0.4940 0.1840 0.5560],'LineStyle',':','Marker','diamond','MarkerSize',9,'MarkerIndices',[1:round(0.01*plot_x4_length):plot_x4_length]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% B. Plotting data without markers
% plot(x1,y1,'LineWidth',1.25,'Color','red','Marker','.','LineStyle','-','MarkerSize',20,'MarkerIndices',[1:round(0.01*plot_x2_length):plot_x1_length]);
% plot(x2,y2,'LineWidth',2,'Color','black','LineStyle','-.','MarkerSize',9,'MarkerIndices',[1:round(0.01*plot_x2_length):plot_x2_length]);
% plot(x3,y3,'LineWidth',3,'Color','blue','LineStyle',':','MarkerSize',9,'MarkerIndices',[1:round(0.01*plot_x3_length):plot_x3_length]);
% 
% plot(x4,y4,'LineWidth',4,'Color','green','LineStyle','--','MarkerSize',9,'MarkerIndices',[1:round(0.01*plot_x4_length):plot_x4_length]);
% plot(x5,y5,'LineWidth',5,'Color','yellow','LineStyle','-.','MarkerSize',9,'MarkerIndices',[1:round(0.01*plot_x5_length):plot_x5_length]);
% plot(x6,y6,'LineWidth',6,'Color','magenta','LineStyle',':','MarkerSize',9,'MarkerIndices',[1:round(0.01*plot_x6_length):plot_x6_length]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C. Self comparison plots
% xlabel ('Displacement (mm)','FontSize', 16)
% ylabel ('Force (N)','FontSize', 16)
% plot(x1,y1,'--ro','LineWidth',2,'Color',[0.8500 0.3250 0.0980],'LineStyle','-','Marker','square','MarkerSize',9,'MarkerIndices',[1:round(0.01*plot_x1_length):plot_x1_length]);
% % plot(x2,y2,'--ro','LineWidth',2,'Color','black','LineStyle','-','Marker','o','MarkerSize',9,'MarkerIndices',[1:round(0.01*plot_x2_length):plot_x2_length]);
% plot(x3,y3,'--ro','LineWidth',2,'Color',[0.4940 0.1840 0.5560],'LineStyle','-','Marker','v','MarkerSize',9,'MarkerIndices',[1:round(0.01*plot_x3_length):plot_x3_length]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ============================= Annotations ==============================%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X1. Zoom into a portion of the graph
% 
% Use the code below to zoom into a particular region of the graph 
% Zoom rectangle
% box_x_left             = 7.8e-3;
% box_y_bottom           = 2.2;
% box_width              = 1e-3;
% box_height             = 0.2;
% zoom_box_x_left_bottom = 0.65;
% zoom_box_y_left_bottom = 0.4;
% zoom_box_width         = 0.2;
% zoom_box_height        = 0.275;
% rectangle('Position',[box_x_left box_y_bottom box_width box_height],'LineStyle','--','LineWidth',1)
% % 
% % Zoom line 1 
% x_line_1=[box_x_left                     (x_end * zoom_box_x_left_bottom) + 0.3e-3] 
% y_line_1=[box_y_bottom                   (y_end * zoom_box_y_left_bottom) - 0.15] 
% line(x_line_1,y_line_1,'LineStyle',':','HandleVisibility','off','Color','black')
% 
% % Zoom line 2
% x_line_2=[box_x_left   + box_width       (zoom_box_x_left_bottom + zoom_box_width)  * x_end + 1.2e-3] 
% y_line_2=[box_y_bottom + box_height      (zoom_box_y_left_bottom + zoom_box_height) * y_end + 0.05] 
% line(x_line_2,y_line_2,'LineStyle',':','HandleVisibility','off','Color','black')

%-----------------------------------------------------------------------------
% Zoomed in part of the graph (1st window)
% axes('position',[zoom_box_x_left_bottom zoom_box_y_left_bottom zoom_box_width zoom_box_height])
% hold on; box on % put box around new pair of axes
% plot(x1,y1,'LineWidth',1,'Color','red','LineStyle','-','Marker','square','MarkerSize',9,'MarkerIndices',[1:round(0.01*plot_x2_length):plot_x1_length])
% plot(x2,y2,'LineWidth',2,'Color','black','LineStyle','-.','Marker','o','MarkerSize',9,'MarkerIndices',[1:round(0.01*plot_x2_length):plot_x2_length])
% plot(x3,y3,'LineWidth',3,'Color','blue','LineStyle','--','Marker','v','MarkerSize',9,'MarkerIndices',[1:round(0.01*plot_x3_length):plot_x3_length])
% ax = gca;
% ax.XLim = [box_x_left box_x_left+box_width];
% ax.YLim = [box_y_bottom box_y_bottom+box_height];

% text(0.0059, 8.233106111949236e+02, 'C', 'Interpreter', 'tex')
% text(0.00595, 9.105681377644195e+02, 'A', 'Interpreter', 'tex')
% text(0.006, 8.697217406279459e+02, 'B', 'Interpreter', 'tex')

%-----------------------------------------------------------------------------
% Zoomed in part of the graph (2nd window)
% axes('position',[.6 .25 .2 .2])
% hold on; box on % put box around new pair of axes
% plot(x1,y1,'LineWidth',2,'Color','red','LineStyle','-','Marker','square','MarkerSize',9,'MarkerIndices',[1:round(0.01*plot_x1_length):plot_x1_length])
% plot(x2,y2,'LineWidth',2,'Color','black','LineStyle','-','Marker','o','MarkerSize',9,'MarkerIndices',[1:round(0.01*plot_x2_length):plot_x2_length])
% plot(x3,y3,'LineWidth',3,'Color','blue','LineStyle','-','Marker','v','MarkerSize',9,'MarkerIndices',[1:round(0.01*plot_x3_length):plot_x3_length])
% ax = gca;
% ax.XLim = [1.5e-3 1.75e-3];
% ax.YLim = [2250 2750];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X2. Text marker annotation 
% text(0.004, 1.025e+03, 'B', 'Interpreter', 'tex')
% text(0.0052, 1.425e+03, 'A', 'Interpreter', 'tex') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X3. Legends
% 1.
% legend(["UAL w/o shear","UAL with shear","NR w/o shear"],'FontSize',16,'position',[.7 .3 .05 .05]);
% 2.
% legend(["UAL","FAL","NR","UAL_old","FAL_old","NR_old"],'FontSize',16,'position',[.7 .55 .05 .05]);
% 3. 
% legend(["Coarse","Fine"],'FontSize',16,'position',[.7 .3 .05 .05]);
% 4. 
% legend(["Analytical tangent","Numerical tangent"],'FontSize',16,'position',[.65 .4 .05 .05]);
% 5. 
legend(["tol=1e-6","tol=1e-8","tol=1e-10","tol=1e-12"],'FontSize',16,'position',[.7 .5 .05 .05]);
% 6.
% legend(["$\Delta \lambda_0 = 1e-3$","$\Delta \lambda_0 = 1e-4$"],'FontSize',14,'position',[.75 .5 .05 .05]);
% 7. 
% legend(["Displacement driven","Traction driven"],'FontSize',14,'position',[.7 .5 .05 .05]);
% 8. 
% legend(["UAL","FAL","NR"],'FontSize',16,'position',[.25 .675 .05 .05]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X4. Choose scientific notation type 
ax = gca;
ax.YAxis.Exponent = 0;
% ax.XAxis.Exponent = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save image
% saveas(gcf,strcat(save_file_name,'.eps'), 'eps')
saveas(gcf,strcat(save_file_name,'.pdf'), 'pdf')
saveas(gcf,strcat(save_file_name,'.png'), 'png')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

