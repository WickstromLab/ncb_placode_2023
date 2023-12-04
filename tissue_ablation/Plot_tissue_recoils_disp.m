function Plot_tissue_recoils_disp (File2);

MaskPath = [File2.pathname, filesep,'Boundary'];
DataPath = [File2.pathname, filesep,'Quantifications'];
ImagesPath = File2.ImPath ;

frames = File2.NFiles.Beads-1 ;

load([DataPath,filesep,'Means1.mat']);
load([DataPath,filesep,'Ang_diff1.mat']);
load([DataPath,filesep,'C1A_dat.mat']);
load([DataPath,filesep,'C1B_dat.mat']);
load([DataPath,filesep,'Can1_dat.mat']);
load([DataPath,filesep,'Can1B_dat.mat']);

figh1=figure;

subplot(2,4,1);
plot(sqrt(Means1(:,1).^2+Means1(:,2).^2));
hold on;
plot([5.5 5.5],[0 3],'--k','Linewidth',0.5);
plot([0 frames],[mean(sqrt(Means1(1:4,1).^2+Means1(1:4,2).^2)) mean(sqrt(Means1(1:4,1).^2+Means1(1:4,2).^2))],'--r','Linewidth',0.5);
title('Fibro A meandisp');
set(gca,'YLim',[0 max(Means1(:))]);
set(gca,'XLim',[0 frames]);

subplot(2,4,2);
plot(sqrt(Means1(:,4).^2+Means1(:,5).^2));
hold on;
plot([5.5 5.5],[0 3],'--k','Linewidth',0.5);
plot([0 frames],[mean(sqrt(Means1(1:4,4).^2+Means1(1:4,5).^2)) mean(sqrt(Means1(1:4,4).^2+Means1(1:4,5).^2))],'--r','Linewidth',0.5);
title('Fibro B meandisp');
set(gca,'YLim',[0 max(Means1(:))]);
set(gca,'XLim',[0 frames]);

subplot(2,4,3);
plot(sqrt(Means1(:,7).^2+Means1(:,8).^2));
hold on;
plot([5.5 5.5],[0 3],'--k','Linewidth',0.5);
plot([0 frames],[mean(sqrt(Means1(1:4,7).^2+Means1(1:4,8).^2)) mean(sqrt(Means1(1:4,7).^2+Means1(1:4,8).^2))],'--r','Linewidth',0.5);
title('Placode meandisp');
set(gca,'YLim',[0 max(Means1(:))]);
set(gca,'XLim',[0 frames]);

subplot(2,4,4);
plot(sqrt(Means1(:,10).^2+Means1(:,11).^2));
hold on;
plot([5.5 5.5],[0 3],'--k','Linewidth',0.5);
plot([0 frames],[mean(sqrt(Means1(1:4,10).^2+Means1(1:4,11).^2)) mean(sqrt(Means1(1:4,10).^2+Means1(1:4,11).^2))],'--r','Linewidth',0.5);
title('Placode CT meandisp');
set(gca,'YLim',[0 max(Means1(:))]);
set(gca,'XLim',[0 frames]);

subplot(2,4,5);
plot(Ang_diff1(:,5));
hold on;
plot([5.5 5.5],[-180 180],'--k','Linewidth',0.5);
plot([0 frames],[90 90],'--r','Linewidth',0.5);
plot([0 frames],[0 0],'--r','Linewidth',0.5);
plot([0 frames],[-90 -90],'--r','Linewidth',0.5);
title('Fibro A angle');
set(gca,'YLim',[-180 180]);
set(gca,'XLim',[0 frames]);

subplot(2,4,6);
plot(Ang_diff1(:,6));
hold on;
plot([5.5 5.5],[-180 180],'--k','Linewidth',0.5);
plot([0 frames],[90 90],'--r','Linewidth',0.5);
plot([0 frames],[0 0],'--r','Linewidth',0.5);
plot([0 frames],[-90 -90],'--r','Linewidth',0.5);
title('Fibro B angle');
set(gca,'YLim',[-180 180]);
set(gca,'XLim',[0 frames]);

subplot(2,4,7);
plot(Ang_diff1(:,7));
hold on;
plot([5.5 5.5],[-180 180],'--k','Linewidth',0.5);
plot([0 frames],[90 90],'--r','Linewidth',0.5);
plot([0 frames],[0 0],'--r','Linewidth',0.5);
plot([0 frames],[-90 -90],'--r','Linewidth',0.5);
title('Placode angle ');
set(gca,'YLim',[-180 180]);
set(gca,'XLim',[0 frames]);

subplot(2,4,8);
plot(Ang_diff1(:,8));
hold on;
plot([5.5 5.5],[-180 180],'--k','Linewidth',0.5);
plot([0 frames],[90 90],'--r','Linewidth',0.5);
plot([0 frames],[0 0],'--r','Linewidth',0.5);
plot([0 frames],[-90 -90],'--r','Linewidth',0.5);
title('Placode CT angle ');
set(gca,'YLim',[-180 180]);
set(gca,'XLim',[0 frames]);

print( figh1, '-djpeg90', '-r300', [ ImagesPath, filesep,'Plots',filesep, 'Velocities_angles_can1.jpg' ] ) ;
close ;

figh1=figure;

subplot(2,4,1);
plot(sqrt(Means1(:,1).^2+Means1(:,2).^2));
hold on;
plot([5.5 5.5],[0 3],'--k','Linewidth',0.5);
plot([0 frames],[mean(sqrt(Means1(1:4,1).^2+Means1(1:4,2).^2)) mean(sqrt(Means1(1:4,1).^2+Means1(1:4,2).^2))],'--r','Linewidth',0.5);
title('C1A meandisp');
set(gca,'YLim',[0 max(Means1(:))]);
set(gca,'XLim',[0 frames]);

subplot(2,4,2);
plot(sqrt(Means1(:,4).^2+Means1(:,5).^2));
hold on;
plot([5.5 5.5],[0 3],'--k','Linewidth',0.5);
plot([0 frames],[mean(sqrt(Means1(1:4,4).^2+Means1(1:4,5).^2)) mean(sqrt(Means1(1:4,4).^2+Means1(1:4,5).^2))],'--r','Linewidth',0.5);
title('C1B meandisp');
set(gca,'YLim',[0 max(Means1(:))]);
set(gca,'XLim',[0 frames]);

subplot(2,4,3);
plot(sqrt(Means1(:,7).^2+Means1(:,8).^2));
hold on;
plot([5.5 5.5],[0 3],'--k','Linewidth',0.5);
plot([0 frames],[mean(sqrt(Means1(1:4,7).^2+Means1(1:4,8).^2)) mean(sqrt(Means1(1:4,7).^2+Means1(1:4,8).^2))],'--r','Linewidth',0.5);
title('Can1 meandisp');
set(gca,'YLim',[0 max(Means1(:))]);
set(gca,'XLim',[0 frames]);

subplot(2,4,4);
plot(sqrt(Means1(:,10).^2+Means1(:,11).^2));
hold on;
plot([5.5 5.5],[0 3],'--k','Linewidth',0.5);
plot([0 frames],[mean(sqrt(Means1(1:4,10).^2+Means1(1:4,11).^2)) mean(sqrt(Means1(1:4,10).^2+Means1(1:4,11).^2))],'--r','Linewidth',0.5);
title('Can1B meandisp');
set(gca,'YLim',[0 max(Means1(:))]);
set(gca,'XLim',[0 frames]);

subplot(2,4,5);
plot(Ang_diff1(:,1));
hold on;
plot([5.5 5.5],[-180 180],'--k','Linewidth',0.5);
plot([0 frames],[90 90],'--r','Linewidth',0.5);
plot([0 frames],[0 0],'--r','Linewidth',0.5);
plot([0 frames],[-90 -90],'--r','Linewidth',0.5);
title('C1A angle');
set(gca,'YLim',[-180 180]);
set(gca,'XLim',[0 frames]);

subplot(2,4,6);
plot(Ang_diff1(:,2));
hold on;
plot([5.5 5.5],[-180 180],'--k','Linewidth',0.5);
plot([0 frames],[90 90],'--r','Linewidth',0.5);
plot([0 frames],[0 0],'--r','Linewidth',0.5);
plot([0 frames],[-90 -90],'--r','Linewidth',0.5);
title('C1B angle');
set(gca,'YLim',[-180 180]);
set(gca,'XLim',[0 frames]);

subplot(2,4,7);
plot(Ang_diff1(:,3));
hold on;
plot([5.5 5.5],[-180 180],'--k','Linewidth',0.5);
plot([0 frames],[90 90],'--r','Linewidth',0.5);
plot([0 frames],[0 0],'--r','Linewidth',0.5);
plot([0 frames],[-90 -90],'--r','Linewidth',0.5);
title('Can1 angle ');
set(gca,'YLim',[-180 180]);
set(gca,'XLim',[0 frames]);

subplot(2,4,8);
plot(Ang_diff1(:,4));
hold on;
plot([5.5 5.5],[-180 180],'--k','Linewidth',0.5);
plot([0 frames],[90 90],'--r','Linewidth',0.5);
plot([0 frames],[0 0],'--r','Linewidth',0.5);
plot([0 frames],[-90 -90],'--r','Linewidth',0.5);
title('Can1B angle ');
set(gca,'YLim',[-180 180]);
set(gca,'XLim',[0 frames]);

print( figh1, '-djpeg90', '-r300', [ ImagesPath, filesep,'Plots',filesep, 'Velocities_angles_can1B.jpg' ] ) ;
close ;



close ;


