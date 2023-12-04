clear all;
close all;

Path=['D:\Ablations_buds\'];
cd(Path);

Npos=13;
positions=[1:5 1:8];
% positions=[1];


DilCancer = strel('octagon',75);
DilCafs = strel('octagon',75);

load('Pos_PIV.mat'); 
% count=1;
% for Pos=1:Npos;
%     folder_names{count}=uigetdir;    
%     count=count+1;
% end

for Pos=1:Npos;

    cd(folder_names{Pos});

load File2;
% load Settings2;
Settings2.cuttime=5;

ImagesPath = File2.ImPath ;
MaskPath = [File2.pathname, filesep,'Boundary'];

mkdir([MaskPath]);
mkdir([ImagesPath,filesep,'Polar_plots']);
mkdir([ImagesPath,filesep,'Vel_Vectors']);
mkdir([ImagesPath,filesep,'Inst_Polar_plots']);
mkdir([ImagesPath,filesep,'InstVel_Vectors']);
mkdir([ImagesPath,filesep,'Plots']);

File2.CutMask_dir=[Path,filesep,CutMask_dir,filesep,'s',num2str(positions(Pos)),filesep,'Mask.tif'];

Settings2.Size_ROIx=512 ;         % This is the size of the ROIs (less than Settings.SizeX and multiple of Settings.resolution).
Settings2.Size_ROIy=512 ;         % This is the size of the ROIs (less than Settings.SizeX and multiple of Settings.resolution).

SizeX = Settings2.Size_ROIx;
SizeY = Settings2.Size_ROIy;

%%Define rois to calculate recoils
Rois_tissue_cuts(File2, Settings2, Path, positions(Pos), DilCancer, DilCafs);

%%Quantify displacement magnitude and angle difference respect to cut
Calculate_ang_diff_tissue_disp_360(File2, Settings2);

%%Plot displacement magnitude and angle difference respect to cut
Plot_tissue_recoils_disp (File2);

clear File2;
close all;
end;




