function Calculate_ang_diff_tissue_disp_360(File2, Settings2 );

VelPath = File2.DispPath ;                  % Location of the traction files.
DataPath = [File2.pathname, filesep,'Quantifications'];
MaskPath = [File2.pathname, filesep,'Boundary'];
ImagesPath = File2.ImPath ;
mkdir([ImagesPath,filesep,'Vel_Vectors']);

hg = fspecial('gaussian', 2, 2); %%parameter for gaussian filter
se=strel('diamond',1);

spacing = 12;

scale_factor = 6;
scale_factor2 = 15;
scale_factor3 = 25;

scale_len = 3;

mkdir(DataPath);

AllCancerImages = double(TIFFStack( File2.name )) ;

frames = File2.NFiles.Beads-1 ;
Edges=(Settings2.Size_ROIx-File2.TractionSize(1).i*(Settings2.Resolution*(1-Settings2.Overlap)))/2;


load([MaskPath,filesep,'Cut_mask.mat']);
load([MaskPath,filesep,'Mask.mat']);
load([MaskPath,filesep,'Mask1.mat']);
load([MaskPath,filesep,'Mask1B.mat']);
load([MaskPath,filesep,'Theta1.mat']);
load([MaskPath,filesep,'Cent1.mat']);
load([MaskPath,filesep,'Cent1B.mat']);
load([MaskPath,filesep,'C1A.mat']);
load([MaskPath,filesep,'C1B.mat']);
load([MaskPath,filesep,'CentC1A.mat']);
load([MaskPath,filesep,'CentC1B.mat']);
load([MaskPath,filesep,'ThetaC1A.mat']);
load([MaskPath,filesep,'ThetaC1B.mat']);


Ref1=atan2(ThetaC1A(2), ThetaC1A(1)) - atan2(Theta1(2), Theta1(1));

Mask_LR=imresize(Mask,[File2.TractionSize(1).i File2.TractionSize(1).j]);
Cut_mask_LR=imresize(Cut_mask,[File2.TractionSize(1).i File2.TractionSize(1).j]);
Mask1_LR=imresize(Mask1,[File2.TractionSize(1).i File2.TractionSize(1).j]);
Mask1B_LR=imresize(Mask1B,[File2.TractionSize(1).i File2.TractionSize(1).j]);
C1A_LR=imresize(C1A,[File2.TractionSize(1).i File2.TractionSize(1).j]);
C1B_LR=imresize(C1B,[File2.TractionSize(1).i File2.TractionSize(1).j]);

Mask_LR=round(Mask_LR);
Cut_mask_LR=round(Cut_mask_LR);
Mask1_LR=round(Mask1_LR);
Mask1B_LR=round(Mask1B_LR);
C1A_LR=round(C1A_LR);
C1B_LR=round(C1B_LR);

Cut_mask_LR=imdilate(Cut_mask_LR,se);
Cut_mask_LR=round(Cut_mask_LR);

Indices=find(Cut_mask_LR==1);
Mask1_ind=find(Mask1_LR==1);
Mask1B_ind=find(Mask1B_LR==1);
C1A_i=find(C1A_LR==1);
C1B_i=find(C1B_LR==1);

Means1=nan(frames,12);
Ang_diff1=nan(frames,8);
C1A_dat=nan(numel(C1A_i),frames,2);
C1B_dat=nan(numel(C1B_i),frames,2);
Can1_dat=nan(numel(Mask1_ind),frames,2);
Can1B_dat=nan(numel(Mask1B_ind),frames,2);

for k=1:frames;
    k
    %%%%CAFs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    TempPC = AllCancerImages(:,:,k);
    Im=TempPC;
    Im = imfilter(Im,hg,'replicate'); %% nice trick to remove noise
    
    Im=Im (Edges+1:end-Edges,Edges+1:end-Edges);
  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Loading Cafs velocities
    TempData=load([VelPath,filesep, 'Displacement_', num2str(k), '.dat']);
    
    VelX = reshape(TempData(:,3),sqrt(size(TempData,1)),sqrt(size(TempData,1)));%File.Velocitxx,File.Velocityy);
    VelY = reshape(TempData(:,4),sqrt(size(TempData,1)),sqrt(size(TempData,1)));%File.Velocitxx,File.Velocityy);
    VelX = VelX-nanmean(VelX(:));
    VelY = VelY-nanmean(VelY(:));
    VelXY=sqrt(VelX.^2+VelY.^2);
    
    VelX_HR=imresize(VelX,size(Im));
    VelY_HR=imresize(VelY,size(Im));
    VelXY_HR=imresize(VelXY,size(Im));
        
    VelX(Indices)=nan;
    VelY(Indices)=nan;
    VelXY(Indices)=nan;
    
    %Loading cancer velocities
    %Cafs 1 quantification;
    MeanC1AX=nanmean(VelX(C1A_i));
    MeanC1AY=nanmean(VelY(C1A_i));
    MeanC1AXY=nanmean(VelXY(C1A_i));
    
    Means1(k,1)=MeanC1AX;
    Means1(k,2)=MeanC1AY;
    Means1(k,3)=MeanC1AXY;
    
    %Angles
    C1AX=VelX(C1A_i);
    C1AY=VelY(C1A_i);
    C1AXY=VelXY(C1A_i);
    
    TempA=Angle_2D_sign([C1AX C1AY], Theta1);
    C1AMean=Angle_2D_sign([MeanC1AX MeanC1AY],Theta1);
    
    if Ref1>0;
    TempA = -TempA;  
    C1AMean=-C1AMean;
    end;
    
    C1A_dat(:,k,1) = TempA;
    C1A_dat(:,k,2)=sqrt(C1AX.^2+C1AY.^2);
    
    Tempcos=cos(TempA);
    Tempsin=sin(TempA);
    TempAng=rad2deg(atan2(nanmean(Tempsin),nanmean(Tempcos)));
        
    Ang_diff1(k,1)=TempAng;
    Ang_diff1(k,5)=rad2deg(C1AMean);
    
    %Cafs C1B quantification
    
    MeanC1BX=nanmean(VelX(C1B_i));
    MeanC1BY=nanmean(VelY(C1B_i));
    MeanC1BXY=nanmean(VelXY(C1B_i));
    
    Means1(k,4)=MeanC1BX;
    Means1(k,5)=MeanC1BY;
    Means1(k,6)=MeanC1BXY;
    
    C1BX=VelX(C1B_i);
    C1BY=VelY(C1B_i);
    C1BXY=VelXY(C1B_i);
    
     TempA=-Angle_2D_sign([C1BX C1BY], Theta1);
    C1BMean=-Angle_2D_sign([MeanC1BX MeanC1BY],Theta1);
    
    if Ref1>0;
    TempA = -TempA;  
    C1BMean=-C1BMean;
    end;
    
    C1B_dat(:,k,1) = TempA;
    C1B_dat(:,k,2)=sqrt(C1BX.^2+C1BY.^2);
     
    Tempcos=cos(TempA);
    Tempsin=sin(TempA);
    TempAng=rad2deg(atan2(nanmean(Tempsin),nanmean(Tempcos)));
    
    Ang_diff1(k,2)=TempAng;
    Ang_diff1(k,6)=rad2deg(C1BMean);
    
    %Cancer cell quantification
    
    MeanX1=nanmean(VelX(Mask1_ind));
    MeanY1=nanmean(VelY(Mask1_ind));
    MeanXY1=nanmean(VelXY(Mask1_ind));
       
    Means1(k,7)=MeanX1;
    Means1(k,8)=MeanY1;
    Means1(k,9)=MeanXY1;
    
    Can1X=VelX(Mask1_ind);
    Can1Y=VelY(Mask1_ind);
    Can1XY=VelXY(Mask1_ind);
    
     TempA=Angle_2D_sign([Can1X Can1Y], Theta1);
    Can1Mean=Angle_2D_sign([MeanX1 MeanY1], Theta1);
    
    Can1_dat(:,k,1)=TempA;
    Can1_dat(:,k,2)=sqrt(Can1X.^2+Can1Y.^2);
        
    Tempcos=cos(TempA);
    Tempsin=sin(TempA);
    TempAng=rad2deg(atan2(nanmean(Tempsin),nanmean(Tempcos)));
    
    Ang_diff1(k,3)=TempAng;
    Ang_diff1(k,7)=rad2deg(Can1Mean);
   
    % Cancer cells control
    MeanX1B=nanmean(VelX(Mask1B_ind));
    MeanY1B=nanmean(VelY(Mask1B_ind));
    MeanXY1B=nanmean(VelXY(Mask1B_ind));
       
    Means1(k,10)=MeanX1B;
    Means1(k,11)=MeanY1B;
    Means1(k,12)=MeanXY1B;
    
    Can1BX=VelX(Mask1B_ind);
    Can1BY=VelY(Mask1B_ind);
    Can1BXY=VelXY(Mask1B_ind);
    
     TempA=Angle_2D_sign([Can1BX Can1BY], Theta1);
    Can1BMean=Angle_2D_sign([MeanX1B MeanY1B], Theta1);
    
    Can1B_dat(:,k,1)=TempA;
    Can1B_dat(:,k,2)=sqrt(Can1BX.^2+Can1BY.^2);
        
    Tempcos=cos(TempA);
    Tempsin=sin(TempA);
    TempAng=rad2deg(atan2(nanmean(Tempsin),nanmean(Tempcos)));
    
    Ang_diff1(k,4)=TempAng;
    Ang_diff1(k,8)=rad2deg(Can1BMean);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%Plot Can 1    
    x=(1:size(Im,1));
    y=(1:size(Im,2));
    
    [X,Y] = meshgrid(x,y);
    
    figh1 = figure ;
    set( figh1, 'visible', 'off' ) ;
    
    imshow(imadjust(uint16(Im))); axis image;
    %             colormap(gray);
    axis 'manual';
    hold on;
    
    xx=(1:spacing:size(Mask,1));
    yy=(1:spacing:size(Mask,2));
    
    u = zeros(size(xx));
    v = zeros(size(yy));
    
    [U,V] = meshgrid(u,v);
    
    %%%Vel CAFs
    
    for i=1:length(xx)
        for j=1:length(yy)
            U(i,j) = VelX_HR(xx(i),yy(j));
            V(i,j) = VelY_HR(xx(i),yy(j));
        end
    end
    
    AA=zeros(size(Im));
    AA(xx,yy)=1;
    
    
    VelX_HR=VelX_HR.*Cut_mask;
    VelY_HR=VelY_HR.*Cut_mask;   
    
    XX=VelX_HR.*AA*scale_factor;
    YY=VelY_HR.*AA*scale_factor;
    
    U = U.*scale_factor;
    V = V.*scale_factor;
    
    %     quiver(X,Y,XX,YY,0,'Linewidth',1,'Color','k');
    quiver(xx,yy,U,V,0,'Linewidth',1,'Color','y');
    
    PP=bwboundaries(C1A);
    plot(PP{1}(:,2),PP{1}(:,1),'--r','linewidth',2);
    PP=bwboundaries(C1B);
    plot(PP{1}(:,2),PP{1}(:,1),'--r','linewidth',2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    quiver(CentC1A(1),CentC1A(2),MeanC1AX*scale_factor2,MeanC1AY*scale_factor2,0,'Linewidth',3,'Color','r','MaxHeadSize',5);
    quiver(CentC1B(1),CentC1B(2),MeanC1BX*scale_factor2,MeanC1BY*scale_factor2,0,'Linewidth',3,'Color','r','MaxHeadSize',5);
   
    Centroid=regionprops(Cut_mask,'Centroid');
%     quiver(Centroid.Centroid(1),Centroid.Centroid(2),ThetaC1A(1)*scale_factor3,ThetaC1A(2)*scale_factor3,0,'Linewidth',2,'Color','b','MaxHeadSize',5);
%     quiver(Centroid.Centroid(1),Centroid.Centroid(2),ThetaC1B(1)*scale_factor3,ThetaC1B(2)*scale_factor3,0,'Linewidth',2,'Color','m','MaxHeadSize',5);
%     quiver(Centroid.Centroid(1),Centroid.Centroid(2),Theta1(1)*scale_factor3,Theta1(2)*scale_factor3,0,'Linewidth',2,'Color','y','MaxHeadSize',5);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PP=bwboundaries(Cut_mask);
    plot(PP{1}(:,2),PP{1}(:,1),'g','linewidth',2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PP=bwboundaries(Mask);
    plot(PP{1}(:,2),PP{1}(:,1),'r','linewidth',2);
    if numel(PP)>1;
    plot(PP{2}(:,2),PP{2}(:,1),'r','linewidth',2);
    end;
    PP=bwboundaries(Mask1);
    plot(PP{1}(:,2),PP{1}(:,1),'--r','linewidth',2);
    PP=bwboundaries(Mask1B);
    plot(PP{1}(:,2),PP{1}(:,1),'--r','linewidth',2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    quiver(Cent1(1),Cent1(2),MeanX1*scale_factor3,MeanY1*scale_factor3,0,'Linewidth',3,'Color','r','MaxHeadSize',5);
    quiver(Cent1B(1),Cent1B(2),MeanX1B*scale_factor3,MeanY1B*scale_factor3,0,'Linewidth',3,'Color','r','MaxHeadSize',5);
   
    %make a scale arrow
    %     quiver(max(x)-max(x)*0.12,max(y)-max(y)*0.05,scale_len*scale_factor,0,0,'Linewidth',1,'Color','y');
    set(gcf,'Color','k');
    set(gcf, 'InvertHardcopy', 'off');
    axis off
    %
    print( figh1, '-djpeg90', '-r300', [ ImagesPath, filesep,'Vel_Vectors',filesep, 'Vel_1_w',num2str(scale_len),'PaSB_v2b_t',num2str(k),'.jpg' ] ) ;
    
    close ;
    
    
end;

save([DataPath,filesep,'Means1.mat'],'Means1');
save([DataPath,filesep,'Ang_diff1.mat'],'Ang_diff1');
save([DataPath,filesep,'C1A_dat.mat'],'C1A_dat');
save([DataPath,filesep,'C1B_dat.mat'],'C1B_dat');
save([DataPath,filesep,'Can1_dat.mat'],'Can1_dat');
save([DataPath,filesep,'Can1B_dat.mat'],'Can1B_dat');




