function File2 = make_kymo_PIVablation_disp_Clem( Settings2, File2 );

ImagesPath = [File2.ImPath,filesep,'Displacement'] ;
mkdir(ImagesPath);


File2.DispKymosPath=[File2.pathname,filesep,'Disp_Kymos'];
mkdir(File2.DispKymosPath);
Cuttime=Settings2.cuttime;
overlap = Settings2.Overlap ;
filt_size_trac = 2 ;
filt_pow_trac = 0.1;

pad_filter = round( filt_size_trac/(2*(1-overlap)) ) ;
Edges=(Settings2.Resolution*Settings2.Overlap)/2;

DispPath = File2.DispPath ;                  % Location of the traction files.
Name = File2.name ;                % Location of the cropped files.
AllImages = double(TIFFStack( File2.name )) ;
AllImages=AllImages(Edges+1:end-Edges,Edges+1:end-Edges,:);
% AllImages=AllImages(1:File2.TractionSize(1).i*(Settings2.Resolution*(1-Settings2.Overlap)),1:File2.TractionSize(1).i*(Settings2.Resolution*(1-Settings2.Overlap)),:);

SizeX = Settings2.Size_ROIx;
SizeY = Settings2.Size_ROIy;
Step=Settings2.timestep_disp;
nBeadsFile = size(AllImages,3) ;

se=strel('diamond',3*(Settings2.Resolution*(1-Settings2.Overlap)));
se2=strel('diamond',2*(Settings2.Resolution*(1-Settings2.Overlap)));

Mask=imread(File2.Maskname);
Mask(Mask>0)=1;    
Mask=Mask(Edges+1:end-Edges,Edges+1:end-Edges);
Mask=imdilate(Mask,se);
Mask=imerode(Mask,se2);

% Mask=Mask(1:File2.TractionSize(1).i*(Settings2.Resolution*(1-Settings2.Overlap)),1:File2.TractionSize(1).i*(Settings2.Resolution*(1-Settings2.Overlap)));
Mask2=imresize(Mask,[File2.TractionSize(1).i File2.TractionSize(1).j ]);
Mask2=round(Mask2);
figh1=figure;
imagesc(Mask2);
print( figh1, '-djpeg90', '-r300', [ ImagesPath, filesep, 'Mask.jpg' ] ) ;
close;

I=find(Mask2);

%%%%%%%%%%%%%%%%%
Dist=round(bwdist(Mask2));
Mask3=imfill(Mask2,'holes');
Signs=ones(size(Mask3));
Signs(find(Mask3))=-1;
Dist=Dist.*Signs;
Dist(I)=nan;


% Edge=abs(min(Dist(:)))+1;

[s,s2,Krv,Coord,Normal] = curv(Mask3);%compute normal direction
%here i prepare the variables to make the field of normal
%direction vectors
Vecfieldx=zeros(size(Mask2));Front=Vecfieldx;Vecfieldy=Vecfieldx;
Tempx=Vecfieldx;Vecfieldx(:)=NaN;
Tempy=Vecfieldy;Vecfieldy(:)=NaN;
%we introduce the values calculated 4 lines above in the matrices
Tempy(sub2ind(size(Vecfieldx),Coord(2,:),Coord(1,:)))=Normal(:,2);
Tempx(sub2ind(size(Vecfieldx),Coord(2,:),Coord(1,:)))=Normal(:,1);
Vecfieldx=Tempx;Vecfieldy=Tempy;

%first we assign to each point the closest vector
%of the leading edge
Front(sub2ind(size(Mask2),Coord(2,:),Coord(1,:)))=1;
[DF,LF]=bwdist(Front);
[RF,CF]=ind2sub(size(Front),LF);
DF=round(DF);
%     DF(Dista==0)=0;
for ni=1:max2(DF)
    Tempx(find(DF==ni))=Tempx(LF(find(DF==ni)));
    Tempy(find(DF==ni))=Tempy(LF(find(DF==ni)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% then we filter them as a function of the distance
for n=1:round(max2(DF))
    [im,jm]=ind2sub(size(DF),find(DF==n));
    for m=1:size(im,1)
        in=im(m);
        jn=jm(m);
        temp = Tempx(max(in-DF(in,jn),1):min(in+DF(in,jn),size(DF,1)),...
            max(jn-DF(in,jn),1):min(jn+DF(in,jn),size(DF,2)));
        Vecfieldx(in,jn)=nanmean(nanmean(temp));
        temp = Tempy(max(in-DF(in,jn),1):min(in+DF(in,jn),size(DF,1)),...
            max(jn-DF(in,jn),1):min(jn+DF(in,jn),size(DF,2)));
        Vecfieldy(in,jn)=nanmean(nanmean(temp));
    end
end

VecMod=sqrt( Vecfieldx.^2+ Vecfieldy.^2);
Vecfieldy= Vecfieldy./VecMod;
Vecfieldx= Vecfieldx./VecMod;
Vecfieldx(isnan(Vecfieldx))=0;
Vecfieldy(isnan(Vecfieldy))=0;

Radial2=[];
Radial3=[];
Modul=[];
Erro3=[];
Numel=[];
StdErro=[];
Ang2=[];
Ang3=[];

Peaks=[];

MeanDisp=[];
tt=1;
for k = 1:Step:nBeadsFile
    k
    
%     PC=AllImages(:,:,k);
% %     Holes=find(PC<Th);
%     Ex=ones(size(PC));
%     Ex(Holes)=nan;
    
    TempData=load([DispPath,filesep, 'Displacement_', num2str(k), '.dat']);
    
    VelX = reshape(TempData(:,3),sqrt(size(TempData,1)),sqrt(size(TempData,1)));%File.Velocitxx,File.Velocityy);
    VelY = reshape(TempData(:,4),sqrt(size(TempData,1)),sqrt(size(TempData,1)));%File.Velocitxx,File.Velocityy);

    VelX=VelX.*Settings2.PS;
    VelY=VelY.*Settings2.PS;
    %%%%Filter

% VelX_filt = zeros( size(VelX) + 2*pad_filter ) ;         % Displacements in x, filtered.
% VelY_filt = zeros( size(VelY) + 2*pad_filter ) ;         % Displacements in y, filtered.
% 
% % Add the specular image of the borders of the displacement field before filtering.
% VelX_filt( 1+pad_filter:end-pad_filter, 1+pad_filter:end-pad_filter ) = VelX ;
% VelY_filt( 1+pad_filter:end-pad_filter, 1+pad_filter:end-pad_filter ) = VelY ;
% 
% VelX_filt( 1:pad_filter, : ) = VelX_filt( 1+2*pad_filter:-1:2+pad_filter, : ) ;
% VelY_filt( 1:pad_filter, : ) = VelY_filt( 1+2*pad_filter:-1:2+pad_filter, : ) ;
% 
% VelX_filt( end-pad_filter+1:end, : ) = VelX_filt( end-pad_filter-1:-1:end-2*pad_filter, : ) ;
% VelY_filt( end-pad_filter+1:end, : ) = VelY_filt( end-pad_filter-1:-1:end-2*pad_filter, : ) ;
% 
% VelX_filt( :, 1:pad_filter ) = VelX_filt( :, 1+2*pad_filter:-1:2+pad_filter ) ;
% VelY_filt( :, 1:pad_filter ) = VelY_filt( :, 1+2*pad_filter:-1:2+pad_filter ) ;
% 
% VelX_filt( :, end-pad_filter+1:end ) = VelX_filt( :, end-pad_filter-1:-1:end-2*pad_filter ) ;
% VelY_filt( :, end-pad_filter+1:end ) = VelY_filt( :, end-pad_filter-1:-1:end-2*pad_filter ) ;
% 
% VelX_filt = inpaint_nans( VelX ) ;
% VelY_filt = inpaint_nans( VelY ) ;
% 
% VelX_filt = filter_2D_Disp( VelX_filt, filt_size_trac * [ 1/( 1 - overlap), 1/( 1 - overlap) ], filt_pow_trac ) ;
% VelY_filt = filter_2D_Disp( VelY_filt, filt_size_trac * [ 1/( 1 - overlap), 1/( 1 - overlap) ], filt_pow_trac ) ;

% VelX = VelX_filt( 1+pad_filter:end-pad_filter, 1+pad_filter:end-pad_filter ) ;
% VelY = VelY_filt( 1+pad_filter:end-pad_filter, 1+pad_filter:end-pad_filter ) ;

% OUTLIER REMOVAL
% VelX = removeOutliers2D( VelX ) ;
% VelY = removeOutliers2D( VelY ) ;
% 
% VelX = VelX{1} ;
% VelY = VelY{1} ;

%%%%%%%%%%
%  VelX=imresize(VelX,size(Dist),'bilinear');
% VelY=imresize(VelY,size(Dist),'bilinear');
% 
VelX(I)=nan;
VelY(I)=nan;

VelX=VelX-nanmean(VelX(:));
VelY=VelY-nanmean(VelY(:));

    VelXY=sqrt(VelX.^2+VelY.^2);
  
    MeanDisp(tt)=nanmean(VelXY(:));
    
    %calculation of magnitudes of interest of the tractions
    scalX=VelX.*Vecfieldx;
    scalY=VelY.*Vecfieldy;
    roX=VelX.*(-Vecfieldy);
    roY=VelY.*Vecfieldx;
    
    count=1;
    for nu=min(min(Dist)):1:max(max(Dist));
        Pos=find(Dist==nu);
        Radial2(tt,count)=nanmean(abs(scalX(Pos)+scalY(Pos)));
        Radial3(tt,count)=nanmean(scalX(Pos)+scalY(Pos));
        Modul(tt,count)=nanmean(VelXY(Pos));
        Erro3(tt,count)=nanstd(scalX(Pos)+scalY(Pos));
        Numel(tt,count)=numel(Dist(Pos));
        StdErro(tt,count)=Erro3(tt,count)/sqrt(size(find(Pos),1));
        Ang2(tt,count)=nanmean(abs(roX(Pos)+roY(Pos)));
        Ang3(tt,count)=nanmean(roX(Pos)+roY(Pos));
        
        if nu==0;
        Edge=count;
        end;
        
        count=count+1;
    end;
   
    if Edge>5;
        Width=5;
    else;
        Width=Edge-1;
    end;
    
    [M Max_i]=max(abs(Radial3(tt,Edge-Width:Edge)));
    Peaks(tt,1)=Radial3(tt,Edge-Max_i);
    
    [M Max_i]=max(abs(Radial3(tt,Edge:Edge+5)));
    Peaks(tt,2)=Radial3(tt,Edge+Max_i-1);
    
    tt=tt+1;
end;

figh1=figure;
imagesc(Radial3);
hold on;
set(gca,'Clim',[-max(abs(Radial3(:))) max(abs(Radial3(:)))]);
colormap(jet);
h = colorbar;
ylabel(h, 'Radial Disp','FontSize',11,'Rotation',-90,'Position',[2.9, 0 0])
xlabel('Distance (pixels)');
ylabel('Time (s)');
% plot((Edge).*(ones(1,nBeadsFile)),1:nBeadsFile,'k','LineWidth',2);
print( figh1, '-djpeg90', '-r300', [File2.DispKymosPath,'\RadialDisp.jpg' ] ) ;

figh6=figure;
imagesc(Modul);
hold on;
colormap(jet);
h = colorbar;
ylabel(h, 'Disp Magnitude','FontSize',11,'Rotation',-90,'Position',[2.9, 0 0])
xlabel('Distance (pixels)');
ylabel('Time (h)');
% plot((Edge).*(ones(1,nBeadsFile)),1:nBeadsFile,'k','LineWidth',2);
print( figh6, '-djpeg90', '-r300', [File2.DispKymosPath,'\ModDisp.jpg' ] ) ;


figh2=figure;
imagesc(Ang3);
hold on;
set(gca,'Clim',[-max(abs(Radial3(:))) max(abs(Radial3(:)))]);
colormap(jet);
h = colorbar;
ylabel(h, 'Tangential Disp','FontSize',11,'Rotation',-90,'Position',[2.9, 0 0])
xlabel('Distance (pixels)');
ylabel('Time (h)');
% plot((Edge).*(ones(1,nBeadsFile)),1:nBeadsFile,'k','LineWidth',2);
print( figh2, '-djpeg90', '-r300', [File2.DispKymosPath,'\TangentialDisp.jpg' ] ) ;

% figh3=figure;
% plot(Radial3(1,:));
% set(gca,'Ylim',[-max(abs(Radial3(1,:)))-0.5 max(abs(Radial3(1,:)))+0.5]);
% % hold on;
% % plot((Edge).*(ones(1,7)),-3:3,'k');
% print( figh3, '-djpeg90', '-r300', [File2.DispKymosPath,'\RadialDisp_cuttime.jpg' ] ) ;

figh4=figure;
plot(MeanDisp);
print( figh4, '-djpeg90', '-r300', [File2.DispKymosPath,'\MeanDisp.jpg' ] ) ;

figh5=figure;
plot(Radial3(end,:));
set(gca,'Ylim',[-max(abs(Radial3(end,:)))-0.5 max(abs(Radial3(end,:)))+0.5]);
% hold on;
% plot((Edge).*(ones(1,7)),-3:3,'k');
print( figh5, '-djpeg90', '-r300', [File2.DispKymosPath,'\RadialDisp_end.jpg' ] ) ;

figh8=figure;
plot(abs(Peaks(:,1)));
hold on;
plot(abs(Peaks(:,2)));
ylabel('Recoil (um)');
xlabel('Time (s)');

legend({'Internal recoil', 'External recoil'},'Location','best');
print( figh8, '-djpeg90', '-r300', [File2.DispKymosPath,'\Internal_external_recoilds.jpg' ] ) ;

figh9=figure;
plot(-Peaks(:,1)+Peaks(:,2));
ylabel('Recoil (um)');
xlabel('Time (s)');
print( figh9, '-djpeg90', '-r300', [File2.DispKymosPath,'\Total displacement.jpg' ] ) ;

close all;

save([File2.DispKymosPath,'\Radial2.mat'],'Radial2');
save([File2.DispKymosPath,'\Radial3.mat'],'Radial3');
save([File2.DispKymosPath,'\Modul.mat'],'Modul');
save([File2.DispKymosPath,'\Erro3.mat'],'Erro3');
save([File2.DispKymosPath,'\Numel.mat'],'Numel');
save([File2.DispKymosPath,'\StdErro.mat'],'StdErro');
save([File2.DispKymosPath,'\Ang2.mat'],'Ang2');
save([File2.DispKymosPath,'\Ang3.mat'],'Ang3');
save([File2.DispKymosPath,'\MeanDisp.mat'],'MeanDisp');
save([File2.DispKymosPath,'\Peaks.mat'],'Peaks');
save([File2.DispKymosPath,'\Edge.mat'],'Edge');
