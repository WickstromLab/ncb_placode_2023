function Rois_tissue_cuts(File2, Settings2, DilCancer, DilCafs);

MaskPath = [File2.pathname, filesep,'Boundary'];

Del= strel('diamond',2);

se=strel('diamond',1);
hg = fspecial('gaussian', 2, 2); %%parameter for gaussian filter

Edges=(Settings2.Size_ROIx-File2.TractionSize(1).i*(Settings2.Resolution*(1-Settings2.Overlap)))/2;

Cut_mask=imread([File2.CutMask_dir]);
Cut_mask=double(Cut_mask(Edges+1:end-Edges,Edges+1:end-Edges));
Cut_mask(find(Cut_mask))=1;
Cent=regionprops(Cut_mask,'Centroid');
Cent=round(Cent.Centroid);
save([MaskPath,filesep,'Cent.mat'],'Cent');
save([MaskPath,filesep,'Cut_mask.mat'],'Cut_mask');

CutBorder=bwboundaries(Cut_mask);

CanMask=round(imdilate (Cut_mask,DilCancer));
DelMask=round(imdilate (Cut_mask,Del));

Indices=find(DelMask);
save([MaskPath,filesep,'Indices.mat'],'Indices');

load([MaskPath,filesep,'preMask.mat']);
Mask=preMask(Edges+1:end-Edges,Edges+1:end-Edges);
save([MaskPath,filesep,'Mask.mat'],'Mask');

%%%%Caf mask for quantifications
CafMask=round(imdilate (Cut_mask,DilCafs));

Cut_mask2=Cut_mask;

for ee=1:10;
    Temp=imerode(Cut_mask2,se);
    CC=bwconncomp(Temp);
    if CC.NumObjects~=1; break
    else
        Cut_mask2=Temp;
    end
end

[Y X]=find(Cut_mask2);
% p = polyfit(X,Y,1);
p = linortfit2(X,Y);
pTheta=p;
Theta = atand(pTheta(1));
if Theta<0;
    Theta=Theta + 180;
end;

ThetaCan(1,1:2)=[cos(deg2rad(Theta)) sin(deg2rad(Theta))];
ThetaCan(2,1:2)=[-cos(deg2rad(Theta)) -sin(deg2rad(Theta))];
ThetaCAF(1,1:2)=[sin(deg2rad(Theta)) -cos(deg2rad(Theta))];
ThetaCAF(2,1:2)=[-sin(deg2rad(Theta)) cos(deg2rad(Theta))];

Cent1=regionprops(Mask,'Centroid');
Cent1=round(Cent1.Centroid);

dThetaCan=[];
VectMask1=Cent-Cent1;
VectMask1=VectMask1./sqrt(VectMask1(1).^2+VectMask1(2).^2);
dThetaCan(1)=rad2deg(acos(dot(VectMask1,ThetaCan(1,1:2))));
dThetaCan(2)=rad2deg(acos(dot(VectMask1,ThetaCan(2,1:2))));
Theta1=ThetaCan(find(dThetaCan==min(dThetaCan)),1:2);

% save([MaskPath,filesep,'Mask1.mat'],'Mask1');
save([MaskPath,filesep,'Theta1.mat'],'Theta1');
save([MaskPath,filesep,'Cent1.mat'],'Cent1');

if abs(p(1))==inf;
ParMask=zeros(size(Mask));
ParMask(:,round(mean(X)))=1;
else
xx=1:0.001:size(Mask,1);
yy=round(polyval(p,xx));
xx=round(xx);
yy_ind=find(yy>0 & yy<=(size(Mask,1)));

ParMask=zeros(size(Mask));
for kk=1:numel(yy_ind);
    ParMask((yy(yy_ind(kk))),(xx(yy_ind(kk))))=1;
end;
end;

perp_Theta=linortfit2(X,size(Mask,1)-Y);
perp_Theta = atand(perp_Theta(1));

Temp=zeros(size(ParMask));
Temp((ParMask+double(Cut_mask))==2)=1;
[row col]=find(bwmorph(Temp,'endpoints'));

Distances=sqrt((row-Cent1(2)).^2+(col-Cent1(1)).^2);
End=find(Distances==max(Distances));

se_perp1=strel('line',size(Mask,1)*2,perp_Theta-90);
se_perp2=strel('line',size(Mask,1)*2,perp_Theta+90);
Temp_perp=zeros(size(Mask));
Temp_perp(row(End),col(End))=1;
Temp_perp=imdilate(Temp_perp,se_perp1);
Temp_perp=imdilate(Temp_perp,se_perp2);

ParMask=ParMask+Temp_perp;
ParMask(ParMask==2)=1;

ParMaskdil=~(imdilate(ParMask,se));
CC=bwconncomp(ParMaskdil);
B0=zeros(size(Mask));
for cc=1:numel(CC.PixelIdxList);
    B0(CC.PixelIdxList{cc}) = -cc; %create a mask with the biggest region
end;

Caf_mask=B0.*double(CafMask);
Caf_mask(Indices)=0;

B0=zeros(size(Mask,1),size(Mask,2),4);
for cc=1:4;
    TempInd=find(Caf_mask==-cc);
    TempB0=zeros(size(Mask));
    TempB0(TempInd)=1;
    B0(:,:,cc)=TempB0;
end;
    
    C1=(B0+Mask);
    for cc=1:4;
        Max1(cc)=max(max(C1(:,:,cc)));
    end;
    C1Ind=find(Max1==2);
    
    if isempty (C1Ind);
        
        C1=(B0+Cut_mask);
    for cc=1:numel(CC.PixelIdxList);
        Max1(cc)=numel(find(C1(:,:,cc)==2));
    end;
    [B I]=sort(Max1);
    C1A=B0(:,:,I(end));
    C1B=B0(:,:,I(end-1));

    else;
    C1A=B0(:,:,C1Ind(1));
    C1B=B0(:,:,C1Ind(2));    
    end;
    

    
    C1A(find(Mask))=0;
    C1B(find(Mask))=0;
    
        CC=bwconncomp(C1A);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [biggest,idx] = max(numPixels);
    C1A=zeros(size(Mask));
    C1A(CC.PixelIdxList{idx}) = 1;
        
    CC=bwconncomp(C1B);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [biggest,idx] = max(numPixels);
    C1B=zeros(size(Mask));
    C1B(CC.PixelIdxList{idx}) = 1;
    
    CentC1A=regionprops(C1A,'Centroid');
    CentC1B=regionprops(C1B,'Centroid');
    
    CentC1A=round(CentC1A.Centroid);
    CentC1B=round(CentC1B.Centroid);
    
    C1A_i=find(C1A==1);
    C1B_i=find(C1B==1);
    
    dThetaCAF=[];
    VectC1A=Cent-CentC1A;
    VectC1A=VectC1A./sqrt(VectC1A(1).^2+VectC1A(2).^2);
    dThetaCAF(1)=rad2deg(acos(dot(VectC1A,ThetaCAF(1,1:2))));
    dThetaCAF(2)=rad2deg(acos(dot(VectC1A,ThetaCAF(2,1:2))));
    ThetaC1A=ThetaCAF(find(dThetaCAF==min(dThetaCAF)),1:2);
    ThetaC1B=ThetaCAF(find(dThetaCAF==max(dThetaCAF)),1:2);
    
    save([MaskPath,filesep,'C1A.mat'],'C1A');
    save([MaskPath,filesep,'C1B.mat'],'C1B');
    
    save([MaskPath,filesep,'CentC1A.mat'],'CentC1A');
    save([MaskPath,filesep,'CentC1B.mat'],'CentC1B');
    
    save([MaskPath,filesep,'ThetaC1A.mat'],'ThetaC1A');
    save([MaskPath,filesep,'ThetaC1B.mat'],'ThetaC1B');
    
    
    %%%%Cancer cell mask for quantifications
MM=bwconncomp(Mask);

MMs=zeros([size(Mask) numel(MM.PixelIdxList)]);
Intersect=zeros([size(Mask) numel(MM.PixelIdxList)]);
II=[];
for mm=1:numel(MM.PixelIdxList);
    MM_ind=MM.PixelIdxList{mm};
    Temp=zeros(size(Mask));
    Temp(MM_ind)=1;
    MMs(:,:,mm)=Temp;
    Intersect(:,:,mm)=Temp+double(CanMask);
    II(mm)=max(max(Intersect(:,:,mm)+ParMask));
end;

M1=find(II==3);

if isempty (M1);
Mask1_ind=find(Mask==1);
    
else;
Mask1_ind=find(Intersect(:,:,M1(1))==2);
    
end;

Mask1=zeros(size(Mask));
Mask1(Mask1_ind)=1;

CC=bwconncomp(Mask1);
numPixels = cellfun(@numel,CC.PixelIdxList);
[biggest,idx] = max(numPixels);
Mask1=zeros(size(Mask));
Mask1(CC.PixelIdxList{idx}) = 1;

Cent1=regionprops(Mask1,'Centroid');
Cent1=round(Cent1.Centroid);

dThetaCan=[];
VectMask1=Cent-Cent1;
VectMask1=VectMask1./sqrt(VectMask1(1).^2+VectMask1(2).^2);
dThetaCan(1)=rad2deg(acos(dot(VectMask1,ThetaCan(1,1:2))));
dThetaCan(2)=rad2deg(acos(dot(VectMask1,ThetaCan(2,1:2))));
Theta1=ThetaCan(find(dThetaCan==min(dThetaCan)),1:2);

save([MaskPath,filesep,'Mask1.mat'],'Mask1');
save([MaskPath,filesep,'Theta1.mat'],'Theta1');
save([MaskPath,filesep,'Cent1.mat'],'Cent1');

[Row Col]=find(MMs(:,:,1));
Dist=sqrt((Row-Cent1(2)).^2+(Col-Cent1(1)).^2);
Mask1B=zeros(size(Mask1));
[M I]=max(Dist);
Mask1B(Row(I),Col(I))=1;

Sdif=nan(1,1);
for mm=1:1:100;
    Mask1B=round(imdilate(Mask1B,se));
    Mask1B(find(~MMs(:,:,1)))=0;
    Sdif=numel(find(Mask1))-numel(find(Mask1B));
    if Sdif<=0;
        break
    end;
end;

Mask1B_ind=find(Mask1B);
Cent1B=regionprops(Mask1B,'Centroid');
Cent1B=round(Cent1B.Centroid);

save([MaskPath,filesep,'Mask1B.mat'],'Mask1B');
save([MaskPath,filesep,'Cent1B.mat'],'Cent1B');
