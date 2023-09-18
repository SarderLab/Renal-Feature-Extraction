function v=feature_extraction_inner_glomerulus(q,segmented_gloms,image_dir,min_object_size)

% Store features in this vector
v=zeros(1,315);
% Search a radius of 2 pixels outside nuclear boundary for boundary
% adjacency features
nucpixradius=2;
% Read in component segmentation and image
composite=imread(fullfile(segmented_gloms(q).folder,segmented_gloms(q).name))>0;
I=imread(fullfile(image_dir(q).folder,image_dir(q).name));

mes_mask=composite(:,:,1);
white_mask=composite(:,:,2);
nuc_mask=composite(:,:,3);

boundary_mask=mes_mask|white_mask|nuc_mask;

% Generate inverted glomerular distance transform for glomerular distance
% transform feature generation
gdist=bwdist(~boundary_mask);
gdist=(-1*(gdist))+max(gdist(:));
gdist(~boundary_mask)=0;

% Determine glomerular area
gArea=sum(boundary_mask(:));
 % Determine glomerular boundary
gOutline=bwperim(boundary_mask);

% Determine glomerular centroid
[r,c]=find(boundary_mask);
rMean=round(mean(r));
cMean=round(mean(c));



gdist2=zeros(size(boundary_mask));

gdist2(~boundary_mask)=1;
gdist2=bwdist(gdist2);
gdist2=(gdist2-max(max(gdist2)))*-1;
gdist2(~boundary_mask)=0;

nuc_dist_bound=double(gdist2.*double(nuc_mask));
lum_dist_bound=double(gdist2.*double(white_mask));
mes_dist_bound=double(gdist2.*double(mes_mask));
% figure(1),imshow(I)
% figure(2),imagesc(mes_dist_bound),axis image,xticks([]),yticks([]),colormap('jet')
lv=quantile(nonzeros(lum_dist_bound(:)),[.1:.1:1]);
nv=quantile(nonzeros(nuc_dist_bound(:)),[.1:.1:1]);
mv=quantile(nonzeros(mes_dist_bound(:)),[.1:.1:1]);
% nuc_dist_display=I;
% for i=1:10
%     r_1=gdist2>=floor(mv(i))-3&gdist2<=ceil(mv(i))+3;
% 
%     if i==41
%         nuc_dist_display=imoverlay(nuc_dist_display,r_1,[0,1,0]);
%     else
%     	nuc_dist_display=imoverlay(nuc_dist_display,r_1,[0,1,0]);
%     end
% end

% figure(3),imshow(nuc_dist_display),pause,return

%figure,imshow(nuc_mask),hold on,pause

[y,x]=find(nuc_mask);

[theta,rho]=cart2pol(x-rMean,y-cMean);

N_n=histcounts(rho,[0:100:1000,1300]);
T_n=histcounts(theta,[-pi:(2*pi/20):pi]);

[y,x]=find(white_mask);

[theta,rho]=cart2pol(x-rMean,y-cMean);

N_l=histcounts(rho,[0:100:1000,1300]);


[y,x]=find(mes_mask);

[theta,rho]=cart2pol(x-rMean,y-cMean);

N_m=histcounts(rho,[0:100:1000,1300]);

% Get PAS+ precursor region from HSV transform to evaluate PAS+
% intra-structural distance features
% HSV=rgb2hsv(I);
% saturation=imadjust(im2double(HSV(:,:,2)),[],[],2);
% mes=imbinarize(saturation,graythresh(saturation));

% mes=mesInt>0.5;
mes=mes_mask;
mes=bwareaopen(mes,min_object_size);
% mes(~boundary_mask)=0;
% I2=I;
% I2(repmat(mes,[1,1,3]))=0;
% I(~repmat(mes,[1,1,3]))=0;
% I2(~repmat(boundary_mask,[1,1,3]))=0;
% figure(1),subplot(121),imshow(I)
% subplot(122),imshow(I2),pause,return
%PAS+ distance transform image
mes2=mes;
mes2=~bwareaopen(~mes2,min_object_size);

mdt=bwdist(~mes2);

%Manually selected distance transform cuts for PAS+ component
m_ext1=mdt>0&mdt<=10;
m_ext2=mdt>10&mdt<=20;
m_ext3=mdt>20&mdt<1000;

% Get histogram data from each various glomerular component distance
% transform
edges=[1:2:80,2000];
N1=histcounts(mdt(mdt(:)>0),edges);
ldt=bwdist(~white_mask);
edges=[1:1:60,2000];
N2=histcounts(ldt(ldt(:)>0),edges);

ndt=bwdist(~nuc_mask);
edges=[1:1:20,2000];
N3=histcounts(ndt(ndt(:)>0),edges);

edges=[2:25:600,20000];
N4=histcounts(gdist(gdist(:)>0),edges);

% Create grayscale representation of image to determine textural features
grayIm=rgb2gray(I);
grayIm(~mes_mask)=NaN;

%    Determine textural and compartment containment features
[ratiosM,s1,mes_num]=getCompRatios(composite,grayIm,min_object_size);

% Re-orient the segmentation channels so that the function 'getCompRatios'
% knows which segmentation is the primary compartment to be examined
composite=cat(3,white_mask,mes_mask,nuc_mask);
composite(~repmat(boundary_mask,[1,1,3]))=0;

% Repeat the steps above for luminal compartments
grayIm=rgb2gray(I);
grayIm(~white_mask)=NaN;

[ratiosL,s2,lum_num]=getCompRatios(composite,grayIm,min_object_size);

% Re-orient the segmentation channels so that the function 'getNucRatios'
% knows which segmentation is the primary compartment to be examined
composite=cat(3,nuc_mask,white_mask,mes_mask);
composite(~repmat(boundary_mask,[1,1,3]))=0;
grayIm=rgb2gray(I);
grayIm(~nuc_mask)=NaN;

% Get nuclear ratios
[ratiosN,s3,nuc_num]=getNucRatios(composite,nucpixradius,grayIm);

% Get distance features between the glomerular periphery, center, and
% between compartments
distsN=getCompDists(nuc_mask,gOutline,[rMean,cMean]);
distsM=getCompDists(mes_mask,gOutline,[rMean,cMean]);
distsL=getCompDists(white_mask,gOutline,[rMean,cMean]);

%Calculate lumen compartmentalization features
v(1,1:3)=mean(ratiosL(:,1:3));
v(1,4)=sum(ratiosL(:,4));
v(1,5)=mean(ratiosL(:,4));
v(1,6)=median(ratiosL(:,4));

%Unpack luminal texture features
v(1,7:10)=[s1(1,1).Contrast,s1(1,2).Correlation,s1(1,3).Energy,s1(1,4).Homogeneity];

%Calculate PAS+ morphological features
v(1,11:13)=mean(ratiosM(:,1:3));
v(1,12)=mean(ratiosM(:,2));
v(1,13)=mean(ratiosM(:,3));
v(1,14)=sum(ratiosM(:,4));
v(1,15)=mean(ratiosM(:,4));
v(1,16)=median(ratiosM(:,4));

% Calculate PAS+ texture
v(1,17:20)=[s2(1,1).Contrast,s2(1,2).Correlation,s2(1,3).Energy,s2(1,4).Homogeneity];

%Calculate nuclear comparmentalization features
v(1,21)=mean(ratiosN(:,1));
v(1,22)=mean(ratiosN(:,2));
v(1,23)=mean(ratiosN(:,3));

v(1,24)=sum(ratiosN(:,4));
v(1,25)=mean(ratiosN(:,4));
v(1,26)=mode(ratiosN(:,4));
% Nuclear texture
v(1,27:30)=[s3(1,1).Contrast,s3(1,2).Correlation,s3(1,3).Energy,s3(1,4).Homogeneity];

%Calculate luminal distance features
v(1,31:37)=mean(distsL);

%Calculate mesangial distance features
v(1,38:44)=mean(distsM);
%Calculate nuclear distance features 
v(1,45:51)=mean(distsN);

%Glomerular area
v(1,52)=gArea;
% Object numbers in each component
v(1,53)=mes_num;
v(1,54)=lum_num;
v(1,55)=nuc_num;
% Hand-selected distance transform cut-associated features
v(1,56)=sum(sum(m_ext1));
v(1,57)=sum(sum(m_ext2));
v(1,58)=sum(sum(m_ext3));
if sum(sum(m_ext2))==0
    v(1,59)=0;
else
    v(1,59)=max(max(mdt(m_ext2)));
end
v(1,60)=max(max(bwlabel(m_ext1)));
v(1,61)=max(max(bwlabel(m_ext2)));

v(1,62)=mean(mean(mdt(m_ext1>0)));
v(1,63)=mean(mean(mdt(m_ext2>0)));


v(1,64)=median(mdt(m_ext1>0));
v(1,65)=median(mdt(m_ext2>0));

stats=regionprops(m_ext1,'Area');
v(1,66)=mean([stats.Area]);
v(1,67)=median([stats.Area]);
if ~isempty([stats.Area])
    v(1,68)=max([stats.Area]);
end

stats2=regionprops(m_ext2,'Area');
v(1,69)=mean([stats2.Area]);
v(1,70)=median([stats2.Area]);

% PAS+ intrastructural distance features
v(1,71:71+39)=N1;


% Luminal intrastructural distance features
v(1,111:111+59)=N2;
% Nuclear intrastructural distance features
v(1,171:171+19)=N3;
% Glomerular intrastructural distance features
v(1,191:191+23)=N4;

[d1,d2,z]=size(I);

mes_int=im2double(I);
mes_int(~repmat(mes2,[1,1,3]))=NaN;
mes_int=reshape(mes_int,[d1*d2,3]);

lum_int=im2double(I);
lum_int(~repmat(white_mask,[1,1,3]))=NaN;
lum_int=reshape(lum_int,[d1*d2,3]);

nuc_int=im2double(I);
nuc_int(~repmat(nuc_mask,[1,1,3]))=NaN;
nuc_int=reshape(nuc_int,[d1*d2,3]);

%Color features
v(1,215)=mean(mes_int(:,1),'omitnan');
v(1,216)=mean(mes_int(:,2),'omitnan');
v(1,217)=mean(mes_int(:,3),'omitnan');
v(1,218)=std(mes_int(:,1),[],'omitnan');
v(1,219)=std(mes_int(:,2),[],'omitnan');
v(1,220)=std(mes_int(:,3),[],'omitnan');


v(1,221)=mean(lum_int(:,1),'omitnan');
v(1,222)=mean(lum_int(:,2),'omitnan');
v(1,223)=mean(lum_int(:,3),'omitnan');
v(1,224)=std(lum_int(:,1),[],'omitnan');
v(1,225)=std(lum_int(:,2),[],'omitnan');
v(1,226)=std(lum_int(:,3),[],'omitnan');

v(1,227)=mean(nuc_int(:,1),'omitnan');
v(1,228)=mean(nuc_int(:,2),'omitnan');
v(1,229)=mean(nuc_int(:,3),'omitnan');
v(1,230)=std(nuc_int(:,1),[],'omitnan');
v(1,231)=std(nuc_int(:,2),[],'omitnan');
v(1,232)=std(nuc_int(:,3),[],'omitnan');
v(1,233:243)=N_n;
v(1,244:254)=N_l;
v(1,255:265)=N_m;
v(1,266:285)=T_n;
v(1,286:295)=nv;
v(1,296:305)=mv;
v(1,306:315)=lv;
end