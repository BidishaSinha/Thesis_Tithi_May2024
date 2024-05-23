clear;clc;close all
wannasave = 'f' ;% write 't' to save
rr=1; % No. of ROIs to follow


readrois='t'; % write t or f
for a =1:1; % was changed to 1:7 to generate pvalue image
    ThC=50*a;
    for b=7:7; % was changed to 1:7 to generate pvalue image
        ThP=50*b;

c=0;
for j=1:1 % change it when doing all cells
    for jjj=1:4
clearvars -except c ThC ThP a b readrois j jj jjj xstart1 xend1 ystart1 yend1 Npall ni rr wannasave Ndall Adall Ndall_rab4_5 Adall_rab4_5 in1 in2 pv rh countsim calsim
if jjj==1 jj=0;
elseif jjj==2 jj=5;
    elseif jjj==3 jj=10;
        elseif jjj==4 jj=30;
end
c=c+1;
dirname1=['cell_' num2str(j) ' ca_' num2str(jj) ' min.tif']; %uigetfile('*.tif');   %'H:\09122023\analysis\cell_1 ca_0 min.tif';
dirname2=['cell_' num2str(j) ' pi_' num2str(jj) ' min.tif'];%uigetfile('*.tif');   %'H:\09122023\analysis\cell_1 pi_0 min.tif';
dir = [dirname1;dirname2];
%olddir=cd(dirname);%list=dir('*.tif');
for n=1:2; %frames
I1(:,:,n)=imread(dir(n,:));
filename{n}=dir(n,:);
end
% Imean=mean(I1(:,:,1));
 Imean=0.5.*I1(:,:,2)+0.5.*I1(:,:,1);
if readrois=='f'
figure (1); imshow(Imean, []); colormap jet; title ("Average of all frames")
r1 = drawrectangle('Label','Bck','Color',[1 0 0]);
    xstart=floor(r1.Position(2));xstart1(c)=xstart;
    xend=xstart+floor(r1.Position(4));xend1(c)=xend;
    ystart=floor(r1.Position(1));ystart1(c)=ystart;
    yend=ystart+floor(r1.Position(3));yend1(c)=yend;
else
    load('rois4timepoints.mat');xstart=xstart1(c);xend=xend1(c);ystart=ystart1(c);yend=yend1(c);

end

close
frames=2;
for n=1:2;%n=1 will read calcium, n=2 will read Piezo
   
    I2=I1(:,:,n);
    II=I2(xstart:xend,ystart:yend);
    I=II;

    fig2=figure(10);
    
    subplot(2,4,4*n-3);imshow(I, [100 0.75*max(I(:))]);colormap jet; colorbar
    subplot(2,4,4*n-2);imshow(I, [100 0.75*max(I(:))]);colormap jet; colorbar
    if n==1
    T = adaptthresh(I, 0.6,'Statistic', 'median');
    BW2 = imbinarize(I,T);
    elseif n==2
        T = adaptthresh(I, 0.6,'Statistic', 'median');
    BW2=imbinarize(I,T);
    end
    % Convert image to binary image, specifying the threshold value.
 
    se = strel('disk',1);BW2b=imfill(BW2,'holes');se2 = strel('disk',1);
BW3=imerode(BW2b, se);BW4=imdilate(BW3, se2);BW5=imfill(BW4, 'holes');
cc = bwconncomp(BW5); 
stats = regionprops(cc,'Area','Eccentricity','PixelIdxList'); 
if n==1
 idx = find(([stats.Area] > 10 )&([stats.Area] < ThC ));nm="calcium";
elseif n==2   
 idx = find(([stats.Area] > 10 )&([stats.Area] < ThP )); nm="Piezo";
end
BW = ismember(labelmatrix(cc),idx); LL=bwlabeln(BW);
cc2 = bwconncomp(BW); 
BWall(:,:,n)=BW;% BWall saves all the particles detected
stats = regionprops(LL,'Area','Eccentricity','PixelIdxList', 'BoundingBox');
stats_AAA = regionprops('table',LL,'Area','Eccentricity','PixelIdxList', 'BoundingBox');
statsall{n}=stats; % saves all the stats
%
np=max(LL(:)); clear int
Np(n)=np; area=(xend-xstart)*(yend-ystart);
Npall(c,n)=np;
Ndall(c,n)=np/(area*0.065*0.065);
Adall(c,n)=sum(stats_AAA.Area)/area;

subplot(2,4,4*n-1);imshow(LL, []);colormap jet;colorbar;
%   pause(2)
end

BW1111=BWall(:,:,1)+BWall(:,:,2);
[m1,n1] = size(BW1111);
%simulating calcium spots now randomly placed
sizsim=floor(sqrt(floor(Adall(c,1)/(Ndall(c,1)*0.065*0.065))));
ss=floor(sizsim/2);
BWsim=zeros(m1, n1);
for k = ss+1 : m1-ss
    for l = 1+ss : n1-ss
        if rand()<(Adall(c,1)/(2*ss+1)/(2*ss+1))
        BWsim(k-ss:k+ss,l-ss:l+ss) = 1;
        end
    end
end
BW1sim=BWall(:,:,2)+BWsim;
 countsim(c)=sum(BW1sim(BW1sim==2))/2/area; % Overlap area fraction
calsim(c)=sum(BWsim(BWsim==1))/area; % simulated cakcium area fraction

count1 = 0;
for k = 1 : m1
    for l = 1 : n1
        if BW1111(k,l) == 2
            count1 = count1 + 1;
            BW2222(k,l) = 1;
        else
            BW2222(k,l) = 0;
        end
    end
end
 figure(10)

set(gcf, "Position", [100 100 1000 800]);
Irgb(:,:,1)=uint16(BWall(:,:,1));Irgb(:,:,2)=uint16(BWall(:,:,2));Irgb(:,:,3)=uint16(BWsim);colorbar;
subplot(2,4,8);imshow(60000*Irgb, []);title(['cell_' num2str(j) ' capi' num2str(jj)]); colorbar;
 figure(11)
 fig=imshow(60000*Irgb, [])
 title(['cell_ ' num2str(j) ' capi' num2str(jj) ' min.tif']);
saveas(fig2, ['ALL50cell_' num2str(j) ' capi' num2str(jj) ' min.tif'], 'tiffn')
%   saveas(fig, ['RESALL50cell_' num2str(j) ' capi' num2str(jj) ' min.tif'], 'tiffn')
imwrite(60000*Irgb, ['ImRESALL50cell_' num2str(j) ' capi' num2str(jj) ' min.tif'])

stats_AAA111 = regionprops('table',BW2222,'Area','Eccentricity','PixelIdxList', 'BoundingBox');
[x111,t111] = size(stats_AAA111);
Ndall_rab4_5(c)=x111/(area*0.065*0.065); % overlapping
Adall_rab4_5(c)=sum(stats_AAA111.Area)/area;%overlapping
% % %cd(olddir);name2=list(n).name;change name accordingly the number end-10;

% II1=I1(:,:,1);II2=I1(:,:,2);
% for ij=1:x111
%     id1=stats_AAA111.PixelIdxList;
%     in1(ij,c)=sum(II1(id1{ij}));in2(ij,c)=sum(II2(id1{ij}));
% end
% if x111 > 10
% [p, rho]=corrcoef(in1, in2);
% pv(c)=p(2,1);rh(c)=rho(2,1);
% end

name=[ dirname1(end-10:end-4) '_ALL_' num2str(xstart, '%04d') num2str(ystart,'%04d' ) num2str(xend,'%04d') num2str(yend,'%04d')];
if wannasave == 't'
save(name)
end

imshow(BW1111, [])
    end
end
%%
if readrois=='f'
save('rois4timepoints.mat','xstart1', 'xend1', 'ystart1', 'yend1')
end
figure(2)
% subplot(2,2,1);plot(Adall(:,2),100*(Adall_rab4_5)'./Adall(:,1), 'o')
subplot(2,2,1);scatter(Adall(:,2), Adall(:,1), 'o', 'filled')
hold on 
subplot(2,2,1);scatter( Adall(:,2),calsim, 'o', 'filled')
xlabel("Area fraction (Piezo peaks)"); ylabel( "Area fraction (Calcium peaks)");
legend(["Data" "Random regions simulated"], 'Position',[0.17 0.93 0.27 0.063])
% set(legend1,...
%     'Position',[0.17 0.93 0.27 0.063]);

subplot(2,2,2);scatter(Adall(:,2),100*(Adall_rab4_5)'./Adall(:,1), 'o', 'filled')
xlabel("Area fraction (Piezo peaks)");ylabel( "% Area of Calcium peaks in overlap");
hold on
subplot(2,2,2);scatter(Adall(:,2),100*(countsim)'./Adall(:,1), 'o', 'filled')
legend(["Data" "Random regions simulated"], 'Position',[0.6 0.93 0.27 0.063])

subplot(2,2,3);imshow(BWsim, []); title("Example Simulated data")

subplot(2,2,4);imshow(BW1111, []);title("Corresponding real data")
% bar(ff,mean(Ndall))
%[p15, h15]=ranksum(Ndall(:,1), Ndall(:,5))
%[p25, h25]=ranksum(Ndall(:,2), Ndall(:,5))
%%
% 
% for i=1:length(Ndall)
%     plot(ff,Ndall(i,:));
%     hold on;
% end
% hold off

pv(a,b,1)=ranksum(countsim, (Adall_rab4_5)')
pv(a,b,2)=ThC;
pv(a,b,3)=ThP;
mean(mean(BWall(:, :,1)))
mean(BWsim(:))

    end
end

% if you want to go through many threholds uncomment last line
% figure(12); imagesc(pv(:,:,1)); colormap jet;color bar