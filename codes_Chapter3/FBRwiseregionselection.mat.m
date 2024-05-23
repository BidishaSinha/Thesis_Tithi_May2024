clear;clc;
% dirname1=  'E:\MBCD\2 hr\control\region wise\001\set\';%Contrl
%  dirname2= 'E:\MBCD\2 hr\control\region wise\001\set\fl\';
filename= 'K:\tf correlation with tirf cholesterol\26052022\mbcd\plate 1\001\BAAAna001-01Cell_01.mat';
%753.95556 : Imin
% 2546.84666 : Imax
% 24.64406 :slope
 inmin=0;
conv=1;nfiles=200;
% olddir=cd('G:\tithi\atp_depleted\06122020\control\analysis_with_one_cell\001\');
A1=imread('K:\tf correlation with tirf cholesterol\26052022\mbcd\plate 1\001\0010001.tif');
A2=imread('K:\tf correlation with tirf cholesterol\26052022\mbcd\plate 1\001\002_translate.tif');
load(filename);

%%
ijxt=0;ijyt=0;
A11=A1(xstart+ijyt:xend+ijyt, ystart+ijxt:yend+ijxt);
A22=A2(xstart:xend , ystart:yend);
figure (1);imshowpair(A11,A22, 'ColorChannels', 'green-magenta')
figure (2);imshow(A11-A22, []); colormap jet


%%



% olddir=cd(dirname1);
list1=dir('*.tif');
% A1=imread(list1(1).name);
% cd(dirname2);
%list2=dir('*.tif');
%A2=imread(list2(1).name);


% cd(olddir)

%% Checking if late FBRs are still falling inside. chnage 90 to whatever
M=FBR(1,3);N=FBR(1,4);
%index11 = fbrTen(1,1)
x=FBR(1,1); y=FBR(1,2);
aa=A11(x:x+M, y:y+M);
figure, imshow(A11, []);colormap jet
    h = imrect(gca, [y x 24 24]);
    %%
    tic
%     cd(dirname1)
    Im1=double(zeros(xend-xstart+1, yend-ystart+1,nfiles));
    A1a1=imread(list1(1).name);[xsize, ysize]=size(A1a1);
    A1b=double(zeros(xend-xstart+1, yend-ystart+1));
    A1a=double(zeros(xsize, ysize));
    nfiles=1;% make it 2048
    for n= 1:nfiles
        A1a=imread(list1(n).name);
        A1b=(A1a(xstart+ijyt:xend+ijyt, ystart+ijxt:yend+ijxt)-inmin)./conv;
        Im1(:,:,n)=A1b;
    end
    toc
    %%
%         cd(dirname2)
%     Im2=double(zeros(xend-xstart+1, yend-ystart+1,nfiles));
%     
%     for n= 1:nfiles
%         A1a=imread(list2(n).name);
%         A1c=(A1a(xstart:xend, ystart:yend)-inmin)./conv;
%         Im2(:,:,n)=A1c;
%     end
    %% Extracting patches
    
%     cd(olddir)
    [nFBR, zz]=size(fbrTen);
%     data = []
counter11 = 0;
    for i=1:nFBR
    index = fbrTen(i,1);
    x=FBR(index,1); y=FBR(index,2);
    aa=A1b(x:x+M-1, y:y+N-1);
%     figure, imshow(A1b, []);colormap jet;
%     h1 = imrect(gca, [y x M N]);
    x111(i,1) = x;
    y111(i,1) = y;
%     figure(1);
    
            f=[ 'Roi' num2str(i, '%03.f') 'x' num2str(x, '%03.f') num2str(y, '%03.f') '_' num2str(M, '%02.f') '.mat'];
%             filename2=['Patches3' '/' f];
            counter11 = counter11+1;
%             filename=convertCharsToStrings(f);
%             h5create(filename2,'/dataset2',[M N 1])% make it [M N 2048]
            mydata = Im1(x:x+M-1, y:y+M-1,:);
            mydata1 = mydata(:,:,1);
            Res=mean(mydata1(:));
            data(counter11,:)=Res;
%             h5write(filename2, '/dataset2', mydata)
    
    end
%%
imshow(A22, [0 1200]);
hold on
for iii = 1:length(x111)
   xn = x111(iii,1);
   yn = y111(iii,1);
   rectangle('Position',[yn xn M N],'EdgeColor','y','LineWidth',2)
%    h1 = imrect(gca, [yn xn M N]);
end
saveas (gcf,'If1.jpg')
imshow(A11,[]);
hold on
for iii = 1:length(x111)
   xn = x111(iii,1);
   yn = y111(iii,1);
   rectangle('Position',[yn xn M N],'EdgeColor','y','LineWidth',2)
%    h1 = imrect(gca, [yn xn M N]);
end
saveas (gcf,'IRM2.jpg')
imshow(A22, [0 1200]);colormap jet
hold on
for iii = 1:length(x111)
   xn = x111(iii,1);
   yn = y111(iii,1);
   rectangle('Position',[yn xn M N],'EdgeColor','y','LineWidth',2)
%    h1 = imrect(gca, [yn xn M N]);
end
saveas (gcf,'lut.jpg')

%%
data1=data';
save('data_new_before.mat','data');
