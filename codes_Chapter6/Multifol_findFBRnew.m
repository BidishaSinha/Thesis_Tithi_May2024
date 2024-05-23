clear;close all;
name='MDC';
reportname="region specific correlation"; 
loadten= 'n';drawnuc = 'y';lcb='n';
rep='n';
M=4;N=4;fs=19.91;order=1;


Inmin=11005.22;
conv=284.90;
counter=0;
%% FOR MULTIPLE FOLDERS with repeats
clear dfolders
answer=inputdlg('no of folders to work with?');nf=str2num(answer{1});
dfolders{1}=[];
dirname='F:\iicb\plate 1\control';
for nnf=1:nf
dfolders{nnf}=convertCharsToStrings(uigetdir(dirname));
end
[rows, nndir]=size(dfolders);
answer=inputdlg('no of repeats in each folder?');rf=str2num(answer{1});
it=nndir*rf;
%% GENERATING FILENAMES
clear fn
c=0;fn=string.empty;
for i1 = 1:nndir
    for i2=1:rf
        c=c+1;
        fn(c)=dfolders{i1};
        namef(c)=string(strcat(num2str(i1,'%.3d'), '-',num2str(i2,'%.2d')));
%         nameS(c)= strcat("PSDPCOVorder32_kap0p6mu0_",namef(c));
        nameSpC(c)=strcat(fn(c),"\set\pcov\allpsd3_cell_", num2str(i1,'%.3d'),"_", num2str(i2,'%.2d'),".mat");
        nameSpB(c)=strcat(fn(c),"\set\pcov\allpsd3_Backg_", num2str(i1,'%.3d'),"_", num2str(i2,'%.2d'),".mat");
        nameSX(c)=strcat(fn(c),"\set\pcov\allpsd2X.dat");
    end
end
%%

cfol=c;
ncells=0;   
for nn = 1:cfol
clearvars -except order d_name rf nn R fs reportname ncells loadten drawnuc rep lcb M N Inmin conv fn namef nameS nameSpC nameSpB nameSX Allfilenames dirname

 d_name=strcat(num2str(fn(nn), '%.2d'), '\set');  %% directory where files are stored
name=namef(nn);
fnSpC=nameSpC(nn);fnSpB=nameSpB(nn);fnSX=nameSX(nn);
old_dir=cd(d_name);
list=dir('*.tif');
Nims=length(list);     %% Number of images
Nims=2048;
mrf=mod(nn,rf);
if mrf>0
Nstart=(mod(nn,rf)-1)*Nims+1;
else
Nstart=(rf-1)*Nims+1;
end
Nend=Nstart+Nims-1;

A=imread(list(Nstart).name);
[row col] = size(A);
cd(old_dir)

answer = inputdlg('no of cells in this frame you want to follow for this folder?');
subcells= str2num(answer{1});
for isub=1:subcells


cd(d_name);list=dir('*.tif'); clear Im2;

%% DRAW ROIs, GET FBRs and PCOV-FBR
counter=0;
while counter == 0
%     old_dir=cd(d_name);list=dir('*.tif');
    figure 
    A=imread(list(Nstart).name);
    counter=0;
    cd(old_dir)
    imagesc(A); title (strcat('OUTER REC', num2str(isub)))
    r1 = drawrectangle('Label','Bck','Color',[1 0 0]);
    xstart=floor(r1.Position(2));
    xend=xstart+floor(r1.Position(4));
    ystart=floor(r1.Position(1));
    yend=ystart+floor(r1.Position(3));
    A=A(xstart:xend,ystart:yend);
if drawnuc == 'y'
   while drawnuc == 'y'     
    if lcb=='y' 
    %.. Load Cell boundary
    % from BW.mat saved inside pcov
    s=[d_name, '\', 'pcov','\','BW.mat'];
    load(s)
    else
    imagesc(A); title (strcat('Draw Boundary Cell:', num2str(isub)))
    h = imfreehand;
    BW = createMask(h);
    end
    BW2=imcomplement(BW);
    %..Draw Nucleus
    imagesc(A); title ("Nucleus")
    h = imfreehand;
    BWN = createMask(h);
    drawnuc = questdlg('Want to redraw?','Choose which','y','n','x');
   end
end
drawnuc='y';
old_dir=cd(d_name);Im2=[]; Im=[];A=[];Irgb=[]; IrgbM=[];
cc=0;
for i=Nstart:Nend
    cc=cc+1;
    list(i).name;A=double(imread(list(i).name));
    Im2(:,:,cc)=(A(xstart:xend, ystart:yend)-Inmin)./conv;
    Im(:,:,cc)=A(xstart:xend, ystart:yend);
end
A=imread(list(Nstart).name);
A=A(xstart:xend, ystart:yend);
cd(old_dir)
%% Input tolerance and get composite
Imax=max(Im, [], 3); maxImax=max(Imax(:));maxImaxstr=num2str(maxImax(:));
Imin=min(Im, [], 3);minImin=min(Imin(:));minIminstr=num2str(minImin(:));

Inmax= maxImax; % CHANGE HERE
Inmin= minImin; % CHANGE HERE
Tolma=maxImax/10.0;
Tolmi=minImin/10.0;
Tmax=Inmax-Tolma;
Tmin=Inmin+Tolmi; % WORK HERE


figure (1); 
subplot(2,2,1);imshow(Imax, []); 
subplot(2,2,2);imshow(Imin, []); 
subplot(2,2,3);histogram(Imax);hold on
histogram(Imin);hold off
xlabel('Intensity)'); ylabel('No. of pixels'); title(strcat('min: ', minIminstr,'  max: ',maxImaxstr));
hold on
histogram(Imin)
hold off

ImaxB=Imax;IminB=Imin; 
ImaxB(ImaxB<Tmax)=0;IminB(IminB>Tmin)=0;
ImaxB=ImaxB; IminB=IminB;%A=double(A);
%Irgb(:,:,1)=uint16(A+IminB);Irgb(:,:,2)=uint16(A+ImaxB);Irgb(:,:,3)=uint16(A);
Irgb(:,:,1)=uint16(A);Irgb(:,:,2)=uint16(A);Irgb(:,:,3)=uint16(A);


%Add the masks:
BWC=BW2+BWN;
AM=10000*bsxfun(@times, A, cast(BWC,class(A)));


IrgbM(:,:,1)=AM+uint16(IminB*10000);
IrgbM(:,:,2)=AM+uint16(ImaxB*10000);
IrgbM(:,:,3)=AM;

imshow(IrgbM);


%% FINDING FBRs (1)

clear FBR; clear SDfbr Imstd SDs2 SDt SDs
Imstd=std(single(Im2), [],3);
Imax=max(Im, [], 3);
G=Imax;
BinI= imbinarize(ImaxB*100,0.5); % max pixels as logical 1
se=strel('disk',2);
BinI2=imdilate(BinI, se);
BinI3=imfill(BinI2, 8, 'holes');
BinI4=imerode(BinI3, se);        % last three lines: filling regions inside 
                                 % max pixels with 1

AMb=imbinarize(uint16(IminB*100)+AM,0.5); % min pixels as logical 1
Fbin=BinI4+AMb;                   % adding both
FbinM=imbinarize(Fbin,0.5);
FBinM1=imclose(FbinM, 8);Inv=imcomplement(FBinM1);
Inv2 = bwareaopen(Inv,300);FBinM2=imcomplement(Inv2);
imshowpair(BinI,FBinM2,'Montage');
counter=0;
% Now place all possible FBRs and check
FBinM3=FBinM2;%imshow(IrgbM, []);
[row col]=size(A);
F4=double(zeros(row,col));
close;

%% OVERLAY FBRs
for ii=2:row-M
    for jj= 2:col-N
        FBRc=FBinM3(ii:ii+M-1,jj:jj+N-1);
        log=mean(FBRc(:));
        if log == 0
            counter=counter+1;
            FBR(counter,1)=ii;FBR(counter,2)=jj;
            FBR(counter,3)=M;FBR(counter,4)=N;
            FBinM3(ii:ii+M-1,jj)=1;
            FBinM3(ii:ii+M-1,jj+N-1)=1;
            FBinM3(ii,jj:jj+N-1)=1;
            FBinM3(ii+M-1, jj:jj+N-1)=1;
            F4(ii:ii+M-1,jj)=1;
            F4(ii:ii+M-1,jj+N)=1;
            F4(ii,jj:jj+N-1)=1;
            F4(ii+M-1, jj:jj+N-1)=1;
            
            
            SDfbr=Imstd(ii:ii+M-1,jj:jj+N-1);
            SDt(counter,2)=mean(SDfbr(:));% Sdtime
            SDt(counter,3)=median(SDfbr(:));% Sdtime
            SDt(counter,4)=std(SDfbr(:)); % Sd of Sd time (intra FBRS)
            SDt(counter,1)=counter;% FBR number
            Imfbr=Im2(ii:ii+M-1,jj:jj+N-1,1:20);
            SDs2=zeros(20,1);
            for nnn=1:20
                Imfbr=double(Im2(ii:ii+M-1,jj:jj+N-1,nnn));
                SDs2(nnn,1)=std(Imfbr(:));
            end
            SDs(counter,2)=mean(SDs2(:)); %Sdspace averaged over 20 frames
            SDs(counter,3)=median(SDs2(:));
            SDs(counter,4)=std(SDs2(:));
            SDs(counter,1)=counter;% FBR number         
        end
    end
end
% counter
IC=[];IrgbM=[];
IC(:,:,1)=A*10; %/max(A(:));
IC(:,:,2)=double(A*10)+F4.*50000.0;%/max(A(:)));
IC(:,:,3)=A*10;%/max(A(:));%/max(A(:)));
IrgbM(:,:,1)=AM+uint16(IminB*10000);
IrgbM(:,:,2)=AM+uint16(ImaxB*10000);
IrgbM(:,:,3)=AM+uint16(F4*60000);%max(AM(:));

imshow(IrgbM);
if counter ==0
    fmsgb = msgbox('No FBR found redraw again');
end
end
%% Save data
  filename= strcat('AAAna', name, 'Cell_',num2str(isub, '%.2d'));
  Allfilenames{ncells+1}=filename;
  save(filename,'maxImax','minImin', 'Tolma','Tolmi','conv','Inmin','A','Imstd','IrgbM','xstart','xend','ystart','yend','N','M','FBR','FBRc', 'FBinM3', 'SDt', 'SDs', 'dirname','d_name', 'namef', 'nameSpB','nameSpC','nameSX','fnSpB', 'fnSpC', 'fnSX', 'name', 'isub','BWC','BWN','BW2')
ncells=ncells+1;
end
end

% IF continuing after a longtime
%  Allfilenames=uigetfile({'*.mat'}, 'Select all matfiles to work woth','MultiSelect','On');
%  [rr ncells]=size(Allfilenames);loadten='n';

%%
for nn=1:ncells
 clear mpxx mpxxBack FBR
 name=Allfilenames{nn};filename= name
 %%
 load(filename);
 clear fbrRS fbrTen fbrV fbrAc fbrC rsAllf eeAllf 
 
  %fnSpC=nameSpCn(ceil(nn/2))
%  fnSpB=nameSpBn(ceil(nn/2))
%  fnSX=nameSXn(ceil(nn/2))


 Allpsd2=load(fnSpC).allpsd3; %% FIX THIS
 Allpsd=Allpsd2(xstart:xend, ystart:yend, 2:91); Allpsd2=[];
 B2=load(fnSpB).allpsd3; %% FIX THIS
 B2(:,:,1)=[];
 B=B2(:,:,1:90);B(B == 0)= NaN;Bck2=mean(mean(B,1, 'omitnan'), 'omitnan');
 B2=[];B=[];
 Bck=reshape(Bck2,90,1);Bck2=[];
 xdata2=load(fnSX);
 xdata1=squeeze(xdata2);
 xdata= xdata1(2:91);
 xx=reshape(xdata, 90, 1);xdata=[];xdata1=[];xdata2=[];
 
 [rf, cf]=size(FBR); c=0;
 
  for k=1:rf
     xf1=FBR(k,1);xf2=FBR(k,1)+FBR(k,3)-1;
     yf1=FBR(k,2);yf2=FBR(k,2)+FBR(k,4)-1;
     PSDfbr=Allpsd(xf1:xf2, yf1:yf2, :);sizePSD=size(PSDfbr);
     mPSDfbr(1:sizePSD(3),k)=reshape(mean(PSDfbr, [1 2]), 1, sizePSD(3));     
  end   
 %% PREPARE FOR FITTING
 yalls=mPSDfbr;
 yallsB=yalls-Bck-0.2;

         if loadten == 'n'
         [FfbrAc, FfbrC, FfbrRS, FfbrTen, FfbrV]=fitPSD(yallsB, xx, rf);  
         end
     
     for k=1:rf
     criR=0.999; criT= 0.05;criA=9.5; 
     if loadten == 'n'
                 if FfbrRS(k) > criR && FfbrTen(k) > criT && FfbrAc(k) < criA
                     c=c+1;
                  fbrRS(c,1)=k;fbrTen(c,1)=k;fbrAc(c,1)=k;fbrV(c,1)=k;fbrC(c,1)=k;
                  fbrRS(c,2)=FfbrRS(k);fbrTen(c,2)=FfbrTen(k);fbrAc(c,2)=FfbrAc(k);fbrV(c,2)=FfbrV(k);fbrC(c,2)=FfbrC(k);
                  fbrRS(c,3)=0;fbrTen(c,3)=0;fbrAc(c,3)=0;fbrV(c,3)=0;fbrC(c,3)=0;
                  fbrRS(c,4)=FfbrRS(k);fbrTen(c,4)=FfbrTen(k);fbrAc(c,4)=FfbrAc(k);fbrV(c,4)=FfbrV(k);fbrC(c,4)=FfbrC(k);
                 end
     end
     end
     
     c
     mAc=mean(fbrAc(:,2));mdAc=median(fbrAc(:,2));sAc=std(fbrAc(:,2));
     mTen=mean(fbrTen(:,2));mdTen=median(fbrTen(:,2));sTen=std(fbrTen(:,2));
     mV=mean(fbrV(:,2));mdV=median(fbrV(:,2));sV=std(fbrV(:,2));
     mC=mean(fbrC(:,2));mdC=median(fbrC(:,2));sC=std(fbrC(:,2));

     %% Save data
    save(strcat("B",filename),'A','Imstd','IrgbM','xstart','xend','ystart','yend', 'FBinM3', 'SDt', 'SDs', 'dirname','d_name', 'namef', 'nameSpB','nameSpC','nameSX','fnSpB', 'fnSpC', 'fnSX', 'name', 'isub','N','M','FBR','FBRc','fbrRS','fbrTen','fbrV','fbrAc','fbrC', 'FBinM3', 'SDt', 'SDs', 'xx', 'yallsB','BWC','BWN','BW2')
end