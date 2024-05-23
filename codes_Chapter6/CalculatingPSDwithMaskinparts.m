%% inputs

d_name='H:\iicb\plate 1\treat\005\set'; % PLEASE KEEP SINGLE INVETED COMMAS NOT DOUBLE
old_dir=cd(d_name);list=dir('*.tif');L=length(list); %cd(old_dir)

listings=dir([d_name,'\*.tif']);
nfiles=length(listings);
A1=imread([d_name,'\',listings(1).name]);
A2=imread([d_name,'\',listings(nfiles).name]);
Afl=A1+A2;

[N, M]=size(A1);
%% Make a mask 
cd(old_dir)
nobj= 2; % No. of objects to mask out separately
clear BW
clear BWs
for it=1:nobj
    figure(1)
    imagesc(Afl); title ("Draw Cell Boundary")
    h = imfreehand;
    BWi = createMask(h); BWs(:,:,it)=BWi;
    imshow(BWi,[]);
end
BW=max(BWs, [], 3);
imshow(BW)
clear BWs
clear BWi

%% Make a background mask
nobjb= 1; % No. of objects to mask out separately
clear BWbck
clear BWs
for it=1:nobjb
 imagesc(Afl); title ("Draw Background ROI")
    h = imrect;
    BWbcki = createMask(h); BWbcks(:,:,it)=BWbcki;
    imshow(BWbcki,[]);
end
    BWbck=max(BWbcks, [], 3);
    imshow(BWbck,[]);
clear BWbcki
clear BWbcks
cd(old_dir)

%% Save the files
d_name2=strcat(d_name,'\pcov'); mkdir(d_name2)
save(strcat(d_name,'\pcov\BW.mat'), 'BW', 'BWbck'); % change backslash in case of PC
%%
%psdIMcovMaskpixelread(dname, directoryname, BW, sps, conv, order, wannasave, rem)
tic
% cd(old_dir)
sps = 19.91; cf = 284.90;
wannasave = 'y'; order = 32;
dname=d_name;
directoryname=d_name;
rem1= 'cell';
psdIMcovMaskpixelreadJG1024Long(d_name,d_name, BW, sps, cf, order, wannasave, rem1) 
toc
%% 
%psdIMcovMaskpixelread(dname, directoryname, BW, sps, conv, order, wannasave, rem)
cd(old_dir)
sps = 19.91; cf = 284.90;
wannasave = 'y'; order = 32;
dname=d_name;
directoryname=d_name;
psdIMcovMaskpixelreadJG1024Long(d_name,d_name, BWbck, sps, cf, order, wannasave, 'Backg')