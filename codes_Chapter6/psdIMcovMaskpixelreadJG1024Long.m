function psdIMcovMaskpixelreadJG1024Long(dname, directoryname, BW, sps, conv, order, wannasave, rem)
% function psdIM2(dname, directoryname, sps, conv, order, wannasave)
% dname -> name of director in which images are stored
% sps -> sampling rate
% Maskfilename is the file of a binary image of same dimension. 1 for pixels which need to be
% fit..
% Bckfilename is the file of a binary image of same dimension. 1 for pixels
% which need to analyszed as background)
% wannasave -> 'y' or 'n' -> yes to save datafile  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=cputime
colormap(gray)
if nargin ~=8
	'usage is:'
	'psdIMcovMaskpixelread(dname, directoryname, BW, sps, conv, order, wannasave, rem)'
	return
end

if ~exist(dname,'dir')
	sprintf('%s: directory does not exist',dname)
	return;
end 
%%
old_dir=cd(dname)
s5 = sprintf('%03d',order);
%% ..............................READING First FILE........................
listings=dir([dname,'\*.tif']);
nfiles=length(listings);L=nfiles;
%%
A1=imread([dname,'\',listings(1).name]);


[row, col]=size(A1);
Nr=row-2; Mr=col-2;
lens=length(dname);
strd=dname(lens-6:lens-4);
    reps=nfiles/2048;
    %%
for itra = 1:reps
%% .............................Check Size......................
dim=Nr*Mr;
lim=1047553;
if dim<lim
    allpsd3=zeros(Nr,Mr,1025); 
end


%% ...................READING ALL FILEs. storing in imgdata(:, :, nfiles) part by part.......................
clear xb
clear yb

xb = floor((Nr-1)/512)+1;
yb = floor((Mr-1)/512)+1;
NMall="x";
NMall=[];


for it1=1:xb
    for it2=1:yb
        if it1==xb 
            N2=Nr;
        else N2=512;
        end    
        if it2==yb 
            M2=Mr;
        else M2=512;
        end
        N1=(it1-1)*512+1;
        M1=(it2-1)*512+1;
        N2;
        M2;
        NMall=[NMall; N1, N2, M1, M2];
        clear imgdata
        
        mrf=mod(itra,reps);
        if mrf>0
        Nstart=(mod(itra,reps)-1)*2048+1;
        else
        Nstart=(reps-1)*2048+1;
        end
        Nend=Nstart+2048-1;
        

for ii=Nstart:Nend
	A=imread([dname,'\',listings(ii).name])./conv; %320 is for covertin to h
	B=A;
	imgdata(:,:,ii)=B(N1:N2,M1:M2) ;
    BW2=BW(N1:N2,M1:M2);
end
[row, col]=size(B(N1:N2,M1:M2));
N=row; M=col;

%...................PSD...x: ff, y=psds....................

%%
% L = nfiles;                     % Length of signal
% NFFT = 2^nextpow2(L); % Next power of 2 from length of y
% 
% ff = sps/2*linspace(0,1,NFFT/2+1);
% psds=zeros(NFFT/2+1,M*N);
psds2=zeros(1025,M*N);
% allpsd=zeros(N,M,NFFT/2+1);
allpsd2=zeros(N,M,1025);

c=0;
for ii=1:M
	for jj=1:N
        c=c+1;
        check1=BW2(jj, ii);
        if check1 
            
            x=squeeze(imgdata(jj,ii,:));
            xx=detrend(x, 'constant'); %vari(c)=var(xx); 
           [F1, f]=pcov(double(xx), order, 2048, sps);
            psds2(:,c)= F1;
            allpsd2(jj,ii,:)=F1;    
        end
	end
end
if dim<lim
    allpsd3(N1:N2, M1:M2, :)=allpsd2;
else
filename=['allpsd2_' sprintf('%d',N1) '_' sprintf('%d',N2) '_' sprintf('%d',M1) '_' sprintf('%d',M2)  rem '.mat'];
save(filename,'allpsd2', '-v7.3');
end

clear allpsd2
    end
end
%%
cd(old_dir)
if wannasave == 'y'
    if dim<lim
    filename=[directoryname '\pcov\allpsd3_' rem '_' strd '_' sprintf('%.2d', itra) '.mat']
    % filename=['allpsd3_' sprintf('%d',Nr) '_' sprintf('%d',Mr)  rem '.mat'];
    save(filename,'allpsd3', '-v7.3');
    end
filename2=[directoryname '\pcov\allpsd2X.dat'];
% filename2=['allpsd2X.dat'];
dlmwrite(filename2,f);

end

% % fts(:,1)=ff;fts(:,2)= mean(psds,2);
% fts2(:,1)=f;fts2(:,2)= mean(psds2,2);
% 
% %...................plot PSD and check variance match.......................
% figure(4)
% % meanvari=mean(vari)
% % loglog(fts(:,1),fts(:,2), 'o'); hold on;
% loglog(fts2(:,1),fts2(:,2)); hold on;
% xlabel('Frequency');
% ylabel('PSD (nm2-sec)')

end
% %........................................................................
NMall
e=cputime-t
%% everying must end even functions!
cd(old_dir)
pwd

