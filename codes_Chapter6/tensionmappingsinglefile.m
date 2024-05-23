%% When you have single files tensionmapping_loop(v)
backname='allpsd3_back_1.mat';
B1=load(backname);%background file%
B2=B1.allpsd3,'-v7.3';B2(:,:,1)=[];
B=B2(:,:,1:90);B(B == 0)= NaN;Bck2=mean(mean(B,1, 'omitnan'), 'omitnan');
Bck=reshape(Bck2,90,1);
cellname='allpsd3_cell_1.mat';
Aall=load(cellname);%cell file%
alls=Aall.allpsd3;
xxxdata='allpsd2X.dat';
xdata1=squeeze(load(xxxdata));
xdata= xdata1(2:91);
xx=reshape(xdata, 90, 1);
alls(:,:,1)=[];
yalls=alls(:,:,1:90);
[row,col,len]=size(alls);
aAll = zeros(row,col); bAll = zeros(row,col); cAll = zeros(row,col); dAll = zeros(row,col); 
eAll = zeros(row,col); eeAll = zeros(row,col); rsAll = zeros(row,col);test = zeros(row,col)
run1=uint16(col/14);
for ii = 1:row
    kk=ii;
    parfor jj = 1:col
        if alls(kk,jj,2) ~= 0
            test(kk,jj)=1;
 
        A=reshape(yalls(kk, jj, :),90,1);
        ydata=A-Bck-0.2;%backgroundsubstraction -B(2:size(A),2)-B(2:size(A),2)
        yy=reshape(ydata,90,1);
      
     %....fitting.....
 
        w=1./(xx);
        [xData, yData, weights] = prepareCurveData( xx, yy, w );
        g=fittype('psdsimfunc6(AT,cytvis,ff,sig,x)');
        topts = fitoptions( g );
       % opts.Display = 'Off';
        topts.Lower = [1,1,0,0];%changed here
        topts.StartPoint = [1,20000,0,1];
        topts.Upper = [10,100000000,1,1000];
        topts.Weights = weights;
%         opt.Method='NonlinearLeastSquares';
        
        [f,g]=fit(xx, yy, g, topts);            
        a=f.AT;b=f.cytvis;c=f.ff;ee=f.sig;
        
        aAll(kk, jj)=a;%%temp
        bAll(kk,jj)=b;%%cytvis
        cAll(kk, jj)=c;%%confinement
        eeAll(kk, jj)=ee;%%tension
        rsAll(kk, jj)=g.rsquare;%% rsquare
        s3=['ee =', num2str(ee)];
        end
        end
        
       
    end
% end
clear alls; clear alls2;clear A; clear Bck;
dAll(aAll>0)=0.6;
filename=sprintf('Singlepixelcheck_Fitparas_kap0p6mu0_1.mat');
save(filename);   
%Change the clims as required. Change them for mapping each parameter in the input value of parameter
%Map_PSD(aAll,'Temperature',v,1,10)
%Map_PSD(bAll/100,'cytoplasmicViscosity',v,1,10000)
%Map_PSD(cAll,'Confinement',v,0,1)
%Map_PSD(eeAll*50,'tension',v,0,1000)
%Map_PSD(rsAll,'Rsqr',v,0,1)
            