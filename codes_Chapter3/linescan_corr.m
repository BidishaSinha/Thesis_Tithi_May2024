%To choose a point, move your cursor to the desired location and press 
% either a mouse button or a key on the keyboard. Press the Return key 
% to stop before all n points are selected.

% eeAll=table2array(eeAll1);
clear;clc
clearvars -except eeAll
olddir=cd('I:\tension map_ 02112021\plate 2')
name='fixed_before';
I=imread('003_10000.tif');c=0;

%loading tension matrix of sme time point as eeAll and that of control
%timepoint as eeAllold
% m = matfile('eeAll.mat');
l=5; % half-length of linescan : length will be 2l+1
w=5; % half width of linescan : width will be 2w+1
TF=exist('eeAll_001');
if TF == 0
    m = matfile('eeAll_001.mat');
    eeAll = m.eeAll;
end
% % m2=  matfile('Singlepixelcheck_Fitparas_kap0p6mu0_old.mat');

% % eeAllold = m2.eeAll;
% % % m = matfile('eeAll(control).mat');
% % % eeAllold = m.eeAll;

figure(1)
%                                 imshow(eeAll*50, [1 200]); colormap jet
                               imshow(I, []); colormap jet;
%  [py,px] = ginput

load('fixed_after.mat', 'px', 'py');
np=max(size(px))
% eeAll2=eeAll.*50;px=int16(px);py=int16(py);
eeAll2=eeAll;px=int16(px);py=int16(py);
% eeAll3=eeAllold.*50;
rp=zeros(np,2);xsc=zeros(np, 2*l+1);ysc=zeros(np, 2*l+1);fscx=zeros(np, 2*l+1);fscy=zeros(np, 2*l+1);
Tscx=zeros(np, 2*l+1);Tscy=zeros(np, 2*l+1);
for it=1:np
    z1=eeAll2(px(it)-l:px(it)+l, py(it)-w:py(it)+w);
    zf1=I(px(it)-l:px(it)+l, py(it)-w:py(it)+w);
    Tscx(it, :)=mean(z1');fscx(it, :)=mean(zf1');
    z2=eeAll2(px(it)-w:px(it)+w, py(it)-l:py(it)+l);
    zf2=I(px(it)-w:px(it)+w, py(it)-l:py(it)+l);
    Tscy(it, :)=mean(z2); fscy(it, :)=mean(zf2);
    
    yyaxis left;plot([-1*l:l], fscx(it,:), '--o');ylim([0 500]);
    yyaxis right;plot([-1*l:l], Tscx(it,:),'--o');ylim([0 500]); 
    
    pause(0.2);
    yyaxis left;plot([-1*l:l], fscy(it,:), '--o');ylim([0 500]);ylabel('Fluorescence');
    yyaxis right;  plot([-1*l:l], Tscy(it,:),'--o');ylim([0 500]); ylabel('Tension');
    pause(0.2);
%     corr for every linescan
    f1=fscx(it,:);f2=fscy(it,:);T1=Tscx(it,:);T2=Tscy(it,:);
%     [rho,pval] = corrcoef(double([f1(T1>1) f2(T2>1)]), double([T1(T1>1) T2(T2>1)]));
z1=z1(:); zf1=zf1(:); z2=z2(:); zf2=zf2(:); 
    [rho,pval] = corrcoef(double([zf1(z1>1)' zf2(z2>1)']), double([z1(z1>1)', z2(z2>1)']));
    rp(it, 1)= rho(1, 2);rp(it, 2)= pval(1, 2); 
end

szx=max(size(fscx(Tscx>1)))*min(size(fscx(Tscx>1))); fx1=fscx(Tscx>1);Tx1=Tscx(Tscx>1);
szy=max(size(fscy(Tscy>1)))*min(size(fscy(Tscy>1))); fy1=fscy(Tscy>1);Ty1=Tscy(Tscy>1);
% figure(2);loglog([squeeze(fscx(Tscx>1)); squeeze(fscy(Tscy>1))]', [squeeze(Tscx(Tscx>1)); squeeze(Tscy(Tscy>1))]', 'o')
figure(2);loglog([reshape(fx1, [szx,1]); reshape(fy1, [szy,1])]', [reshape(Tx1, [szx,1]); reshape(Ty1, [szy,1])]', 'o')

[rho,pval] = corrcoef(double([fscx fscy]), double([Tscx Tscy]))
[rho2,pval2] = corrcoef(double([fscx(Tscx>1); fscy(Tscy>1)]), double([Tscx(Tscx>1); Tscy(Tscy>1)]))
% [rho3,pval3] = corrcoef(double([fscx(oTscx>1); fscy(oTscy>1)]), double([oTscx(oTscx>1); oTscy(oTscy>1)]))
mean(rp)
save(name)
cd(olddir)
name2=[name 'rp']
save(name2, 'rp')
