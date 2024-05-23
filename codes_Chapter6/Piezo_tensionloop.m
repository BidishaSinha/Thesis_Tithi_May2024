clear
olddir=cd('G:\Shared drives\Tithi Work\IICB_Chinky\polarity analysis\plate 2');
X = [0 3 6 9 15 26 30];
X1 = [0 4 7 10 16 27 31];
% olddir=cd('G:\Shared drives\Tithi Work\IICB_Chinky\polarity analysis\plate 4');
% X = [0 1 5 15 30];
% X1 = [0 2 6 16 31];
% olddir=cd('G:\Shared drives\Tithi Work\IICB_Chinky\polarity analysis\plate 1');
% X = [0 1 ];
% X1 = [0 2 ];
% olddir=cd('G:\Shared drives\Tithi Work\IICB_Chinky\polarity analysis\plate 3');
% X = [0 1 ];
% X1 = [0 2 ];
F=dir('*.tif');
M=dir('*.mat');
nfiles=length(F);
for n=1:nfiles

load(M(n).name)
 Tten=Ten;
% Tten=Ten(191:370,155:315);
[row, col]=size(Tten);
I1=double(imread(F(n).name));
I=I1(1:row,1:col);
% I=I1(191:370,155:315);
% imshow(I, []);colormap jet; colorbar
figure (1)
imshow(Tten, [0 250]);colormap jet; colorbar
threshI=1550;
Tten2=Tten(Tten>1 & I>threshI);
I2=I(Tten>1 & I>threshI);
pn(n)=length(I2(I2<1650 & I2>1500));%(I2<8650));
pn2(n)=mean(I2);%(I2(I2<1650 & I2>1500));
% loglog(I2, Tten2, 'o')
mten(n)=mean(Tten2);%(I2<8650));
nn=600;
Ten=Tten2;
eaa2=I2;
t=zeros(nn,1);Tn=threshI;
step=2.55;%(max(Ten)-Tn)/n;
j=eaa2;
for i= 1:nn
    t(i)=mean(Ten(j>step^(i-1)+Tn & j<(step^i)+Tn));
    msd(i)= std(Ten(j>step^(i-1)+Tn & j<(step^i)+Tn));
    ea(i)=mean(eaa2(j>step^(i-1)+Tn & j<(step^i)+Tn));
    emsd(i)= std(eaa2(j>step^(i-1)+Tn & j<(step^i)+Tn));
%     t(i)=mean(Ten(Ten>step*(i-1)+Tn & Ten<(step*i)+Tn));
%     ea(i)=mean(eaa2(Ten>step*(i-1)+Tn & Ten<(step*i)+Tn));
end
c=0;
for i=1:nn
    if ~isnan(ea(i)) && ~isnan(t(i))
        c=c+1;
        eaNew(c)=ea(i); 
        tNew(c)=t(i); 
        msdNew(c)=msd(i);
        emsdNew(c)=emsd(i);
    end
end

figure (2)
plot(eaNew,tNew, '-s','MarkerSize',5,...
    'LineWidth',5) %'MarkerEdgeColor','red','MarkerFaceColor','red',
% set(gca, 'XScale', 'log')
xlabel ('Intensity (a.u.)', 'fontsize', 20)
ylabel('Tension (pN/\mum)','fontsize', 20)
hold on

figure(5); [f,xi] = ksdensity(I(Tten>1)); 
semilogx(xi,f);hold on; ylabel("Probability"); xlabel("Tension (pN/\mum)")

end
hold off

figure (4) 
p=pn;
plot(X1,(p-min(p))/(max(p)-min(p)),'-s','MarkerSize',5, 'linewidth',1.5)
hold on
plot(X,(mten-min(mten))/(max(mten)-min(mten)),'-s','MarkerSize',5, 'linewidth',1.5)
hold off
xlabel ("Time (min)")
ylabel ("Normalized Piezo (blue); Tension (red)")

% 
% % delten=mten(1:end-1)-mten(2:end);
% % piezo=pn2(1:end-1);
% % [fitresult, gof] = createFit(piezo, delten)
% delpiezo=eaNew(1:end-1)-eaNew(2:end);
% delten=(tNew(1:end-1)-tNew(2:end))./delpiezo;
% delpiezo=eaNew(1:end-1)-eaNew(2:end);
% piezo=eaNew(1:end-1);
% [fitresult, gof] = createFit(piezo, delten)
% 
cd(olddir)