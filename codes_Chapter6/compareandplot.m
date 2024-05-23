
load('chemokine.mat')
allchemokineten=Par{:,3};
allchemokineSDt=Par{:,8};
load('control.mat')
allcontrolten=Par{:,3};
allcontrolSDt=Par{:,8};
ranksum(allcontrolten(:,2),allchemokineten(:,2))
median(allcontrolten(:,2))
median(allchemokineten(:,2))
mean(allcontrolten(:,2))
mean(allchemokineten(:,2))
std(allcontrolten(:,2))
std(allchemokineten(:,2))

histogram(log10(allcontrolten(:,2)),'Normalization', 'probability', 'FaceAlpha',0.5, 'FaceColor', 'r', 'BinWidth',0.2 )
hold on
histogram(log10(allchemokineten(:,2)),'Normalization', 'probability', 'FaceAlpha',0.5, 'FaceColor', 'g', 'BinWidth',0.2 )
hold off
ax = gca;
ax.FontSize = 15; 
xlabel('log_1_0(Tension values in 10^-^6 N)','FontSize', 20);
ylabel('Probability', 'FontSize', 20)
%%
A=allcontrolten(:,2); B=allchemokineten(:,2);
lenA=length(allcontrolten(:,2)); lenB=length(allchemokineten(:,2));
lenmax=max(lenA, lenB);
A1 = padarray(A,(lenmax-lenA),0/0, 'post');B1 = padarray(B,lenmax-lenB,0/0, 'post');
figure(2)



