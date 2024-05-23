clear
load('chemokine.mat')
allchemokineten=Par{:,3};
allchemokineSDt=Par{:,8};
load('control.mat')
allcontrolten=Par{:,3};
allcontrolSDt=Par{:,8};
comp1=allcontrolSDt(:,2);
comp2=allchemokineSDt(:,2);
ranksum(comp1,comp2)
median(comp1)
median(comp2)
mean(comp1)
mean(comp2)
std(comp1)
std(comp2)

histogram(comp1,'Normalization', 'probability', 'FaceAlpha',0.5, 'FaceColor', 'r', 'BinWidth',0.5 )
hold on
histogram(comp2,'Normalization', 'probability', 'FaceAlpha',0.5, 'FaceColor', 'g', 'BinWidth',0.5 )
hold off
ax = gca;
ax.FontSize = 15; 
xlabel('SD in nm','FontSize', 20);
ylabel('Probability', 'FontSize', 20)
%%
A=comp1; B=comp2;
lenA=length(comp1); lenB=length(comp2);
lenmax=max(lenA, lenB);
A1 = padarray(A,(lenmax-lenA),0/0, 'post');B1 = padarray(B,lenmax-lenB,0/0, 'post');
figure(2)



