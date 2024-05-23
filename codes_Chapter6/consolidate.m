clear
d='D:\analysis_jurkat\chemokine\Field ';
counter=0;
counterfol=0;
for f=1:4
    for c=1:2
dirname= [d num2str(f) ' Cell ' num2str(c)]
if isfolder(dirname)
    counterfol=counterfol+1;
old_dir=cd(dirname);
cellname= ['BAAAna00' num2str(f) '-01Cell_0' num2str(c)];
if isfile([cellname '.mat'])
counter=counter+1;
load(cellname);
Par{counter,3}=fbrTen;
Par{counter,4}=fbrAc;
Par{counter,5}=fbrV;
Par{counter,6}=fbrC;
Par{counter,7}=fbrRS;
Par{counter,8}=SDt;
Par{counter,9}=SDs;
Par{counter,1}=dirname;
Par{counter,2}=cellname;
end
end
    end
end

cd('D:\analysis_jurkat')
save('chemokine','Par')