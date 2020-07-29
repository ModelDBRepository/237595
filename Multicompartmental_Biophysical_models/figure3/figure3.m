
figure;
subplot(2,2,1);boxplot(Hipp(:,3),Hipp(:,1),'Notch','on','Labels',{'sublinear','supralinear'});title('Hippocampus Length');
subplot(2,2,2);boxplot(PFCN(:,3),PFCN(:,1),'Notch','on','Labels',{'sublinear','supralinear'});title('PFC Length');
subplot(2,2,3);boxplot(Hipp(:,2),Hipp(:,1),'Notch','on','Labels',{'sublinear','supralinear'});title('Hippocampus Diameter');
subplot(2,2,4);boxplot(PFCN(:,2),PFCN(:,1),'Notch','on','Labels',{'sublinear','supralinear'});title('PFC Diameter');

figure;
subplot(1,2,1);  boxplot(ALL(:,10),ALL(:,1),'Notch','on','Labels',{'sublinear','supralinear'}); title('HPC and PFC Volume'); %volume all

subplot(1,2,2);  boxplot(ALL(:,9),ALL(:,1),'Notch','on','Labels',{'sublinear','supralinear'});title('HPC and PFC Resistance'); %Rin all

    
     % CAUSAL MANIPULATION
 load('July17_CAUSAL_MANIPULATION.mat');
figure;

subplot(1,3,1);boxplot([causalmanipSUBLINEAR(:,1),causalmanipSUPRA(:,1)],'Notch','on','Labels',{'sublinear','supralinear'});title('control');
subplot(1,3,2);boxplot([causalmanipSUBLINEAR(:,2),causalmanipSUPRA(:,2)],'Notch','on','Labels',{'sublinear','supralinear'});title('supralinear morphology');
subplot(1,3,3);boxplot([causalmanipSUBLINEAR(:,3),causalmanipSUPRA(:,3)],'Notch','on','Labels',{'sublinear','supralinear'});title('sublinear morphology');


% statistics
[h,p] = ttest2(PFC1(:,2),PFC0(:,2),'Vartype','unequal') %PFC diameter
[h,p] = ttest2(PFC1(:,3),PFC0(:,3),'Vartype','unequal') %PFC length

[h,p] = ttest2(Hipp1(:,2),Hipp0(:,2),'Vartype','unequal') %Hipp diameter
[h,p] = ttest2(Hipp1(:,3),Hipp0(:,3),'Vartype','unequal') %Hipp length


[h,p] = ttest2(ALL1(:,10),ALL0(:,10),'Vartype','unequal'); % HPC PFC volume
[h,p] = ttest2(ALL1(:,9),ALL0(:,9),'Vartype','unequal'); % HPC PFC resistance


