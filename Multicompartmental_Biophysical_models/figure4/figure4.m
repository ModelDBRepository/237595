%Generate Figure 4

load('figure4.mat');

A=figure;

subplot(2,2,1); 


shadedErrorBar(synapsessparse_10,meanfiringsparse_10_hipp,stdfiringsparse_10_hipp,'lineprops','-b');
hold on;
shadedErrorBar(synapsescluster_10,meanfiringcluster_10_hipp,stdfiringcluster_10_hipp,'lineprops','-m');
 title('Hippocampus');

subplot(2,2,2);

shadedErrorBar(synapsessparse_10,meanfiringsparse_10_pfc,stdfiringsparse_10_pfc,'lineprops','-b');
hold on;
 shadedErrorBar(synapsescluster_10,meanfiringcluster_10_pfc,stdfiringcluster_10_pfc,'lineprops','-m');

 title('PFC');
 
 subplot(2,2,3); 

 
 shadedErrorBar(synapsessparse_10,diamtwo_ia0_meanfiringsparse_10_hipp,diamtwo_ia0_stdfiringsparse_10_hipp,'lineprops','-b');
hold on;
shadedErrorBar(synapsescluster_10,diamtwo_ia0_meanfiringcluster_10_hipp,diamtwo_ia0_stdfiringcluster_10_hipp,'lineprops','-m');
 title('Hippocampus');

subplot(2,2,4);

shadedErrorBar(synapsessparse_10,diamtwo_ia0_meanfiringsparse_10_pfc,diamtwo_ia0_stdfiringsparse_10_pfc,'lineprops','-b');
hold on;
 shadedErrorBar(synapsescluster_10,diamtwo_ia0_meanfiringcluster_10_pfc,diamtwo_ia0_stdfiringcluster_10_pfc,'lineprops','-m');

 title('PFC');