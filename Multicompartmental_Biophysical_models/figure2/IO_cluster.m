
% Investigating dendritic nonlinearities in a Realistic neocotical PV model,Alexandra Tzilivaki 2015 
%====================================================================================================





disp('Running mmatlab_job.m script...');

% Get the Task ID of this job
ID= getenv('LPARAMS');
%0okkk
my_vars = strsplit(ID, ' ');
area=my_vars(1);
cells=str2double(my_vars(2));
p=my_vars{3};
%f=str2double(my_vars(4));

maindir='/home/cluster/aleka/MainPath/Desktop/PV_Interneuron_Alexandra/Detailed_Alexandra/experiment/IO/';
cd(maindir);


if strcmp(area,'PFC');
    path_load=sprintf('cluster/PFC/PFC_%d/new/%s',cells,p);
    
     if ~exist(path_load)
      'File does not exist. Exiting now';
        exit
    end
else
path_load=sprintf('cluster/Hipp/Somogyi_%d/new/%s',cells,p);
 if ~exist(path_load)
      'File does not exist. Exiting now';
        exit
   
 end

end
cd(path_load);    
ndendrites =0:cells; %number of dendrites

% if strcmp(area,'Hipp');
    
    if cells==215
 l=[0 1 2 9 13 15 18 21 33 34 35 38 40 44 48 57 59 65 67 69 70 73 82 93 98 108 111 115 116 117 121 127 134 135 137 138 139 141 145 146 148 149 154 157 160 164 171 182 173 187 190 195 198 210]+1;  
    end
        if cells==62
     l=[0 5 11 12 13 14 15 16 21 27 40 4]+1;
        end
     if cells==57
        l=[0 3 9 13 26 41 48]+1;
     end
        if cells==199
        l=[0 6 40 66 89 90 91 92 111 112 126 161 179]+1;
        end
    if cells==71
        l=[0 1 14 15 19 27 33 32 37 41 47 48 54]+1;
    end
if cells==47
    l=[0 9 16 28 31 36 41]+1;
end

if cells==64
    l=[0 12 13 14 23 24 29 32 43 48 50 64]+1;
end

if cells==67
   l=[0 6 14 17 27 34 52 53 55 56 65]+1;
end



%if strcmp(p,'control2')
    
    vector=(1:1:20);
    synapses=1:1:20;
% else
% vector=(1:1:40); %vector used for expected rensponses.
% synapses =1:1:40; %number of synapses to each dendrite per run
% end

soma_voltage = {};
% expecteddend={};
% maxexpecteddend={};
% maxdendritesvoltage={};
%===============================================================================================================


for dends=1:length(ndendrites)
    
    for syns = 1:length(synapses)
        
        soma_voltage{syns,dends} = load(sprintf('recdendrite_%d_Synapses_%d.txt',ndendrites(dends),synapses(syns)));%dendritic voltage
        % calculate expected rensponses==============================
        if syns == 1
            
            expecteddend(syns,dends)=soma_voltage(syns,dends); %dendritic voltage when synapses =1
            maxexpecteddend_1(syns,dends)=max(cell2mat(expecteddend(syns,dends)));
            
        else
            
            maxexpecteddend_1(syns,dends)= -54 -  (-54 - (maxexpecteddend_1(1,dends))) .*(syns);
        end
        %calculate realrensponses====================================================
        
        maxdendritesvoltage_1(syns,dends)=max(soma_voltage{syns,dends});
        %======================================================================================================
        
    end
end

g=gt(maxdendritesvoltage_1,maxexpecteddend_1); %find supra       
g(:,l)=[];            %delete dendrites with false diameter!




 [sel, r] = max( g~=0, [], 1 ); %find first row of each column
 eval(sprintf('suprathreshold_%d_%s = r(r~=1)',cells,p)); 



cd(maindir);
if strcmp(area,'PFC');
savepath=('Data/PFC/Analysisnew/');
else
savepath=('Data/Hipp/Analysisnew/');
end
    cd(savepath);

filename=sprintf('results%d_%s.mat',cells,p);
%save(filename,'supra','maxdendritesvoltage_1','maxexpecteddend_1');
save(filename);