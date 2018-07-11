
subs = [2 3 4 5 6 7 8 9 10 11 12 13 14];
TR = 0.75;
dur = 1.2;
BIDSdir = '~/projects/rrg-akhanf/jdekrake/ANNA_OBJCAT_MTL_3T';

%% Import data from text file.
filename = 'sub_key.tsv';
delimiter = '\t';
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);
fclose(fileID);

subkey1 = [dataArray{1:end-1}];

% Clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans;

for s = subs
    for r = 1:7
        onset = xlsread(sprintf('~/projects/rrg-akhanf/sudesnac/single_trial_GLM/single_trial_GLM_basics/itemwise_regressors/s%d_r%d_itemwise.xls',s,r));
        sub = subkey1(find(subkey1==string(s))-1);
        onset = onset'/TR;
        
        duration = ones(length(onset),1) * dur;
        trial_type = cell(length(onset),1);
        
        T = table(onset,duration,trial_type);
        writetable(T,'test.txt','Delimiter','tab');
        system(sprintf('mv test.txt %s/%s/func/%s_task-ContRecog_run-%02d_events.tsv',BIDSdir,sub,sub,r));
    end
end
