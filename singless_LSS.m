function singless_LSS(bids_dir,fmriprep_foldername,output,sub,expstart_vol)



bids_dir = '/home/jordandekraker/graham/projects/rrg-akhanf/hyang336/Anna_LSS/BIDS_stable';
fmriprep_foldername = '/home/jordandekraker/graham/projects/rrg-akhanf/hyang336/Anna_LSS/BIDS_stable/derivatives/fmriprep_1.0.7/fmriprep';
output = 'test';
sub = 'sub-BH14';
expstart_vol = 1;

%updated 20190227, added another input variable to indicate the volume when the
%experiement starts (i.e. if there are 4 dummy scans, the experiment starts at the 5th
%TR/trigger/volume). In this version every participant in every run has to have the same number of
%dummy scans. 
%so far putting most output in temp dir

%20190603, added input var to specify fmriprep folder name,
%assuming it's in derivative under BIDS folder
sub_dir=strcat(output,'/singletrial_GLM/',sub);
%only works when you run the whole script
git_file=mfilename('fullpath');
s=dbstack();
folderpath=erase(git_file,s(1).name);
%% step 1 generate alltrial regressor and noise regressor
        %assume BIDS folder structure
        %temp_dir now under sub_dir
        mkdir (sub_dir,'output');
        mkdir (sub_dir,'temp');
        output_dir=strcat(sub_dir,'/output/');
        temp_dir=strcat(sub_dir,'/temp/');
        copyfile(strcat(folderpath,'template_AvsB.mat'), temp_dir); %copy the template matlabbatch into temp dir

            %get #runs from each derivatives/fmriprep_1.0.7/fmriprep/subject/func/ folder
            runkey=fullfile(strcat(bids_dir,'/derivatives/',fmriprep_foldername,'/fmriprep/',sub,'/func/'),'*run-*_bold_space-MNI152NLin2009cAsym_preproc.nii.gz');
            runfile=dir(runkey);
            substr=struct();
            substr.run=extractfield(runfile,'name');
            substr.id=sub;
            task=regexp(substr.run{1},'task-\w*_','match');%this will return something like "task-blablabla_"
            %unzip the nii.gz files into the temp directory
            gunzip(strcat(bids_dir,'/derivatives/',fmriprep_foldername,'/fmriprep/',sub,'/func/',substr.run),temp_dir);
            ninfo=niftiinfo(strcat(bids_dir,'/derivatives/',fmriprep_foldername,'/fmriprep/',sub,'/func/',substr.run));
            TR=ninfo.PixelDimensions(4)
            %load the nii files, primarily to get the number of time points
            substr.runexp=spm_vol(strcat(temp_dir,erase(substr.run,'.gz')));
            
            
            %loop through runs for each participants, store output in the
            %structure named "substr"
            for j=1:length(substr.run)
%% **moved into the run for-loop on 9/17/2018 14:41, did not change behavior**
                temp=struct(); 
                temp.matlabbatch = cell(1);
                %% load the template in the software path
                t=load(strcat(temp_dir,'template_AvsB.mat'));
                t.matlabbatch{1}.spm.stats.fmri_spec.timing.RT=TR;
                temp.matlabbatch{1} = t.matlabbatch{1};
                %make run-specific dir
                mkdir(temp_dir,strcat('run_', num2str(j)));
                run_temp=strcat(temp_dir,strcat('run_', num2str(j)));
                %change the output dir in matlabbatch to subject-specific and run-specific temp
                %dir
                temp.matlabbatch{1}.spm.stats.fmri_spec.dir=cellstr(run_temp);
                run=regexp(substr.run{j},'run-\d\d_','match');%find each run to load the events.tsv
                substr.runevent{j}=load_tsv(bids_dir,sub,run{1},task{1});%store the loaded event files in sub.runevent
                conf_name=strcat(bids_dir,'/derivatives/',fmriprep_foldername,'/fmriprep/',sub,'/func/',sub,'_',task{1},run{1},'bold_','confounds.tsv');%use run{1} since it's iteratively defined
                substr.runconf{j}=tdfread(conf_name,'tab');
                %% step 1 generate alltrialregressor (convolve with hrf)
                    %point the ...sess.scans to the correct file
                    %specify each time slice in the nii file
                    slice=(expstart_vol:length(substr.runexp{j}));
                    slice=cellstr(num2str(slice'));
                    slice=cellfun(@strtrim,slice,'UniformOutput',false);%get rid of the white spaces
                    comma=repmat(',',(length(substr.runexp{j})-expstart_vol+1),1);
                    comma=cellstr(comma);
                    prefix={substr.runexp{j}(expstart_vol:end).fname};
                    prefix=prefix';
                    sliceinfo=cellfun(@strcat,prefix,comma,slice,'UniformOutput',false);%use the first fname since they are all the same
                    temp.matlabbatch{1}.spm.stats.fmri_spec.sess.scans=[];
                    temp.matlabbatch{1}.spm.stats.fmri_spec.sess.scans=sliceinfo;
                temp.matlabbatch{1}.spm.stats.fmri_spec.sess.cond.onset = cell2mat(substr.runevent{1,j}(:,1));
                temp.matlabbatch{1}.spm.stats.fmri_spec.sess.cond.duration = cell2mat(substr.runevent{1,j}(1,2));
                save(strcat(run_temp,'/','step1.mat'),'temp');
                %delete('SPM.mat');   
                
                spm_jobman('run',temp.matlabbatch);%should generate SPM.mat in subj_temp
                
                %% check point 20180609
                %cd(run_temp); %start working in the subject-specific temp dir
                SPM1 = load(strcat(run_temp,'/','SPM.mat'));
                delete(strcat(run_temp,'/','SPM.mat'));
                alltrials = SPM1.SPM.xX.X(:,1);

                X = substr.runconf{j}.X; %6 motion parameters
                Y = substr.runconf{j}.Y;
                Z = substr.runconf{j}.Z;
                RotX = substr.runconf{j}.RotX;
                RotY = substr.runconf{j}.RotY;                
                RotZ = substr.runconf{j}.RotZ;
                motion = X+Y+Z+RotX+RotY+RotZ;

%                 regressors_matlabbatch = cell(length(sub.runevent{1,j}),1);
                %loop through trials
                
                for trial = 1:length(substr.runevent{1,j})
                    temp_trial=struct();
                    temp_trial.matlabbatch = cell(1);
                    %cd(run_temp); %SPM.swd will change working dir

%                     matlabbatch = [];
                    %% step 2 generate singletrialregressor (convolve with hrf)
                    t_trial=load(strcat(run_temp,'/','step1.mat'));%load matlabbatch
                    temp_trial.matlabbatch{1} = t_trial.temp.matlabbatch{1};
                    temp_trial.matlabbatch{1}.spm.stats.fmri_spec.sess.cond.onset = substr.runevent{1,j}{trial,1};
                    temp_trial.matlabbatch{1}.spm.stats.fmri_spec.sess.cond.duration = substr.runevent{1,j}{trial,2};
                    spm_jobman('run',temp_trial.matlabbatch);%generate SPM.mat
                    SPM2 = load(strcat(run_temp,'/','SPM.mat'));
                    delete(strcat(run_temp,'/','SPM.mat'));
                    singletrialregressor = SPM2.SPM.xX.X(:,1);

                    %% step3 combine singletrialregressor and alltrialregressor+noise-singletrialregressor for one trial
                    temp_trial.matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1) = [];
                    temp_trial.matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).name = 'singletrial';
                    temp_trial.matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).val = singletrialregressor;
                    temp_trial.matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).name = 'other';
                    temp_trial.matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).val = alltrials + motion(expstart_vol:end) - singletrialregressor;
                    spm_jobman('run',temp_trial.matlabbatch);

                    %% step 4 estimate betas for one trial
                    SPM3=load(strcat(run_temp,'/','SPM.mat'));
                    delete(strcat(run_temp,'/','SPM.mat'));                 
                    mkdir(run_temp,strcat('trial_', num2str(trial)));
                    SPM3.SPM.swd=strcat(run_temp,'/',strcat('trial_', num2str(trial)));%set output dir, need to change workign dir back in the next loop
                    spm_spm(SPM3.SPM); % lets spm_spm this SPM!!!
                    %movefile('beta_0001.nii',sprintf('outputs/beta_0001_trial-%03d.nii',trial));
%                     regressors_matlabbatch{trial} = temp_trial.matlabbatch;        
                    %delete(strcat(run_temp,'/','SPM.mat'));
%                     parsave(strcat(run_temp,'regressors_matlabbatch.mat'),regressors_matlabbatch);
                end

                
            end
        end 
