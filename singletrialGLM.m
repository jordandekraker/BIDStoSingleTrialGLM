%updated 20180709, so far putting most output in temp dir
%non-BIDS part unfinished, for now only works for BIDS folder
%bids=input('is this BIDS (0 or 1)?');
bids=1;

%only works when you run the whole script
git_file=mfilename('fullpath');
s=dbstack();
folderpath=erase(git_file,s(1).name);
%% step 1 generate alltrial regressor and noise regressor

%load the template in the software path
load(strcat(folderpath,'template_AvsB.mat'));

%raw_dir according to sudesna's graham folder (pre-BIDS)
%bids_dir according to Jordan's graham folder (post-BIDS, post-fmriprep)
switch bids
    case 1
        %assume BIDS folder structure
        %read the BIDS event file
        bids_dir=input('BIDS directory:  ','s');%gets bids directory
        
        parent=strcat(bids_dir,'/derivatives/single_trial_GLM/');
        mkdir (parent,'output');
        mkdir (parent,'temp');
        output_dir=strcat(bids_dir,'/derivatives/single_trial_GLM/output/');
        temp_dir=strcat(bids_dir,'/derivatives/single_trial_GLM/temp/');
        copyfile strcat(folderpath,'template_AvsB.mat') temp_dir; %need to test when the whole script is complete
        %get all subject from BIDS base directory
        filekey=fullfile(bids_dir,'sub-*');
        file=dir(filekey);
        subject_id=extractfield(file,'name');
        
 
        %loop through subjects
        for i=1:length(subject_id)
            mkdir(output_dir,subject_id{i});
            mkdir(temp_dir,subject_id{i});
            subj_output=strcat(output_dir,subject_id{i},'/');
            subj_temp=strcat(temp_dir,subject_id{i},'/');
           
            %get #runs from each derivatives/fmriprep_1.0.7/fmriprep/subject/func/ folder
            runkey=fullfile(strcat(bids_dir,'/derivatives/fmriprep_1.0.7/fmriprep/',subject_id{i},'/func/'),'*run-*_bold_space-MNI152NLin2009cAsym_preproc.nii.gz');
            runfile=dir(runkey);
            sub(i).run=extractfield(runfile,'name');
            sub(i).id=subject_id{i};
            task=regexp(sub(i).run{1},'task-\w*_','match');%this will return something like "task-blablabla_"
            %unzip the nii.gz files into the temp directory
            gunzip(strcat(bids_dir,'/derivatives/fmriprep_1.0.7/fmriprep/',subject_id{i},'/func/',sub(i).run),strcat(temp_dir,subject_id{i}));
            %load the nii files, primarily to get the number of time points
            sub(i).runexp=spm_vol(strcat(temp_dir,subject_id{i},'/',erase(sub(i).run,'.gz')));
            
            
            %loop through runs for each participants, store output in the
            %structure named "sub"
            for j=1:length(sub(i).run)
                %make run-specific dir
                mkdir(subj_temp,strcat('run_', num2str(j)));
                run_temp=strcat(subj_temp,strcat('run_', num2str(j)));
                %change the output dir in matlabbatch to subject-specific and run-specific temp
                %dir
                matlabbatch{1}.spm.stats.fmri_spec.dir=cellstr(run_temp);
                run=regexp(sub(i).run{j},'run-\w*_bold_','match');%find each run to load the events.tsv
                sub(i).runevent{j}=load_tsv(bids_dir,subject_id{i},run{1},task{1});%store the loaded event files in sub.runevent
                sub(i).runconf{j}=load_confounds_2(bids_dir,subject_id{i},run{1},task{1});
                %% step 1 generate alltrialregressor (convolve with hrf)
                    %point the ...sess.scans to the correct file
                    %specify each time slice in the nii file
                    slice=(1:length(sub(i).runexp{j}));
                    slice=cellstr(num2str(slice'));
                    slice=cellfun(@strtrim,slice,'UniformOutput',false);%get rid of the white spaces
                    comma=repmat(',',length(sub(i).runexp{j}),1);
                    comma=cellstr(comma);
                    prefix={sub(i).runexp{j}(:).fname};
                    prefix=prefix';
                    sliceinfo=cellfun(@strcat,prefix,comma,slice,'UniformOutput',false);%use the first fname since they are all the same
                    matlabbatch{1}.spm.stats.fmri_spec.sess.scans=[];
                    matlabbatch{1}.spm.stats.fmri_spec.sess.scans=sliceinfo;
                matlabbatch{1}.spm.stats.fmri_spec.sess.cond.onset = cell2mat(sub(i).runevent{1,j}(:,1));
                matlabbatch{1}.spm.stats.fmri_spec.sess.cond.duration = cell2mat(sub(i).runevent{1,j}(1,2));
                save(strcat(run_temp,'step1.mat'),'matlabbatch');
                %delete('SPM.mat');   
                
                spm_jobman('run',matlabbatch);%should generate SPM.mat in subj_temp
                
                %% check point 20180609
                cd(run_temp); %start working in the subject-specific temp dir
                SPM1 = load('SPM.mat');
                delete('SPM.mat');
                alltrials = SPM1.SPM.xX.X(:,1);

                motion = cell2mat(sub(i).runconf{1, j}(:,30:35)); %6 motion parameters
                motion = sum(motion,2);

                regressors_matlabbatch = cell(length(sub(i).runevent{1,j}),1);
                %loop through trials
                for trial = 1:length(sub(i).runevent{1,j})
                    cd(run_temp); %SPM.swd will change working dir
                    tic
                    matlabbatch = [];
                    %% step 2 generate singletrialregressor (convolve with hrf)
                    load(strcat(run_temp,'step1.mat'))%load matlabbatch
                    matlabbatch{1}.spm.stats.fmri_spec.sess.cond.onset = sub(i).runevent{1,j}{trial,1};
                    matlabbatch{1}.spm.stats.fmri_spec.sess.cond.duration = sub(i).runevent{1,j}{trial,2};
                    spm_jobman('run',matlabbatch);%generate SPM.mat
                    SPM2 = load('SPM.mat');
                    delete('SPM.mat');
                    singletrialregressor = SPM2.SPM.xX.X(:,1);

                    %% step3 combine singletrialregressor and alltrialregressor+noise-singletrialregressor for one trial
                    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1) = [];
                    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).name = 'singletrial';
                    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).val = singletrialregressor;
                    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).name = 'other';
                    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).val = alltrials + motion - singletrialregressor;
                    spm_jobman('run',matlabbatch);

                    %% step 4 estimate betas for one trial
                    SPM3=load('SPM.mat');
                    delete('SPM.mat');                 
                    mkdir(run_temp,strcat('trial_', num2str(trial)));
                    SPM3.SPM.swd=strcat(subj_temp,strcat('trial_', num2str(trial)));%set output dir, need to change workign dir back in the next loop
                    spm_spm(SPM3.SPM); % lets spm_spm this SPM!!!
                    %movefile('beta_0001.nii',sprintf('outputs/beta_0001_trial-%03d.nii',trial));
                    regressors_matlabbatch{trial} = matlabbatch;        
                    %delete('SPM.mat');
                    save(strcat(run_temp,'regressors_matlabbatch.mat'),'regressors_matlabbatch');
                toc
                end

                
            end
        end 
    %% checkpoint 20180709
    %if not BIDS (**still unfinished**)
    case 0
        Subject_list=['BQ6','SH12','LH24','UO28','HE11','HE22','IS09','NL07','QK05','KX20','BH14','FO10','CD18'];
        Subject_number=['s2','s3','s4','s5','s6','s7','s8','s9','s10','s11','s12','s13','s14'];
        %xlsread the raw event file (onsets in TR)
        events=xlsread('C:\Users\haozi\Desktop\Anna_LSS\itemwise regressors\s1_r1_itemwise.xls');
end


