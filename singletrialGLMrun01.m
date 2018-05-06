mkdir outputs

%% step 1 generate alltrial regressor and noise regressor
load('template_AvsB.mat')
run('C:\Users\khanlab-admin\Documents\graham\load_tsv.m')
matlabbatch{1}.spm.stats.fmri_spec.sess.cond.onset = subBH14taskContRecogrun01events(:,1);
matlabbatch{1}.spm.stats.fmri_spec.sess.cond.duration = subBH14taskContRecogrun01events(1,2);
save('step1.mat','matlabbatch');
delete('SPM.mat');
spm_jobman('run',matlabbatch);
step1 = load('SPM.mat');
delete('SPM.mat');
alltrials = step1.SPM.xX.X(:,1);

run('C:\Users\khanlab-admin\Documents\graham\load_confounds.m');
noise = subBH14taskContRecogrun01boldconfounds(:,29:34); %6 motion parameters
noise = sum(noise,2);

regressors_matlabbatch = cell(length(subBH14taskContRecogrun01events),1);
for trial = 1:length(subBH14taskContRecogrun01events)
    tic
    matlabbatch = [];
    %% step 2 generate singletrialregressor
    load('template_AvsB.mat')
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond.onset = subBH14taskContRecogrun01events(trial,1);
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond.duration = subBH14taskContRecogrun01events(1,2);
    spm_jobman('run',matlabbatch);
    step2 = load('SPM.mat');
    delete('SPM.mat');
    singletrialregressor = step2.SPM.xX.X(:,1);

    %% step3 combine singletrialregressor and alltrialregressor+noise-singletrialregressor 
    load('template_AvsB.mat')
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1) = [];
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).name = 'singletrial';
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).val = singletrialregressor;
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).name = 'other';
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).val = alltrials + noise - singletrialregressor;
    spm_jobman('run',matlabbatch);
    
    %% step 4 estimate betas
    load('SPM.mat');
    delete('SPM.mat');
    spm_spm(SPM); % lets spm_spm this SPM!!!
    movefile('beta_0001.nii',sprintf('outputs/beta_0001_trial-%03d.nii',trial));
    regressors_matlabbatch{trial} = matlabbatch;
    matlabbatch = [];
    delete('SPM.mat');
    save('outputs/regressors_matlabbatch.mat','regressors_matlabbatch');
toc
end
%TODO: cleanup
