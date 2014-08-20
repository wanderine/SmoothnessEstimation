
%%

clear all
clc
close all

addpath('D:\Smoothness_estimation\Data\func');
addpath('D:\spm8');

% Loop over subjects
for file_number = 1:1
    
    tic
    
    if file_number < 10
        filename = ['D:\Smoothness_estimation\Data\func\subject_000' num2str(file_number) '\rest_' num2str(file_number) '.nii'];
        filedirectory = ['D:\Smoothness_estimation\Data\func\subject_000' num2str(file_number) '\'];
    elseif file_number < 100
        filename = ['D:\Smoothness_estimation\Data\func\subject_00' num2str(file_number) '\rest_' num2str(file_number) '.nii'];
        filedirectory = ['D:\Smoothness_estimation\Data\func\subject_00' num2str(file_number) '\'];
    elseif file_number < 1000
        filename = ['D:\Smoothness_estimation\Data\func\subject_0' num2str(file_number) '\rest_' num2str(file_number) '.nii'];
        filedirectory = ['D:\Smoothness_estimation\Data\subject_0' num2str(file_number) '\'];
    else
        filename = ['D:\Smoothness_estimation\Data\func\subject_' num2str(file_number) '\rest_' num2str(file_number) '.nii'];
        filedirectory = ['D:\Smoothness_estimation\Data\func\subject_' num2str(file_number) '\'];
    end
    
    V = spm_vol(filename);
    
    % Get number of time points
    st = size(V,1);
    st
    
    % Get repetition time
    a = V.private;
    TR = a.timing.tspace;
    
    
    %% Path containing data
    %--------------------------------------------------------------------------
    data_path = 'D:\Smoothness_estimation\Data\func_only_1\';
    data_path2 = 'D:\Smoothness_estimation\Data\';
    
    %% Initialise SPM defaults
    %--------------------------------------------------------------------------
    spm('Defaults','fMRI');
    
    spm_jobman('initcfg'); % useful in SPM8 only
    
    
    %% WORKING DIRECTORY (useful for .ps only)
    %--------------------------------------------------------------------------
    clear pjobs
    pjobs{1}.util{1}.cdir.directory = cellstr(data_path);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SPATIAL PREPROCESSING
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Select functional and structural scans
    %--------------------------------------------------------------------------
    
    % Select data
    if file_number < 10
        f = spm_select('FPList', fullfile(data_path2,['func\subject_000' num2str(file_number) '\']),'func*') ;
    elseif file_number < 100
        f = spm_select('FPList', fullfile(data_path2,['func\subject_00' num2str(file_number) '\']),'func*') ;
    elseif file_number < 1000
        f = spm_select('FPList', fullfile(data_path2,['func\subject_0' num2str(file_number) '\']),'func*') ;
    else
        f = spm_select('FPList', fullfile(data_path2,['func\subject_' num2str(file_number) '\']),'func*') ;
    end
    
    
    %% REALIGN
    %--------------------------------------------------------------------------
    pjobs{2}.spatial{1}.realign{1}.estwrite.data{1} = cellstr(f);
    
    %% SMOOTHING
    %--------------------------------------------------------------------------
    
    smoothings = 1:16;
    
    for smoothing = 1:16
        pjobs{2}.spatial{2 + smoothing - 1}.smooth.data = editfilenames(f,'prefix','r');
        pjobs{2}.spatial{2 + smoothing - 1}.smooth.fwhm = [smoothings(smoothing) smoothings(smoothing) smoothings(smoothing)];
        pjobs{2}.spatial{2 + smoothing - 1}.smooth.prefix = ['s' num2str(smoothings(smoothing))];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% FINISH
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    save('rest_batch_preprocessing.mat','pjobs');
    %spm_jobman('interactive',jobs); % open a GUI containing all the setup
    error1 = 0;
    try
        spm_jobman('run',pjobs);        % run preprocessing
    catch err
        err
        error1 = 1;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CLASSICAL STATISTICAL ANALYSIS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if error1 == 0
        
        for exp = 1:1
            for smoothing = 1:17
                
                clear jobs
                jobs{1}.util{1}.cdir.directory = cellstr(data_path);
                
                % Remove old SPM.mat
                spm_file = fullfile('D:\Smoothness_estimation\Data\func_only_1\classical','SPM.mat');
                if exist(spm_file,'file')==2
                    %system(['rm' spm_file]);
                    delete('D:\Smoothness_estimation\Data\func_only_1\classical\SPM.mat')
                end
               
                %% OUTPUT DIRECTORY
                %--------------------------------------------------------------------------
                jobs{1}.util{1}.md.basedir = cellstr(data_path);
                jobs{1}.util{1}.md.name = 'classical';
                
                %% MODEL SPECIFICATION AND ESTIMATION
                %--------------------------------------------------------------------------
                jobs{2}.stats{1}.fmri_spec.dir = cellstr(fullfile(data_path,'classical'));
                jobs{2}.stats{1}.fmri_spec.timing.units = 'secs';
                jobs{2}.stats{1}.fmri_spec.timing.RT = TR;
                if smoothing == 17
                    % Use unsmoothed data
                    jobs{2}.stats{1}.fmri_spec.sess.scans = editfilenames(f,'prefix','r');
                else
                    jobs{2}.stats{1}.fmri_spec.sess.scans = editfilenames(f,'prefix',['s' num2str(smoothings(smoothing)) 'r']);
                end
                jobs{2}.stats{1}.fmri_spec.sess.cond.name = 'active';
                %
                if exp == 1  % 10 s boxcar
                    jobs{2}.stats{1}.fmri_spec.sess.cond.onset = 10:20:(st*TR);
                    jobs{2}.stats{1}.fmri_spec.sess.cond.duration = 10;
                elseif exp == 2 % 15 s boxcar
                    jobs{2}.stats{1}.fmri_spec.sess.cond.onset = 15:30:(st*TR);
                    jobs{2}.stats{1}.fmri_spec.sess.cond.duration = 15;
                elseif exp == 3 % 20 s boxcar
                    jobs{2}.stats{1}.fmri_spec.sess.cond.onset = 20:40:(st*TR);
                    jobs{2}.stats{1}.fmri_spec.sess.cond.duration = 20;
                elseif exp == 4 % 30 s boxcar
                    jobs{2}.stats{1}.fmri_spec.sess.cond.onset = 30:60:(st*TR);
                    jobs{2}.stats{1}.fmri_spec.sess.cond.duration = 30;
                elseif exp == 5 % event 1 (2, 6)
                    jobs{2}.stats{1}.fmri_spec.sess.cond.onset = 8:8:(st*TR);
                    jobs{2}.stats{1}.fmri_spec.sess.cond.duration = 2;
                elseif exp == 6 % event 2 (4, 8)
                    jobs{2}.stats{1}.fmri_spec.sess.cond.onset = 12:12:(st*TR);
                    jobs{2}.stats{1}.fmri_spec.sess.cond.duration = 4;
                elseif exp == 7 % event 3 (rand 1)
                    
                    stt = st * TR;
                    last = 1;
                    tot = 0;
                    while 1
                        tot = onsets1(last) + durations1(last);
                        
                        if tot > stt
                            last = last - 1;
                            tot = onsets1(last) + durations1(last);
                            break
                        end
                        last = last + 1;
                        
                    end
                    
                    jobs{2}.stats{1}.fmri_spec.sess.cond.onset = onsets1(1:last);
                    jobs{2}.stats{1}.fmri_spec.sess.cond.duration = durations1(1:last);
                    
                elseif exp == 8 % event 4 (rand 2)
                    stt = st * TR;
                    last = 1;
                    tot = 0;
                    while 1
                        tot = onsets2(last) + durations2(last);
                        
                        if tot > stt
                            last = last - 1;
                            tot = onsets2(last) + durations2(last);
                            break
                        end
                        last = last + 1;
                        
                    end
                    
                    jobs{2}.stats{1}.fmri_spec.sess.cond.onset = onsets2(1:last);
                    jobs{2}.stats{1}.fmri_spec.sess.cond.duration = durations2(1:last);
                end
                
                % Use temporal derivatives
                jobs{2}.stats{1}.fmri_spec.bases.hrf.derivs = [1 0];
                
                % Use global normalization
                %jobs{2}.stats{1}.fmri_spec.global = 'Scaling';
                
                % Don't use global normalization
                jobs{2}.stats{1}.fmri_spec.global = 'None';
                
                % Add motion regressors
                jobs{2}.stats{1}.fmri_spec.sess.multi_reg = {[filedirectory 'rp_func_001.txt']};
                
                jobs{2}.stats{2}.fmri_est.spmmat = cellstr(fullfile(data_path,'classical','SPM.mat'));
                
                %% INFERENCE
                %--------------------------------------------------------------------------
                jobs{2}.stats{3}.con.spmmat = cellstr(fullfile(data_path,'classical','SPM.mat'));
                jobs{2}.stats{3}.con.consess{1}.tcon = struct('name','active > rest','convec', 1,'sessrep','none');
                
                jobs{2}.stats{4}.results.spmmat = cellstr(fullfile(data_path,'classical','SPM.mat'));
                jobs{2}.stats{4}.results.conspec.contrasts = Inf;
                jobs{2}.stats{4}.results.conspec.threshdesc = 'none';
                jobs{2}.stats{4}.results.conspec.thresh = 0.001; % for cluster based threshold
                jobs{2}.stats{4}.results.conspec.extent = 0;
                
                save('rest_batch_analysis.mat','jobs');
                
                error2 = 0;
                try
                    spm_jobman('run',jobs);        % run analysis
                catch err
                    err
                    error2 = 1;
                end
                
                % Move SPM.mat
                spm_file = fullfile('D:\Smoothness_estimation\Data\func_only_1\classical','SPM.mat');
                if exist(spm_file,'file')==2
                    
                    if smoothing == 17
                        % No smoothing
                        smoothmm = num2str(0);
                    else
                        smoothmm = num2str(smoothings(smoothing));
                    end
                    
                    if file_number < 10
                        movefile('D:\Smoothness_estimation\Data\func_only_1\classical\SPM.mat',['D:\Smoothness_estimation\Results\SPM_000' num2str(file_number) '_smoothing' smoothmm 'mm.mat'  ])
                    elseif file_number < 100
                        movefile('D:\Smoothness_estimation\Data\func_only_1\classical\SPM.mat',['D:\Smoothness_estimation\Results\SPM_00' num2str(file_number) '_smoothing' smoothmm 'mm.mat'  ])
                    elseif file_number < 1000
                        movefile('D:\Smoothness_estimation\Data\func_only_1\classical\SPM.mat',['D:\Smoothness_estimation\Results\SPM_0' num2str(file_number) '_smoothing' smoothmm 'mm.mat'  ])
                    else
                        movefile('D:\Smoothness_estimation\Data\func_only_1\classical\SPM.mat',['D:\Smoothness_estimation\Results\SPM_' num2str(file_number) '_smoothing' smoothmm 'mm.mat'  ])
                    end
                    
                end
                
            end
            
        end
        
    end
    
    file_number
    
    
    toc
    
end








