% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%                                                                        %
%                Temporal Filtering and Signal Regression                %
%         (global, white matter and cerebrospinal fluid signals)         %
%                                                                        %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% See the Anzar Abbas' preprocessing pipeline from 2016 for more info.
function filtering_signalreg(subjectdir, TR, fil, gsr, wmcsfr, ...
                             WM_probseg, CSF_probseg)
% Set the input file
input_FIL  = [subjectdir, '_smoothed.nii.gz'];

% Get information about subject and current time
k = strfind(subjectdir,'sub');  
sub = [subjectdir(k:k+5),' - ']; clear k;
time = whatsthetime();
warning('off','all');
%% ------------------------- Temporal Filtering ------------------------ %%
fprintf([time,' ',sub,'Temporal filtering ... '])
% Load data
nii = load_untouch_nii(input_FIL);
f = nii.img;
[X,Y,Z,T] = size(f);
% Getting the dimensions of the functional scan
f = reshape(f,[X*Y*Z,T]); 
f = f';
% Reshaping the scan into a 2D matrix of time over voxels
tc = f;
% Renaming the functional timeseries as tc
SF = 1/TR';
% Calculating the sampling frequency using the TR
x = size(tc,1);
% x is the number of timepoints
tc_mirror = tc(end:-1:1,:);
% This the mirror of the timeseries
tc_extended = [tc_mirror; tc; tc_mirror];
% Extend the timecourse in a periodical way
tc_fft = fft(tc_extended);
% Getting the fast fourier transform
if fil(1) > 0
    high_cutoff = round(size(tc_extended,1)*fil(1)/SF);
    tc_fft(1:high_cutoff,:) = 0;
end
% Conducting the high pass filtering
if fil(2) > 0
    low_cutoff = round(size(tc_extended,1)*fil(2)/SF);
    tc_fft(low_cutoff:end,:) = 0;
end
% Conducting the low pass filtering
tc_filtered = ifft(tc_fft);
% Inverse fast fourier transform
temp = 2*real(tc_filtered);
tc_filtered = temp(x+1:2*x,:);
% Getting the final filtered functional scan
f = tc_filtered';
f = reshape(f,[X,Y,Z,T]);
nii.img = f;
save_untouch_nii(nii,[subjectdir,'_fil.nii.gz']);
% Just to cut down on file size, I'm multiplying the functional
% scan by the brain mask again
fprintf('completed successfully.\n')
clear f1 tc_filtered temp tc_fft low_cutoff high_cutoff
clear tc_extended tc_mirror x

%% ---------------------- Global Signal Regression --------------------- %%
input_GSR  = [subjectdir, '_fil.nii.gz'];
time = whatsthetime();
fprintf([time,' ', sub, 'Global signal regression ... '])
if gsr == 1
    nii = load_untouch_nii(input_GSR);
    f = nii.img;
    % Loading the functional scan into Matlab
    f_gsr = f;
    % Predefining the global signal regressed functional scan
    f_mean = zeros(1,size(f,4));
    % Predefining what will me the mean of the functional scan
    for t = 1:size(f,4)
        temp = f(:,:,:,t);
        temp_sum = sum(temp(:));
        temp_length = length(find(temp(:)));
        temp_mean = temp_sum/temp_length;
        f_mean(t) = temp_mean;
    end
    % Taking the mean of the functional timeseries. This is the
    % functional timecourse.
    f_mean_zsc = zscore(f_mean);
    % Z-scoring the mean
    for x = 1:size(f,1)
        for y = 1:size(f,2)
            for z = 1:size(f,3)
                v = double(squeeze(f(x,y,z,:)));
                beta = (f_mean_zsc * f_mean_zsc') \ ...
                    (f_mean_zsc * v);
                f_gsr(x,y,z,:) = v' - f_mean_zsc * beta;
            end
        end
    end
    % Regressing the global signal from the functional
    % timeseries
    f = f_gsr; 
    nii.img = f;
    save_untouch_nii(nii,[subjectdir,'_gsr_fil.nii.gz']);
    % Saving the global signal regressed functional scan
    fprintf('completed successfully.\n')
else
    fprintf('Skipped\n')
end

%% --------------- White Matter and CSF Signal Regression -------------- %%
input_WMCSF  = [subjectdir, '_gsr_fil.nii.gz'];
time = whatsthetime();
fprintf([time, ' ', sub, 'WM/CSF signal regression ... '])
if wmcsfr == 1
    if gsr == 0
        nii = load_untouch_nii(input_WMCSF);
        f = nii.img;
        % Loading the functional scan
    end
    wm_mask = load_nii([subjectdir, '/anat/', WM_probseg]);
    wm_mask = wm_mask.img;
    % Loading white matter mask
    csf_mask = ...
        load_nii([subjectdir, '/anat/', CSF_probseg]);
    csf_mask = csf_mask.img;
    % Loading the CSF mask
    wmcsf_mask = double(logical(wm_mask + csf_mask));
    % Creating a mask of the white matter and CSF together and
    % making sure that it is a binary matrix
    f_wmcsf = double(f);
    % Predefining the matrix that will hold just the timeseries
    % from the white matter and CSF 
    for t = 1:size(f,4)
        f_wmcsf(:,:,:,t) = f_wmcsf(:,:,:,t) .* wmcsf_mask;
    end
    % Removing all the gray matter voxels from f_wmcsf
    f_wmcsf_r = f_wmcsf;
    % Predefining the matrix that will be the regressed signal
    % from the white matter and CSF
    f_wmcsf_mean = zeros(1,size(f,4));
    % Predefining what will me the mean of the functional scan
    for t = 1:size(f,4)
        temp = f_wmcsf(:,:,:,t);
        temp_sum = sum(temp(:));
        temp_length = length(find(temp(:)));
        temp_mean = temp_sum/temp_length;
        f_wmcsf_mean(t) = temp_mean;
    end
    % Taking the mean of the wmcsf signal. This is the
    % functional timecourse of the wmcsf.
    f_wmcsf_mean_zsc = zscore(f_wmcsf_mean);
    % Z-scoring the mean
    for x = 1:size(f,1)
        for y = 1:size(f,2)
            for z = 1:size(f,3)
                v = double(squeeze(f(x,y,z,:)));
                beta = (f_wmcsf_mean_zsc*f_wmcsf_mean_zsc') ...
                    \ (f_wmcsf_mean_zsc * v);
                f_wmcsf_r(x,y,z,:) = ...
                    v' - f_wmcsf_mean_zsc * beta;
            end
        end
    end
    % Regressing the white matter and CSF signal from the
    % functional timeseries 
    f = f_wmcsf_r;
    nii.img = f;
    save_untouch_nii(nii,[subjectdir,'_wmcsfr_gsr_fil.nii.gz']);
    % Saving the white matter/CSF regressef functional scan
    fprintf('completed successfully.\n')
else
    fprintf('Skipped\n')
end