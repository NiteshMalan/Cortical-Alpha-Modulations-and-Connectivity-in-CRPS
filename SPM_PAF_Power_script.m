% List of open inputs
nrun = 1; % enter the number of runs here
jobfile = {'C:\Users\malann\Documents\MATLAB\SPM_PAF_Power_script_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'EEG');
spm_jobman('run', jobs, inputs{:});
