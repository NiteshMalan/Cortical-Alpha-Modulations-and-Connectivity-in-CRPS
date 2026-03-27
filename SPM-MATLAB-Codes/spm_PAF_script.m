% List of open inputs
nrun = X; % enter the number of runs here
jobfile = {'E:\gifty_from_BS\spm_PAF_script_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'EEG');
spm_jobman('run', jobs, inputs{:});
