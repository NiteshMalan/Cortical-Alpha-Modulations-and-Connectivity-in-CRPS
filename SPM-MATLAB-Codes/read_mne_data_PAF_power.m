% Reading data from Raghavan
% - Needs to have in MNE-Matlab in Matlab Path
% - Brainstorm needs to be running

%brainstorm
%% Original data
% Reading source estimate (`.stc` )
% https://mne.tools/stable/auto_tutorials/inverse/10_stc_class.html
sub_list = {'PASP001','PASP002','PASP003', 'PASP004', 'PASP005', 'PASP007', 'PASP008', 'PASP013', 'PASP015', 'PASP016', 'PASP017', 'PASP022', 'PASP023', 'PASP025', 'PASP026', 'PASP027', 'PASP028', 'PASP029', 'PASP030', 'PASP031', 'PASP032','PASP018', 'PASP006'};

%sub_list = {'PASP001'};

for sb=1: numel(sub_list)
    data_dir_lh = ['E:/GLM_SPM/PAF_Power/' sub_list{sb} '_PAF_power-lh.stc'];
    data_dir_rh = ['E:/GLM_SPM/PAF_Power/' sub_list{sb} '_PAF_power-rh.stc'];
    
    if isfile(data_dir_lh)
        stcL = mne_read_stc_file(data_dir_lh);        
    end
    if isfile(data_dir_lh)
        stcR = mne_read_stc_file(data_dir_rh);        
    end
    
    % Reading source space (`.fif`)
    % https://mne.tools/stable/generated/mne.SourceSpaces.html
    ssp = mne_read_source_spaces('E:/GLM_SPM/fsaverage-ico-5-src.fif');

    %% Importing in Brainstorm
    % Create Subject
    SubjectName = sub_list{sb};
    [~, iSubject] = db_add_subject(SubjectName, [], 0, 0);
    % Import the Anatomy that was used in MNE Python
    % In this case it was the FsAverage template
    sTemplates = bst_get('AnatomyDefaults');
    ix = find(~cellfun(@isempty,(regexpi({sTemplates.Name}, 'fsaverage'))));
    db_set_template(iSubject, sTemplates(ix), 0);
    [sSubject, iSubject] = bst_get('Subject', SubjectName);
    % Set WHITE MATTER surface with 327684 vertices (163842 per hemisphere) 
    % This is used as the source space in MNE-Python (fsaverage-ico-5-src.fif) is 
    % a downsample vertion of the surface with 327684 vertices (163842 per hemisphere) 
    iHdSurface = find(~cellfun(@isempty,(regexpi({sSubject.Surface.Comment}, '.*white_327684V'))));
    sSubject = db_surface_default(iSubject, 'Cortex', iHdSurface, 1);
    hdSurfaceFile = sSubject.Surface(sSubject.iCortex).FileName;
    
    %% Create a surface with 20k vertices (10k per hemisophere) based in the information
    % in the MNE-Python source space file (fsaverage-ico-5-src.fif)
    sHdSurfaceMat = in_tess_bst(hdSurfaceFile);
    % Create a copy to the map surface, this will be 20k vertices 
    ldSurfaceFileFull = strrep(file_fullpath(hdSurfaceFile), '_high.mat', '_low_20k.mat');
    copyfile(file_fullpath(hdSurfaceFile), ldSurfaceFileFull);
    ldSurfaceFile = file_short(ldSurfaceFileFull);
    % Edit vertices in the low resolution surface
    sLdSurfaceMat = in_tess_bst(ldSurfaceFile);
    sLdSurfaceMat.Comment = 'white_20484V';
    sLdSurfaceMat.VertConn = [];
    
    % Vertices from MNE-Python
    % 'stable' MUST BE USED TO KEEP THE SAME ORDER AS IN MNE-PYTHON FILE
    % as this is the order of source map
    [tmp, I, J] = intersect([ssp(1).rr; ssp(2).rr], ...  % Vertices 300k
                            [ssp(1).rr(logical(ssp(1).inuse), :); ssp(2).rr(logical(ssp(2).inuse), :)], ... % Vertices 20k
                            'rows', 'stable');
    % Vertices 
    sLdSurfaceMat.Vertices = sLdSurfaceMat.Vertices(I, :);
    % Index offset right hemisphere
    ixOffsetRight = sum(ssp(1).inuse);
    % Faces
    mneFaces = [ssp(1).use_tris; ssp(2).use_tris + ixOffsetRight];
    sLdSurfaceMat.Faces = mneFaces(:, [2,1,3]);
    % Normals 
    sLdSurfaceMat.VertNormals = sLdSurfaceMat.VertNormals(I,:);
    % Sphere registration
    sLdSurfaceMat.Reg.Sphere.Vertices = sLdSurfaceMat.Reg.Sphere.Vertices(I,:);
    % Curvature
    sLdSurfaceMat.Curvature = sLdSurfaceMat.Curvature(I,:);
    % SulciMap
    sLdSurfaceMat.SulciMap = sLdSurfaceMat.SulciMap(I,:);
    
    % Atlases 
    for iAtlas = 1:length(sLdSurfaceMat.Atlas)
        iRmScout = [];
        for iScout = 1:length(sLdSurfaceMat.Atlas(iAtlas).Scouts)
            % Replace the old vertices index with the new ones
            [a,b,c] = intersect(sLdSurfaceMat.Atlas(iAtlas).Scouts(iScout).Vertices, I);
            sLdSurfaceMat.Atlas(iAtlas).Scouts(iScout).Vertices = reshape(sort(J(c)), 1, []);
            % If scout has no vertex left: tag for deletion
            if isempty(sLdSurfaceMat.Atlas(iAtlas).Scouts(iScout).Vertices)
                iRmScout(end+1) = iScout;
            end
        end
        % Remove empty scouts
        if ~isempty(iRmScout)
            sLdSurfaceMat.Atlas(iAtlas).Scouts(iRmScout) = [];
        end
        % Set scouts seeds
        sLdSurfaceMat.Atlas(iAtlas).Scouts = panel_scout('SetScoutsSeed', sLdSurfaceMat.Atlas(iAtlas).Scouts, sLdSurfaceMat.Vertices);
    end
    % Save new 20k surface
    bst_save(ldSurfaceFileFull, sLdSurfaceMat, [], 1);
    iLdSurface = db_add_surface(iSubject, ldSurfaceFile, sLdSurfaceMat.Comment);
    sSubject = db_surface_default(iSubject, 'Cortex', iLdSurface, 1);
    
    %% Create a Brainstorm souce files and import it
    % Create Study
    iStudy = db_add_condition(sSubject.Name, 'PAF_power');
    % Create Source file 
    sResultMat = db_template('ResultsMat');
    sResultMat.Time          = [stcL.tmin, stcL.tmin];     % One sample
    sResultMat.Comment       = sub_list{sb}; 
    sResultMat.HeadModelType = 'surface';
    sResultMat.ImageGridAmp  = zeros(size(stcL.data,1) + size(stcR.data,1), length(sResultMat.Time)); % [VerticesLR, Time]
    ixOffsetRight = size(stcL.data,1);
    iVertices = double([stcL.vertices + 1; stcL.vertices + 1 + ixOffsetRight]); % +1 as Python is 0-indexed
    sResultMat.ImageGridAmp(:, 1) = [stcL.data(stcL.vertices + 1); stcR.data(stcR.vertices + 1)];
    sResultMat.ImageGridAmp(:, 2) = sResultMat.ImageGridAmp(iVertices, 1);
    % Save source file
    bst_save(bst_fullfile(pwd, 'results_surface_stc_as_bst_sources.mat'), sResultMat);
    
    %%
    % Import source file to database
    MapFile = import_sources(iStudy, ldSurfaceFile, bst_fullfile(pwd, 'results_surface_stc_as_bst_sources.mat'), [], 'BST', sub_list{sb});


end
% The source space seems to be WHITE matter of the the FsAverage 164k vertices (per hemisphere) 
% reduced to 20k vertices (10k per hemisphere)

