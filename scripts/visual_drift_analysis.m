%% Load Neuropixels data
clear;clc;
warning('off');

% Define the path contining the neuropixels .mat files
neuropixels_results_path = 'D:\daniel-master\AllenBrainObservatory\Neuropixels\results\visual_dynamics_DD6\';

% Create a list of all the neuropixels .mat files to be loaded
% each .mat file corresponds to a single neuropixels recorded mouse
neuropixels_mouse_list = dir(neuropixels_results_path);
neuropixels_mouse_list(ismember({neuropixels_mouse_list.name}, {'.', '..'})) = [];
neuropixels_mouse_list = {neuropixels_mouse_list.name}; % should have 58 entries.

% Define the list of visual areas that will be loaded and analysed
brain_areas = {'VISp','VISl','VISal','VISpm','VISrl','VISam','LGd','LP'};

% Loop over the list of .mat files, load each mouse indevidually and bin the neuronal
% activty (e.i., calculate population vectors) for each of the presented visual stimuli.
for mouse = 1:length(neuropixels_mouse_list)
    clc;
    disp(['Loading Neuropixels data:'])
    disp(['Mouse: ',num2str(mouse),'/',num2str(length(neuropixels_mouse_list))])
    
    load([neuropixels_results_path,neuropixels_mouse_list{mouse}]) % load the .mat file
    
    % Store the number of movie repeats presented in each movie.
    % This variable will be used to subsample  mice based on their experimental groups.
    movie_repeats(mouse,:) = cellfun(@size,mean_running_speed_repeats,{1,1});
    
    % Store the number of cells recorded from each visual area.
    % This variable will be used to exclude areas with less than 15 recorded cells.
    neuropixels_cell_count(mouse,:) = cell_num(1,:);
    
    for nat_movie = 1:2
        repeats = [1:movie_repeats(mouse,nat_movie)]-1;
        
        % Define binning parameters for 'Natural movie 1' and 'Shuffled natural movie 1':
        if movie_repeats(mouse,nat_movie) == 10 || movie_repeats(mouse,nat_movie) == 30
            frames_rate = 30; % 30 frames per second
            movie_length = 30; % 30 seconds
            movie_frames = frames_rate * movie_length; % total of 900 frames
            bin_size = 30; % 1 bin = 30 frames = 1 second
            
            % Define bin ID for each frame in the movie:
            % e.g., frames 1-30 = 1, frames 31-60 = 2 ... frames 871-900 = 30
            binned_movie = nan(movie_frames,1);
            bin_edges = 1:bin_size:movie_frames;
            for bin = 1:length(bin_edges)
                binned_movie(bin_edges(bin):bin_edges(bin)+bin_size-1) = bin;
            end
            binned_movie_repeated = repmat(binned_movie,[movie_repeats(mouse,nat_movie),1]);
            
            % Bin the data into 90 bins instead of 30 bins for visualization only
            % (related to figures 6F, Figure S7A and Figure S7E):
            bin_size_tsne = 10; % 10 frames per second
            binned_movie_tsne = ones(movie_frames,1);
            bin_edges_tsne = 1:bin_size_tsne:movie_frames;
            for bin = 1:length(bin_edges_tsne)
                binned_movie_tsne(bin_edges_tsne(bin):bin_edges_tsne(bin)+bin_size-1) = bin;
            end
            binned_movie_repeated_tsne = repmat(binned_movie_tsne,[movie_repeats(mouse,nat_movie),1]);
            
            
            % Define binning parameters for 'Natural movie 3':
        elseif movie_repeats(mouse,nat_movie) == 5
            frames_rate = 30; % 30 frames per second
            movie_length = 120; % 120 seconds
            movie_frames = frames_rate * movie_length; % total of 3600 frames
            bin_size = 120; % 1 bin = 120 frames = 1 second
            
            % Define bin ID for each frame in the movie:
            % e.g., frames 1-30 = 1, frames 31-60 = 2 ... frames 3481-3600 = 30
            binned_movie = nan(movie_frames,1);
            bin_edges= 1:bin_size:movie_frames;
            for bin = 1:length(bin_edges)
                binned_movie(bin_edges(bin):bin_edges(bin)+bin_size-1) = bin;
            end
            binned_movie_repeated = repmat(binned_movie,[movie_repeats(mouse,nat_movie),1]);
        end
        
        % Bin the neuronal activty for each of the presented visual stimuli.
        for area = 1:length(brain_areas) % loop over brain areas
            if ~isempty(informative_rater_mat{nat_movie,area}) % check if are was recorded
                
                % Calculate population vecors based on the bins defined above:
                pop_vector_info_trials = []; % will contain the binned neuronal activit
                % rows - neuron ID
                % columns - Time bin ID
                % layers (3rd dim) - repeat ID
                
                sub = 1;
                for block = 1:2 % loop over blocks
                    for repeat = 1:movie_repeats(mouse,nat_movie) % loop over repeats
                        frames_temp = [1:movie_frames] + (movie_frames*repeats(repeat)); % frames ids of the spesific movie repeat
                        current_repeat = informative_rater_mat{nat_movie,area}(:,frames_temp,block); % subsample the neuronal activity of the specific movie repeat in the speficic block
                        
                        % Average the neuronal activity in the frames assigned to the same temporal bin
                        % (e.i., calculating the population vector of each temporal bin)
                        for bin = 1:length(bin_edges) % loop over time bins
                            pop_vector_info_trials(:,bin,sub) = mean(current_repeat(:,binned_movie_repeated(1:movie_frames) == bin),2,'omitnan');
                        end
                        sub = sub + 1;
                    end
                end
                
                % Store the binned neuronal data for each brain area of each mouse
                neuropixels_population_vectors{mouse,area,nat_movie} = pop_vector_info_trials;
                % rows - mouse ID
                % columns - brain area ID
                % layers (3rd dim) - movie ID
                
                % For tsne visualization only - 90 time bins instead of 30
                if nat_movie == 1 || movie_repeats(mouse,2) == 10
                    pop_vector_info_trials_tsne = [];
                    sub = 1;
                    for block = 1:2
                        for repeat = 1:movie_repeats(mouse,nat_movie)
                            frames_temp = [1:movie_frames] + (movie_frames*repeats(repeat));
                            current_repeat = informative_rater_mat{nat_movie,area}(:,frames_temp,block);
                            for bin = 1:length(bin_edges_tsne)
                                pop_vector_info_trials_tsne(:,bin,sub) = nanmean(current_repeat(:,binned_movie_repeated_tsne(1:movie_frames) == bin),2);
                            end
                            sub = sub + 1;
                        end
                    end
                    neuropixels_population_vectors_tsne{mouse,area,nat_movie} = pop_vector_info_trials_tsne;
                end
                
            else % n cases the iteraded brain area was not recorded define binned neuronal data as empty
                neuropixels_population_vectors{mouse,area,nat_movie} = [];
                if nat_movie == 1
                    neuropixels_population_vectors_tsne{mouse,area,nat_movie} = [];
                end
                
            end
        end
        
    end
    
    % Store for each mouse the average running speed in each movie repeat
    neuropixels_running_speed(mouse,:) = mean_running_speed_repeats; % Related to Figure S2A
    
    % Store for each mouse the average pupil size area in each movie repeat
    if ~isempty(mean_pupil_size_repeats)
        neuropixels_pupil_size(mouse,:) = mean_pupil_size_repeats; % Related to Figure S2C
    else
        neuropixels_pupil_size(mouse,:) = cell(1,2);
    end
    
    % Store for each mouse the binned neuronal actiivty in response to drifting gratings
    neuropixels_drifting_gratings(mouse,:) = valid_units_drifting_gratings; % Related to Figures S4A,F,G
end
clc;
disp(['Loading Neuropixels data:'])
disp('DONE!')

clearvars -except neuropixels_population_vectors neuropixels_drifting_gratings brain_areas...
    neuropixels_running_speed neuropixels_pupil_size neuropixels_cell_count movie_repeats...
    neuropixels_population_vectors_tsne

%% Load excitatory calcium imaging data
for area = 1:6
    % Define the path contining the excitatory calcium imaging .mat files
    results_path = ['E:\daniel_master\AllenBrainObservatory\calcium_imaging\results_files\excitatory4\',brain_areas{area},'\'];
    
    % Create a list of all the excitatory calcium imaging .mat files to be loaded
    % each .mat file corresponds to a single calcium imaging recorded mouse
    mat_list = dir([results_path,'*.mat']);
    mat_list = {mat_list.name};
    
    % Define empty variable that will store data across all mice from the
    % same brain area:
    calcium_population_vectors_across_mice = {}; % Binned activity (events) during natural movies
    calcium_drifting_gratings_across_mice = {}; % Binned activity (events) during drifting gratings
    calcium_spont_population_vectors_across_mice = {}; % Binned activity (events) during spontanous activity
    raw_calcium_population_vectors_across_mice = {}; % Binned activity (raw df/f) during natural movies
    mean_cell_num_across_mice = []; % Number of recorded cells for each mouse in each area
    binned_running_speed_across_mice = {}; % Running speed across days during natural movie 1
    binned_pupil_size_across_mice = {}; % Pupil size across days during natural movie 1
    imaging_depth_all_mice = []; % Imaging depth for each mouse
    sorted_mouse_age = []; % Age in each imaging session for each mouse
    
    for file = 1:length(mat_list) % loop over mice
        clc;
        disp(['Loading calcium imaging excitatory data:'])
        disp(['Area: ',num2str(area),'\',num2str(6)])
        disp(['Mouse: ',num2str(file),'\',num2str(length(mat_list))])
        load([results_path,mat_list{file}]) % load .mat file for current mouse
        
        % Define binning parameters for 'Natural movie 1':
        frames_rate = 30; % imaging frame rate (30 Hz).
        repeats_movie1 = 30; % number of natural movie repeats across days (10 repeats in each session x 3 sessions = 30 movie repeats)
        movie1_length = 30; % length of natural movie 1 (30 seconds)
        movie1_frames = frames_rate * movie1_length; % number frames in a single natural movie 1 repeat (900 frames)
        movie1_bin_size = 30; % number of frames in each bin of neuronal activity (1 bin = 30 frames = 1 sec)
        
        % Define bin ID for each frame in the natural movie 1:
        % e.g., frames 1-30 = 1, frames 31-60 = 2 ... frames 3481-3600 = 30
        binned_movie1 = ones(movie1_frames,1);
        movie1_bin_edges = 1:movie1_bin_size:movie1_frames;
        for bin = 1:length(movie1_bin_edges)
            binned_movie1(movie1_bin_edges(bin):movie1_bin_edges(bin)+movie1_bin_size-1) = bin;
        end
        binned_movie_repeated1 = repmat(binned_movie1,[repeats_movie1,1]);
        
        % Bin the data into 90 bins instead of 30 bins for visualization only
        % (related to figures 6F, Figure S7A and Figure S7E):
        frames_rate = 30;
        repeats_movie2 = 30;
        movie2_length = 30;
        movie2_frames = frames_rate * movie2_length;
        movie2_bin_size = 10; % 10 frames per second
        binned_movie2 = ones(movie2_frames,1);
        movie2_bin_edges = 1:movie2_bin_size:movie2_frames;
        for bin = 1:length(movie2_bin_edges)
            binned_movie2(movie2_bin_edges(bin):movie2_bin_edges(bin)+movie2_bin_size-1) = bin;
        end
        binned_movie_repeated2 = repmat(binned_movie2,[repeats_movie2,1]);
        
        % Define binning parameters for 'Natural movie 3':
        frames_rate = 30; % imaging frame rate (30 Hz).
        repeats_movie3 = 10;
        movie3_length = 120;
        movie3_frames = frames_rate * movie3_length;
        movie3_bin_size = 4*frames_rate;
        binned_movie3 = ones(movie3_frames,1);
        movie3_bin_edges = 1:movie3_bin_size:movie3_frames;
        for bin = 1:length(movie3_bin_edges)
            binned_movie3(movie3_bin_edges(bin):movie3_bin_edges(bin)+movie3_bin_size-1) = bin;
        end
        binned_movie_repeated3 = repmat(binned_movie3,[repeats_movie3,1]);
        
        % Define binning parameters for NM3 VS DG decoding analysis
        % (Related to Fig. S7G):
        frames_rate = 30; % 10 frames per second
        repeats_movie3_decode = 10;
        movie3_bin_size_decode = 90;
        binned_movie3_decode = ones(movie3_frames,1);
        movie3_bin_edges_decode = 1:movie3_bin_size_decode:movie3_frames;
        for bin = 1:length(movie3_bin_edges_decode)
            binned_movie3_decode(movie3_bin_edges_decode(bin):movie3_bin_edges_decode(bin)+movie3_bin_size_decode-1) = bin;
        end
        binned_movie_repeated3_decode = repmat(binned_movie3_decode,[repeats_movie3_decode,1]);
        
        % store bin information for each of the natural movies
        repeats_movie = [repeats_movie1,repeats_movie3,repeats_movie2,repeats_movie3];
        movie_frames = [movie1_frames,movie3_frames,movie2_frames,movie3_frames];
        movie_bin_edges = {movie1_bin_edges,movie3_bin_edges,movie2_bin_edges,movie3_bin_edges_decode}; 
        binned_movie_repeated = {binned_movie_repeated1,binned_movie_repeated3,binned_movie_repeated2,binned_movie_repeated3_decode};
        
        % load neuronal responses for different natural movies 
        natural_movie1_traces = reshape(united_traces_days_events,size(united_traces_days_events,1),27000); % natural movie 1
        natural_movie3_traces = filtered_traces_days_events{1,2}; % natural movie 3
        natural_movie_traces = {natural_movie1_traces,natural_movie3_traces,natural_movie1_traces,natural_movie3_traces}; % store neuronal responses of both movies
        
        
        spont_blocks_repeats = 20;
        spont_length = 30;
        spont_frames = frames_rate * spont_length;
        spont_bin_size = 30;
        binned_spont_blocks = ones(spont_frames,1);
        spont_bin_edges = 1:spont_bin_size:spont_frames;
        for bin = 1:length(spont_bin_edges)
            binned_spont_blocks(spont_bin_edges(bin):spont_bin_edges(bin)+spont_bin_size-1) = bin;
        end
        binned_spont_blocks_repeated = repmat(binned_spont_blocks,[spont_blocks_repeats,1]);
        
        repeats_spont = [repeats_movie1,spont_blocks_repeats];
        spont_frames = [spont_frames,spont_frames];
        bin_edges_spont = {movie1_bin_edges,spont_bin_edges};
        binned_spont_repeated = {binned_movie_repeated1,binned_spont_blocks_repeated};
        
        spont_traces_days = reshape(united_traces_days_spont_events,size(united_traces_days_spont_events,1),27000);
        spont_traces_blocks = reshape(filtered_traces_days_events{3,3},size(filtered_traces_days_events{3,3},1),18000);
        spont_traces = {spont_traces_days,spont_traces_blocks};
        
        
        natural_movie1_running = cell2mat(natural_movie_running_sorted(:,1))';
        natural_movie1_pupil= cell2mat(natural_movie_pupil_sorted(:,1))';
        
        % binning procedure for natural movies and behavioural measurments
        for nat_movie = 1:4
            pop_vector_info_trials = [];
            binned_running_speed = [];
            binned_pupil_size = [];
            sub = 1;
            for repeat = 1:repeats_movie(nat_movie)
                frames_temp = [1:movie_frames(nat_movie)] + (movie_frames(nat_movie)*(repeat-1));
                current_repeat = natural_movie_traces{nat_movie}(:,frames_temp);
                if nat_movie ==1
                    current_repeat_running = natural_movie1_running(:,frames_temp);
                    current_repeat_pupil = natural_movie1_pupil(:,frames_temp);
                end
                for bin = 1:length(movie_bin_edges{nat_movie})
                    pop_vector_info_trials(:,bin,sub) = nanmean(current_repeat(:,binned_movie_repeated{nat_movie}(frames_temp) == bin),2);
                    if nat_movie ==1
                        binned_running_speed(:,bin,sub)  = nanmean(current_repeat_running(:,binned_movie_repeated{nat_movie}(frames_temp) == bin),2);
                        binned_pupil_size(:,bin,sub)  = nanmean(current_repeat_pupil(:,binned_movie_repeated{nat_movie}(frames_temp) == bin),2);
                    end
                end
                sub = sub + 1;
            end
            calcium_population_vectors_across_mice(file,nat_movie) = {pop_vector_info_trials};
            binned_running_speed_across_mice(file,nat_movie) = {binned_running_speed};
            binned_pupil_size_across_mice(file,nat_movie) = {binned_pupil_size};
        end
        
        % spont activity
        for nat_movie = 1:2
            pop_vector_info_trials = [];
            sub = 1;
            for repeat = 1:repeats_spont(nat_movie)
                frames_temp = [1:spont_frames(nat_movie)] + (spont_frames(nat_movie)*(repeat-1));
                current_repeat =spont_traces{nat_movie}(:,frames_temp);
                for bin = 1:length(bin_edges_spont{nat_movie})
                    pop_vector_info_trials(:,bin,sub) = nanmean(current_repeat(:, binned_spont_repeated{nat_movie}(frames_temp) == bin),2);
                end
                sub = sub + 1;
            end
            calcium_spont_population_vectors_across_mice(file,nat_movie) = {pop_vector_info_trials};
        end
        
        
        % binned activity raw df/f
        raw_calcium_population_vectors_across_mice(file) = {raw_pop_vector_info_trials};
        calcium_drifting_gratings_across_mice(file) = filtered_traces_days_events(1,4);
        
        mean_cell_num_across_mice(file,:) = cellfun(@length,cell_registration);
        imaging_depth_all_mice(file) = imaging_depth;
        sorted_mouse_age(file,:) = sort(mouse_age);
    end
    calcium_excitatory_population_vectors{area} = calcium_population_vectors_across_mice;
    calcium_excitatory_drifting_gratings{area} =  calcium_drifting_gratings_across_mice;
    calcium_excitatory_spont_population_vectors{area} = calcium_spont_population_vectors_across_mice;
    calcium_excitatory_cell_count{area} =  mean_cell_num_across_mice;
    calcium_excitatory_imaging_depth{area} = imaging_depth_all_mice;
    calcium_excitatory_running_speed{area} = binned_running_speed_across_mice;
    calcium_excitatory_pupil_size{area} = binned_pupil_size_across_mice;
    calcium_excitatory_population_vectors_raw{area} = raw_calcium_population_vectors_across_mice;
    calcium_excitatory_sorted_mouse_age{area} = sorted_mouse_age;
end

clc;
disp(['Loading calcium imaging excitatory data:'])
disp('DONE!')

clearvars -except neuropixels_population_vectors neuropixels_drifting_gratings brain_areas...
    neuropixels_running_speed neuropixels_pupil_size neuropixels_cell_count movie_repeats...
    calcium_excitatory_population_vectors calcium_excitatory_drifting_gratings...
    calcium_excitatory_spont_population_vectors calcium_excitatory_cell_count...
    calcium_excitatory_imaging_depth calcium_excitatory_running_speed ...
    calcium_excitatory_pupil_size calcium_excitatory_population_vectors_raw...
    calcium_excitatory_sorted_mouse_age neuropixels_population_vectors_tsne


%% Load inhibitory calcium imaging data

for area = [1,2,4]
    results_path = ['E:\daniel_master\AllenBrainObservatory\calcium_imaging\results_files\inhibitory\',brain_areas{area},'\'];
    mat_list = dir([results_path,'*.mat']);
    mat_list = {mat_list.name};
    calcium_population_vectors_across_mice = {};
    mean_cell_num_across_mice = [];
    cre_line_across_mice = [];
    for file =1:length(mat_list)
        clc;
        disp(['Loading calcium imaging inhibitory data:'])
        disp(['Area: ',num2str(area),'\',num2str(3)])
        disp(['Mouse: ',num2str(file),'\',num2str(length(mat_list))])
        load([results_path,mat_list{file}])
        
        if ~isempty(strfind(cre_line,'Sst') == 1)
            cre_line_across_mice(file) = 1; %sst
        elseif ~isempty(strfind(cre_line,'Vip') == 1)
            cre_line_across_mice(file) = 2; %vip
        elseif ~isempty(strfind(cre_line,'Pvalb') == 1)
            cre_line_across_mice(file) = 3; %pvalb
        end
        
        frames_rate = 30;
        repeats_movie1 = 30;
        movie1_length = 30;
        movie1_frames = frames_rate * movie1_length;
        movie1_bin_size = 30;
        binned_movie1 = ones(movie1_frames,1);
        movie1_bin_edges = 1:movie1_bin_size:movie1_frames;
        for bin = 1:length(movie1_bin_edges)
            binned_movie1(movie1_bin_edges(bin):movie1_bin_edges(bin)+movie1_bin_size-1) = bin;
        end
        binned_movie_repeated1 = repmat(binned_movie1,[repeats_movie1,1]);
        
        
        repeats_movie3 = 10;
        movie3_length = 120;
        movie3_frames = frames_rate * movie3_length;
        movie3_bin_size = 120;
        binned_movie3 = ones(movie3_frames,1);
        movie3_bin_edges = 1:movie3_bin_size:movie3_frames;
        for bin = 1:length(movie3_bin_edges)
            binned_movie3(movie3_bin_edges(bin):movie3_bin_edges(bin)+movie3_bin_size-1) = bin;
        end
        binned_movie_repeated3 = repmat(binned_movie3,[repeats_movie3,1]);
        
        repeats_movie = [repeats_movie1,repeats_movie3];
        movie_frames = [movie1_frames,movie3_frames];
        movie_bin_edges = {movie1_bin_edges,movie3_bin_edges};
        binned_movie_repeated = {binned_movie_repeated1,binned_movie_repeated3};
        natural_movie1_traces = reshape(united_traces_days_events,size(united_traces_days_events,1),27000);
        natural_movie3_traces = filtered_traces_days_events{1,2};
        natural_movie_traces = {natural_movie1_traces,natural_movie3_traces};
        
        natural_movie1_running = cell2mat(natural_movie_running_sorted(:,1)')';
        natural_movie3_running = natural_movie_running{1,2}';
        natural_movie_running = {natural_movie1_running,natural_movie3_running};
        
        
        
        for nat_movie = 1:2
            pop_vector_info_trials = [];
            sub = 1;
            for repeat = 1:repeats_movie(nat_movie)
                frames_temp = [1:movie_frames(nat_movie)] + (movie_frames(nat_movie)*(repeat-1));
                current_repeat = natural_movie_traces{nat_movie}(:,frames_temp);
                for bin = 1:length(movie_bin_edges{nat_movie})
                    pop_vector_info_trials(:,bin,sub) = nanmean(current_repeat(:,binned_movie_repeated{nat_movie}(frames_temp) == bin),2);
                end
                sub = sub + 1;
            end
            calcium_population_vectors_across_mice(file,nat_movie) = {pop_vector_info_trials};
        end
        
        mean_cell_num_across_mice(file,:) = cellfun(@length,cell_registration);
        
        
        
    end
    calcium_inhibitory_population_vectors{area} = calcium_population_vectors_across_mice;
    calcium_inhibitory_cell_count{area} =  mean_cell_num_across_mice;
    calcium_inhibitory_cre_line{area} =  cre_line_across_mice;
end


clc;
disp(['Loading calcium imaging inhibitory data:'])
disp('DONE!')

clearvars -except neuropixels_population_vectors neuropixels_drifting_gratings brain_areas...
    neuropixels_running_speed neuropixels_pupil_size neuropixels_cell_count movie_repeats...
    calcium_excitatory_population_vectors calcium_excitatory_drifting_gratings...
    calcium_excitatory_spont_population_vectors calcium_excitatory_cell_count...
    calcium_excitatory_imaging_depth calcium_excitatory_running_speed ...
    calcium_excitatory_pupil_size calcium_excitatory_population_vectors_raw...
    calcium_inhibitory_population_vectors calcium_inhibitory_cell_count...
    calcium_inhibitory_cre_line calcium_excitatory_sorted_mouse_age...
    neuropixels_population_vectors_tsne

brain_areas = {'V1','LM','AL','PM','RL','AM','dLGN','LP','CA1','CA2','CA3','DG'};

%% Define color schemes and load colormaps
colors = [0 0.7 0.8 ;0 0.7 0.6; 0.9 0.8 0.2; 0.9 0.6 0.2; 0.9 0.5 0.7; 1 0.3 0.4;...
    0.4 0.7 0.8; 0.4 0.7 0.6; 0.5 0.5 0.5; 0.4 0.4 0.4; 0.3 0.3 0.3; 0.1 0.1 0.1];
colors2 = [0 0.5 0.6 ;0 0.5 0.4; 0.7 0.6 0; 0.7 0.4 0; 0.7 0.3 0.5; 0.8 0.1 0.2;...
    0.5 0.5 0.5; 0.4 0.4 0.4; 0.3 0.3 0.3; 0.2 0.2 0.2; 0.1 0.1 0.1; 0.0 0.0 0.0];

load('D:\daniel-master\AllenBrainObservatory\Neuropixels\rank_colormap.mat')
load('D:\daniel-master\AllenBrainObservatory\Neuropixels\stat_colormap.mat')
load('D:\daniel-master\AllenBrainObservatory\Advanced\analysis\newmap3.mat')
load('D:\daniel-master\AllenBrainObservatory\Neuropixels\magma_colormap.mat')
load('D:\daniel-master\AllenBrainObservatory\Neuropixels\sig_colormap.mat')
load('D:\daniel-master\AllenBrainObservatory\Neuropixels\new_jet_colormap.mat')

%% Figure 1B - calcium imaging excitatory Cre lines cell counts
% creates a boxplot visualizing the distribution of the average number of cells
% across the three imaging sessions for each mouse in each visual area:

average_cell_count = nan(94,6); % define an empty NaN matrix (94 mice x 6 brain areas)
for area = 1:6 % loop over brain areas
    % current_area - matrix in the size of Mx3 (number of mice x 3 imaging sessions)
    % each entry is the number of recorded cells for each mouse (rows) in each
    % imaging session (columns):
    current_area = calcium_excitatory_cell_count{area};
    average_current_area = mean(current_area,2,'omitnan'); % mean number of cells across sessions for each mouse
    num_mice(area) = length(current_area); % number of recorded mice
    average_cell_count(1:num_mice(area),area) = average_current_area;
end

% visualization of cell counts distributions
figure
figure_boxplot(average_cell_count(:,1:6))
set(gca,'xticklabels',brain_areas(1:6))
xtickangle(90)
ylabel('Cell count')

%% Figure 1D - Single cell tuning examples - seconds 3,9,14
% visualize the neuronal responses of three exmaple calcium imaging cells
% during the presentation of 'Natural movie 1'.

nat_movie = 1; % natural movie 1
area_list = [6,2,2]; % area AM, area AL, area AL
mouse_list = [6,57,44]; % mouse #6, #57 and #44 from each area respectively
cell_list = [40,21,78]; % cells #40, #21 and #78  from each mouse respectively

for current_cell = 1:3 % loop over cells
    % example_cell - binned activity rate (calcium events) for example cell
    % (30 movie repeats x 30 time bins)
    example_cell = squeeze(calcium_excitatory_population_vectors{area_list(current_cell)}{mouse_list(current_cell),nat_movie}(cell_list(current_cell),:,:)*30)';
    
    % gaussian smoothing of neuronal activity for each movie repeat
    smooth_example_cell = [];
    for repeat = 1:size(example_cell,1) % loop over movie repeats
        smooth_example_cell(repeat,:) = imgaussfilt(example_cell(repeat,:),2); % sigma = 2;
    end
    
    % normalize neuronal activity for each movie repeat based on peak firing rate
    norm_example_cell = smooth_example_cell ./ max(smooth_example_cell,[],2);
    
    % visualize normalized activity patterns for each cell
    subplot(2,3,current_cell)
    imagesc(norm_example_cell)
    hold on
    plot(xlim, [10 10]+0.5,'--','linewidth',2,'color','w')
    plot(xlim, [20 20]+0.5,'--','linewidth',2,'color','w')
    title(['Cell #',num2str(current_cell)])
    if current_cell == 1
        ylabel('Movie repeat')
        text(0.425,0.75,'Session 1','Units','normalized','color','w')
        text(0.425,0.4,'Session 2','Units','normalized','color','w')
        text(0.425,0.1,'Session 3','Units','normalized','color','w')
    end
    colormap(newmap3)
    
    
    % calculate mean activity rate across movie repeats for each session
    mean_example_cell = [];
    mean_example_cell(1,:) = mean(example_cell(1:10,:),1,'omitnan'); % average activity for session 1 (repeats 1-10)
    mean_example_cell(2,:) = mean(example_cell(11:20,:),1,'omitnan'); % average activity for session 2 (repeats 11-20)
    mean_example_cell(3,:) = mean(example_cell(21:30,:),1,'omitnan'); % average activity for session 3 (repeats 21-30)
    
    % visualize average activity patterns (tuning curve) for session of each cell
    subplot(2,3,current_cell+3)
    hold on
    plot(mean_example_cell(1,:),'-','markersize',20,'linewidth',1.5,'color',[0.4 0.4 0.4]) % session 1
    plot(mean_example_cell(2,:)','-','markersize',20,'linewidth',1.5,'color',[0.6 0.6 0.6]) % session 2
    plot(mean_example_cell(3,:)','-','markersize',20,'linewidth',1.5,'color',[0.8 0.8 0.8]) % session 3
    
    xlim([1 30])
    ylim([0 18])
    if current_cell == 1
        ylabel('Mean activity rate (events/sec)')
        legend({'Session 1','Session 2','Session 3'})
        legend('boxoff')
    elseif current_cell == 2
        xlabel('Time in movie (sec)')
    end
    
end

%% Figure 1F - neuropixels cell counts
% creates a boxplot visualizing the distribution of the number Neuropixels
% recorded units from mouse in each visual area:

% valid_cell_counts - matrix containg the number of recorded units from
% each mouse and visual area
valid_cell_counts = neuropixels_cell_count(:,:,1);
valid_cell_counts(valid_cell_counts==0) = NaN; % conver 0 values into NaNs
sum(valid_cell_counts(:,1:6)>0)
% visualization of cell counts distributions
figure
figure_boxplot(valid_cell_counts(:,1:6)) % visualize only visual cortical areas (columns 1-6)
set(gca,'xticklabels',brain_areas(1:6))
xtickangle(90)
ylabel('Cell count')

%% Figure 1H -  Single cell tuning examples - seconds 19,24,29
% visualize the neuronal responses of three exmaple Neuropixels units
% during the presentation of 'Natural movie 1'.

nat_movie = 1;
area_list = [2,1,2]; % area AL, area V1, area AL
mouse_list = [53,54,58]; % mouse #53, #54 and #48 from each area respectively
unit_list = [22,39,7]; % units #22, #39 and 7  from each mouse respectively

for unit = 1:3 % loop over units
    % current_units - binned activity rate (spikes) for example unit
    % (60 movie repeats x 30 time bins)
    current_unit = squeeze(neuropixels_population_vectors{mouse_list(unit),area_list(unit),nat_movie}(unit_list(unit),:,:)*30)';
    
    % gaussian smoothing of neuronal activity for each movie repeat
    smooth_current_unit = [];
    for repeat = 1:size(current_unit,1)
        smooth_current_unit(repeat,:) = imgaussfilt(current_unit(repeat,:),2); % sigma = 2;
    end
    
    % normalize neuronal activity for each movie repeat based on peak firing rate
    norm_current_unit = smooth_current_unit ./ max(smooth_current_unit,[],2);
    
    % visualize normalized activity patterns for each unit
    subplot(2,3,unit)
    imagesc(norm_current_unit)
    hold on
    title(['Unit #',num2str(unit)])
    plot(xlim, [30 30]+0.5,...
        'linewidth',2,'color','w')
    
    if unit == 1
        ylabel('Movie repeat')
    elseif unit == 3
        text(0.1,0.9,'Block A','Units','normalized','color','w')
        text(0.1,0.4,'Block B','Units','normalized','color','w')
    end
    colormap(newmap3)
    
    % calculate mean activity rate across movie repeats for each block
    mean_current_unit = [];
    mean_current_unit(1,:) = mean(current_unit(1:30,:),1,'omitnan');
    mean_current_unit(2,:) = mean(current_unit(31:60,:),1,'omitnan');
    
    % visualize average activity patterns (tuning curve) for session of each unit
    subplot(2,3,unit+3)
    hold on
    plot(mean_current_unit(1,:),'-','markersize',20,'linewidth',1.5,'color',[0.4 0.4 0.4])
    plot(mean_current_unit(2,:)','-','markersize',20,'linewidth',1.5,'color',[0.6 0.6 0.6])
    
    xlim([1 30])
    if unit == 1
        ylabel('Mean activity rate (spike/sec)')
        
    elseif unit == 2
        xlabel('Time in movie (sec)')
    elseif unit == 3
        legend({'Block A','Block B'})
        legend('boxoff')
    end
    
end

%% Figure 2A - PV correlation across repeats for a single mouse recorded from area PM
% Main panel: population vector (PV) correlation between the first 10 (out of 30) natural movie 1 repeats
% of the first block, recorded using Neuropixels from area PM of a single representative mouse.

nat_movie = 1; % natural movie 1
area = 4; % area PM
mouse = 53; % mouse #53

% subset the neuronal data of the example mouse
% current_mouse - matrix of the size N by T by R (72 units x 30 time bins x 30 movie repeats)
current_mouse = neuropixels_population_vectors{mouse,area,nat_movie}(:,:,1:30);

% subset the first 10 movie repeats
subset_current_mouse = current_mouse(:,:,1:10);

% reshape 'subset_current_mouse' into 2D matrix over the columns (72 units x 300 time bins)
reshape_current_mouse = reshape(subset_current_mouse,[size(subset_current_mouse,1),...
    size(subset_current_mouse,2)*size(subset_current_mouse,3)]);

% calculate PV correlation across time bins
current_movie_pv_corr = corr(reshape_current_mouse);

figure % visualzing main panel
subplot(1,2,1,'units','normalized','position',[0.1 0.3 0.47 0.6])
imagesc(current_movie_pv_corr)
hold on
for line = 1:9
    plot(xlim,[30.5 30.5] +30*(line-1),'color',newmap3(1,:),'linewidth',1.25)
    plot([30.5 30.5] +30*(line-1),ylim,'color',newmap3(1,:),'linewidth',1.25)
end
set(gca,'xtick',15:30:300,'xticklabel',1:10,'ytick',15:30:300,'yticklabel',1:10)
xlabel('Movie repeat')
ylabel('Movie repeat')
title('Single animal example')
colormap(newmap3)
cb = colorbar;
set(cb,'position',[0.6 0.3 0.04 0.35])
cb.Label.String = 'PV correlation';
cb.FontSize = 12;

% Inset: average PV correlation over all pairs across different movie repeats.
current_mouse = neuropixels_population_vectors{mouse,area,nat_movie}(:,:,1:30);

% calculate the PV correlation between all time bins of one movie repeat
% and the time bins of a different movie repeat (results in a 30 by 30
% matrix for each pair), and stack tham into one 30 by 30 by 870 matrix
all_repeat_pv = [];
sub = 1; % will be used as an increasing index for the 3rd dim
for repeat1 = 1:30
    for repeat2 = 1:30
        if ~(repeat2 == repeat1) % perform only on different pairs of movie repeats
            all_repeat_pv(:,:,sub) = corr(current_mouse(:,:,repeat1),current_mouse(:,:,repeat2));
            sub = sub + 1;
        end
    end
end

average_across_pairs = mean(all_repeat_pv,3,'omitnan'); % average across pairs of movie repeats

% visualzing inset panel
subplot(1,2,2,'units','normalized','position',[0.6 0.7 0.15 0.2])
imagesc(average_across_pairs)
title('Time in movie')
set(gca,'xtick',[],'ytick',[])
colormap(newmap3)

%% Figure 2 B - Mean PV correlation for single PM example mouse
% Mean PV correlation for each pair of movie repeats from the same mouse shown in Fig. 2A.
% For visualization only, the diagonal was set to the maximal value.

nat_movie = 1; % natural movie 1
area = 4; % area PM
mouse = 53; % mouse #53

% current_mouse - subset binned neuronal activiy during first 30 movie repeats (block A)
% for the example mouse (size of 72 units x 30 time bins x 30 movie repeats)
current_mouse = neuropixels_population_vectors{mouse,area,nat_movie}(:,:,1:30);

% calculate the mean PV correlation between all pairs of movie repeats
mean_current_movie_pv_corr = []; % defining an empty variable that will contain the mean PV values

for repeat1 = 1:30 % loop over movie repeats
    for repeat2 = 1:30
        current_pv = corr(current_mouse(:,:,repeat1),current_mouse(:,:,repeat2)); % calculate PV corr between a time bins of a pair of movie repeats (30 by 30 in size)
        mean_current_movie_pv_corr(repeat1,repeat2) = mean(diag(current_pv),'omitnan'); % average across corresponding time bins (diagonal of the matrix)
    end
end

main_diag_ind = boolean(eye(size(mean_current_movie_pv_corr,1))); % indices of main diagonal
mean_current_movie_pv_corr(main_diag_ind) = NaN; % convert main diagonal into NaNs
mean_current_movie_pv_corr(main_diag_ind) = max(mean_current_movie_pv_corr(:)); % replace main diagonal with maximal value

figure % visualization
imagesc(mean_current_movie_pv_corr)
colormap(newmap3)
xlabel('Movie repeat')
ylabel('Movie repeat')
title('Single animal example')
cb = colorbar();
cb.Label.String = 'PV correlation';
cb.FontSize = 12;


%% Figure 2C - Mean PV correlation as fuction of interval for single PM example mouse
% Mean PV correlation as a function of time. Each data point represents the
% mean PV correlation value for a single pair of movie repeats from Fig. 2B.

nat_movie = 1; % natural movie 1
area = 4; % area PM
mouse = 53; % mouse #53

% current_mouse - subset binned neuronal activiy during first 30 movie repeats (block A)
% for the example mouse (size of 72 units x 30 time bins x 30 movie repeats)
current_mouse = neuropixels_population_vectors{mouse,area,nat_movie}(:,:,1:30);

% calculate the mean PV correlation between all pairs of movie repeats
mean_current_movie_pv_corr = []; % defining an empty variable that will contain the mean PV values

for repeat1 = 1:30 % loop over movie repeats
    for repeat2 = 1:30
        current_pv = corr(current_mouse(:,:,repeat1),current_mouse(:,:,repeat2)); % calculate PV corr between a time bins of a pair of movie repeats (30 by 30 in size)
        mean_current_movie_pv_corr(repeat1,repeat2) = mean(diag(current_pv),'omitnan'); % average across corresponding time bins (diagonal of the matrix)
    end
end

% calculating the mean PV correlation as function of elapsed time
mean_pv_elapse_repeat = []; % defing empty variable for mean values
elapse_repeat_values_mat = []; % defining empty variable for raw values
for diagonal = 1:size(mean_current_movie_pv_corr,1)-1 % loop over diagonals
    diag_values = diag(mean_current_movie_pv_corr,diagonal); % PV correlation values of current diagonal
    elapse_repeat_values_mat = [elapse_repeat_values_mat;[diag_values,ones(length(diag_values),1)*diagonal]]; % labeling raw values based on elpased time
    mean_pv_elapse_repeat(diagonal) = mean(diag_values,'omitnan'); % average all PV correlation values of each diagonal
end

figure % visualize PV correlation as function of elapsed time
hold on
scatter(elapse_repeat_values_mat(:,2),elapse_repeat_values_mat(:,1),25,[0.7 0.7 0.7],'filled')
plot(mean_pv_elapse_repeat,'color',newmap3(1,:),'linewidth',4)
text(0.05,0.075,'1 repeat = 30 seconds','Units','normalized','FontSize',11)
ylabel('PV correlation')
xlabel('Elapsed time (# of movie repeats)')
title('Single animal example')
set(gca,'xtick',[1 10 20 29],'ytick',[0.72:0.04:0.92])

%% Figure 2 D - Mean PV correlation matrices across mice and areas
% Mean PV correlation between movie repeats across animals and brain areas.

nat_movie = 1; % natural movie 1
cell_cutoff = 15; % cell count threshold of 15 units
num_repeats = 30; % number of movie repeats in each block

mean_pv_corr_all_areas = {}; % define an empty variable for average pv corr matrices across areas
clim = [];
% calculate for each visual area of each mouse the mean pv correlation between
% pairs of movie repeats (30 by 30 matrix).
for area = 1:6 % loop over brain areas
    % valid_mice - only mice from the 'Functional connectivity' group with at least 15 units
    valid_mice = neuropixels_cell_count(:,area,nat_movie) >= cell_cutoff & movie_repeats(:,nat_movie) == num_repeats;
    
    current_area = neuropixels_population_vectors(valid_mice,area,nat_movie); % subset all mice that met the requirements of 'valid mice'
    mean_pv_corr_all_mice = []; % define an empty variable for average pv corr matrices across mice (size of 30 by 30 by #mice)
    for mouse = 1:length(current_area) % loop over mice
        clc;
        disp(['Calculating PV correlation between movie repeats for Neuropixels data:'])
        disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        current_mouse = current_area{mouse}; % subset a single mouse
        current_mouse_blockA =  current_mouse(:,:,1:30); % subset neuronal activity during the first block (repeats 1-30)
        current_mouse_blockB =  current_mouse(:,:,31:60); % subset neuronal activity during the second block (repeats 31-60)
        
        current_mouse_mean_pv_corr = []; % define an empty variable for average pv corr across mocie repeats
        for repeat1 = 1:30 % loop over movie repeats
            for repeat2 = 1:30
                current_pv_blockA = corr(current_mouse_blockA(:,:,repeat1),current_mouse_blockA(:,:,repeat2)); % PV corr between time bins of pair of movie repeats in block A
                current_mouse_mean_pv_corr(repeat1,repeat2,1) = mean(diag(current_pv_blockA),'omitnan'); % average PV corr across corresponding time bin
                
                current_pv_blockB = corr(current_mouse_blockB(:,:,repeat1),current_mouse_blockB(:,:,repeat2)); % PV corr between time bins of pair of movie repeats in block B
                current_mouse_mean_pv_corr(repeat1,repeat2,2) = mean(diag(current_pv_blockB),'omitnan'); % average PV corr across corresponding time bin
            end
        end
        mean_pv_corr_all_mice(:,:,mouse) = mean(current_mouse_mean_pv_corr,3,'omitnan'); % average across blocks
        
    end
    mean_pv_corr_all_areas{area} =  mean_pv_corr_all_mice; % save matrices for current area
    
    % calculate and save color limits for visualization
    mean_pv_corr_all_mice = mean(mean_pv_corr_all_mice,3,'omitnan'); % average pv corr matrix across mice
    mean_pv_corr_all_mice(boolean(eye(size(mean_pv_corr_all_mice,1)))) = NaN; % set main diagonal to NaNs
    clim(area,:) = minmax(mean_pv_corr_all_mice(:)'); % extract min and max corr values (color limits)
end

% visualization of average PV correlation matrices across areas and mice
figure('units','normalized','position',[0.3 0.3 0.26 0.3])
for area = 1:6
    subplot(2,3,area)
    current_area = mean_pv_corr_all_areas{area}; % sebset pv corr matrices for the specific visual area
    average_current_area = mean(current_area,3,'omitnan'); % average pv corr matrices across mice
    main_diag_ind = boolean(eye(size(average_current_area,1))); % define main diagonal indeces
    average_current_area(main_diag_ind) = NaN; % set main diagonal to NaNs
    average_current_area(main_diag_ind) = max(clim(:)); % set main diagonal to maximal value
    imagesc(average_current_area,minmax(clim(:)'))
    if area == 5
        xlabel('Movie repeat')
    elseif area == 4 || area == 1
        ylabel('Movie repeat')
    elseif area == 6
        cb = colorbar;
        set(cb,'position',[0.925 0.1 0.04 0.75])
        
    end
    title([brain_areas{area},'            N=',num2str(size(mean_pv_corr_all_areas{area},3))])
end
suptitle('Average across animals')
colormap(jet)

%% Figure 2E - Mean PV correlation as function of elapsed time across mice

nat_movie = 1; % natural movie 1
cell_cutoff = 15; % cell count threshold of 15 units
num_repeats = 30; % number of movie repeats in each block

mean_pv_elapsed_repeat_per_area = {}; % define an empty variable for average pv corr across areas

% calculate for each visual area of each mouse the mean pv correlation between
% pairs of movie repeats (30 by 30 matrix).
for area = 1:6 % loop over brain areas
    % valid_mice - only mice from the 'Functional connectivity' group with at least 15 units
    valid_mice = neuropixels_cell_count(:,area,nat_movie) >= cell_cutoff & movie_repeats(:,nat_movie) == num_repeats;
    
    current_area = neuropixels_population_vectors(valid_mice,area,nat_movie); % subset all mice that met the requirements of 'valid mice'
    mean_pv_elapsed_repeat = []; % define an empty variable for average pv corr matrices across mice (size of #mice by 29)
    for mouse = 1:length(current_area) % loop over mice
        clc;
        disp(['Calculating PV correlation between movie repeats for Neuropixels data:'])
        disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        current_mouse = current_area{mouse}; % subset a single mouse
        current_mouse_blockA =  current_mouse(:,:,1:30); % subset neuronal activity during the first block (repeats 1-30)
        current_mouse_blockB =  current_mouse(:,:,31:60); % subset neuronal activity during the second block (repeats 31-60)
        
        current_mouse_mean_pv_corr = []; % define an empty variable for average pv corr across movie repeats
        for repeat1 = 1:30 % loop over movie repeats
            for repeat2 = 1:30 % loop over movie repeats
                current_pv_blockA = corr(current_mouse_blockA(:,:,repeat1),current_mouse_blockA(:,:,repeat2)); % PV corr between time bins of pair of movie repeats in block A
                current_mouse_mean_pv_corr(repeat1,repeat2,1) = mean(diag(current_pv_blockA),'omitnan'); % average PV corr across corresponding time bin
                
                current_pv_blockB = corr(current_mouse_blockB(:,:,repeat1),current_mouse_blockB(:,:,repeat2)); % PV corr between time bins of pair of movie repeats in block B
                current_mouse_mean_pv_corr(repeat1,repeat2,2) = mean(diag(current_pv_blockB),'omitnan'); % average PV corr across corresponding time bin
            end
        end
        mean_pv_corr_across_blocks = mean(current_mouse_mean_pv_corr,3,'omitnan'); % average across blocks
        
        % calculate pv correlation as function of elapsed time
        for diagonal = 1:num_repeats-1 % loop over diagonals
            mean_pv_elapsed_repeat(mouse,diagonal) = mean(diag(mean_pv_corr_across_blocks,diagonal),'omitnan'); % average pv corr values across diagonals
        end
    end
    mean_pv_elapsed_repeat_per_area{area} = mean_pv_elapsed_repeat; % store mean pv corr values across mice for each area
end

% friedman test parameters
pval = []; % define empty variable for the pvalues
df = []; % define empty variable for the degrees of freedom
chi = []; % define empty variable for the chi square values

figure % visualization of mean pv corr as function of time across areas and mice
for area = 1:6 % loop over areas
    current_area = mean_pv_elapsed_repeat_per_area{area}; % subset values of specific visual area
    [pval(area),tbl,stats] = friedman(current_area ,[1],'off'); % perform friedman test for main effect of elapsed time
    df(area) = tbl{2,3}; % store degrees of freedom
    chi(area) = tbl{2,5}; % store chi square values
    
    mean_across_mice = mean(current_area,'omitnan'); % calculate mean pv corr across mice
    std_across_mice = std(current_area,'omitnan'); % calculate standard deviation across mice
    ste_norm_vals = std_across_mice./sqrt(size(current_area,1)); % calculate standard error across mice
    
    hold on
    x = [1:length(mean_across_mice)]';
    y = mean_across_mice';
    dy = ste_across_mice';
    fill([x;flipud(x)],[y-dy;flipud(y+dy)],colors(area,:),'linestyle','none','facealpha',0.4);
    plt(area) = plot(y,'color',colors2(area,:),'linewidth',3);
end
lgd = legend(plt,brain_areas(1:6));
legend('boxoff')
lgd.Position = [0.15 0.225 0.15 0.15];
text(0.05,0.075,'1 repeat = 30 seconds','Units','normalized','FontSize',11)
ylabel('PV correlation')
xlabel('Elapsed time (# of movie repeats)')
title('Average across animals')
set(gca,'xtick',[1 10 20 29],'ytick',[0.66:0.04:0.84])

% correct for multiple comparisons using bonferroni-holm method
corrected_pval = bonf_holm(pval);

% define statistics summary table
VarNames = {'area','df','chi','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas(1:6)'),df(:),chi(:),pval(:),corrected_pval(:),'VariableNames',VarNames);

% display statistics summary table
clc;
disp(['PV correlation as a function of elapsed time'])
disp(['Friedman�s tests with Holm�Bonferroni correction:'])
disp(statistics)

%% Figure 2F - Activity rates of three PM units across movie repeats within a block

nat_movie = 1; % natural movie 1
mouse = 53; % mouse #53
area = 4; % area PM

% subset average neuronal activity of example mouse during block A (size of 72 units x 30 time bins:)
example_mouse = squeeze(mean(neuropixels_population_vectors{mouse,area,nat_movie}(:,:,1:30)*30,2,'omitnan'));

% visualize the neuronal activity of three example units and median across all units
figure('units','normalized','position',[0.3 0.3 0.26 0.38])
hold on
plot(nanmedian(example_mouse),'linewidth',8,'color',newmap3(1,:));
plot(example_mouse(32,:),'linewidth',4,'color',[0.8 0.8 0.8]);
plot(example_mouse(71,:),'linewidth',4,'color',[0.5 0.5 0.5]);
plot(example_mouse(26,:),'linewidth',4,'color',[0.2 0.2 0.2]);

text(0.25, 0.66,['Median activity'],'Units','normalized','color',newmap3(1,:),'Fontsize',12)
text(0.25, 0.61,['across all units'],'Units','normalized','color',newmap3(1,:),'Fontsize',12)
text(0.825, 0.775,['Unit #32'],'Units','normalized','color',[0.8 0.8 0.8],'Fontsize',12)
text(0.825, 0.575,['Unit #71'],'Units','normalized','color',[0.5 0.5 0.5],'Fontsize',12)
text(0.825, 0.15,['Unit #26'],'Units','normalized','color',[0.2 0.2 0.2],'Fontsize',12)

xlabel('Movie repeat')
ylabel('Activity rate (spike/sec)')
title('Single animal example:')
set(gca,'xtick',[1,10,20,30])

%% Figure 2G - responsiveness of three V1 example neurons across movie repeats spanning two blocks

nat_movie = 1; % natural movie 1
mouse = 53; % mouse #53
area = 1; % area V1

% subset neuronal activity of example mouse across two blocks
example_mouse = neuropixels_population_vectors{mouse,area,nat_movie};
units_list = [32,26,20]; % define a list of example units ind

for unit = 1:length(units_list) % loop over units
    current_unit = squeeze(example_mouse(units_list(unit),:,:)*30)'; % subset example unit
    
    % gaussian smoothing of neuronal activity in each movie repeat
    smooth_current_unit = []; % define empty variable for smoothen activity
    for repeat = 1:size(current_unit,1) % loop over movie repeats
        smooth_current_unit(repeat,:) = imgaussfilt(current_unit(repeat,:),2); % sigma = 2
    end
    % normalize neuronal activity based on peak firing rate in each movie repeat
    norm_current_unit = smooth_current_unit ./ max(smooth_current_unit,[],2);
    
    % visualize tuning curve of example unit across movie repeats and blocks
    subplot(2,3,unit)
    imagesc(norm_current_unit)
    hold on
    plot(xlim, [30 30]+0.5,'--','linewidth',2,'color','w')
    if unit ==1
        ylabel('Movie repeat')
        text(0.075, 0.9,['Block A'],'Units','normalized','color','w')
        text(0.075, 0.4,['Block B'],'Units','normalized','color','w')
    elseif unit == 3
        cb = colorbar;
        set(cb,'position',[0.925 0.575 0.04 0.35])
        
    end
    colormap(newmap3)
    title(['Unit #',num2str(units_list(unit))])
    
    % calculate example unit average tuning curve for each block
    mean_current_unit = [];
    mean_current_unit(1,:) = nanmean(current_unit(1:30,:));
    mean_current_unit(2,:) = nanmean(current_unit(31:60,:));
    
    %visualize example unit average tuning curve
    subplot(2,3,unit+3)
    hold on
    plot(mean_current_unit(1,:),'-','markersize',20,'linewidth',1.5,'color',[0.3 0.3 0.3])
    plot(mean_current_unit(2,:)','-','markersize',20,'linewidth',1.5,'color',[0.6 0.6 0.6])
    [r,p] = corr(mean_current_unit(1,:)',mean_current_unit(2,:)'); % tuning curve similarity between blocks
    text(0.1, 0.95,['r = ',num2str(r)],'Units','normalized','color',[0.25 0.25 0.25])
    xlim([1 30])
    if unit ==1
        ylabel('Mean activity rate (spike/sec)')
        legend({'Block A','Block B'})
        legend('boxoff')
        ylim([0 45])
    elseif unit ==2
        xlabel('Time in movie (sec)')
        ylim([0 60])
    else
        ylim([0 24])
    end
    
end

%% Figure 2H - Ensemble rate correlation as function of elapsed time across mice

nat_movie = 1; % natural movie 1
cell_cutoff = 15; % cell count threshold of 15 units
num_repeats = 30; % number of movie repeats in each block

mean_rate_elapsed_repeat_per_area = {};% define an empty variable for ensemble rate corr across areas

% calculate for each visual area of each mouse the ensemble rate correlation between
% pairs of movie repeats and as a function of elapsed time
for area = 1:6 % loop over brain areas
    % valid_mice - only mice from the 'Functional connectivity' group with at least 15 units
    valid_mice = neuropixels_cell_count(:,area,nat_movie) >= cell_cutoff & movie_repeats(:,nat_movie) == num_repeats;
    
    current_area = neuropixels_population_vectors(valid_mice,area,nat_movie); % subset all mice that met the requirements of 'valid mice'
    mean_rate_elapsed_repeat = [];  % define an empty variable for ensemble rate corr matrices across mice (size of #mice by 29)
    for mouse = 1:length(current_area) % loop over mice
        clc;
        disp(['Calculating ensemble rate correlation between movie repeats for Neuropixels data:'])
        disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        current_mouse = current_area{mouse}; % subset the neuronal activity of a single mouse (size of #units by 30 time bins by 60 movie repeats)
        current_mouse_blockA =  squeeze(mean(current_mouse(:,:,1:30),2,'omitnan')); % subset and calculate mean neuronal activity in each movie repeat during the first block (repeats 1-30)
        current_mouse_blockB =  squeeze(mean(current_mouse(:,:,31:60),2,'omitnan')); % subset and calculate mean neuronal activity in each movie repeat during the second block (repeats 31-60)
        
        % calculate ensemble rate correlation between pairs of movie
        % repeats in each block
        rate_corr = [];
        rate_corr(:,:,1) = corr(current_mouse_blockA); % ensemble rate correlation for block A
        rate_corr(:,:,2) = corr(current_mouse_blockB); % ensemble rate correlation for block B
        
        mean_rate_corr_across_blocks = mean(rate_corr,3,'omitnan'); % average ensemble rate correlation across blocks
        
        % calculate ensemble rate correlation as function of elapsed time
        for diagonal = 1:num_repeats-1 % loop over diagonals
            mean_rate_elapsed_repeat(mouse,diagonal) = mean(diag(mean_rate_corr_across_blocks,diagonal),'omitnan'); % average pv corr values across diagonals
        end
    end
    mean_rate_elapsed_repeat_per_area{area} = mean_rate_elapsed_repeat; % store ensemble rate corr values across mice for each area
end

% friedman test parameters
pval = []; % define empty variable for the pvalues
df = []; % define empty variable for the degrees of freedom
chi = []; % define empty variable for the chi square values


figure % visualization of ensemble rate corr as function of time across areas and mice
for area = 1:6 % loop over areas
    current_area = mean_rate_elapsed_repeat_per_area{area}; % subset values of specific visual area
    [pval(area),tbl,stats] = friedman(current_area ,[1],'off'); % perform friedman test for main effect of elapsed time
    df(area) = tbl{2,3}; % save degrees of freedom
    chi(area) = tbl{2,5}; % save chi square values
    
    mean_across_mice = mean(current_area,'omitnan'); % calculate mean ensemble rate corr across mice
    std_across_mice = std(current_area,'omitnan'); % calculate standard deviation across mice
    ste_across_mice = std_across_mice./sqrt(size(current_area,1)); % calculate standard error across mice
    
    hold on
    x = [1:length(mean_across_mice)]';
    y = mean_across_mice';
    dy = ste_across_mice';
    fill([x;flipud(x)],[y-dy;flipud(y+dy)],colors(area,:),'linestyle','none','facealpha',0.4);
    plt(area) = plot(y,'color',colors2(area,:),'linewidth',3);
end
lgd = legend(plt,brain_areas(1:6));
legend('boxoff')
lgd.Position = [0.15 0.225 0.15 0.15];
text(0.05,0.075,'1 repeat = 30 seconds','Units','normalized','FontSize',11)
ylabel('Ensemble rate correlation')
xlabel('Elapsed time (# of movie repeats)')
title('Average across animals')
set(gca,'xtick',[1 10 20 29],'ytick',[0.86:0.04:0.98])

% correct for multiple comparisons using bonferroni-holm method
corrected_pval = bonf_holm(pval);

% define statistics summary table
VarNames = {'area','df','chi','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas(1:6)'),df(:),chi(:),pval(:),corrected_pval(:),'VariableNames',VarNames);

% display statistics summary table
clc;
disp(['Ensemble rate correlation as a function of elapsed time'])
disp(['Friedman�s tests with Holm�Bonferroni correction:'])
disp(statistics)

%% Figure 2I - Tuning curve correlation as function of elapsed time across mice

nat_movie = 1; % natural movie 1
cell_cutoff = 15; % cell count threshold of 15 units
num_repeats = 30; % number of movie repeats in each block

mean_tuning_elapsed_repeat_per_area = {};% define an empty variable for tuning curve corr across areas

% calculate for each visual area of each mouse the tuning curve correlation between
% pairs of movie repeats as function of elapsed time
for area = 1:6 % loop over brain areas
    % valid_mice - only mice from the 'Functional connectivity' group with at least 15 units
    valid_mice = neuropixels_cell_count(:,area,nat_movie) >= cell_cutoff & movie_repeats(:,nat_movie) == num_repeats;
    
    current_area = neuropixels_population_vectors(valid_mice,area,nat_movie); % subset all mice that met the requirements of 'valid mice'
    mean_tuning_elapsed_repeat = [];% define an empty variable for tuning curve corr matrices across mice (size of #mice by 29)
    for mouse = 1:length(current_area) % loop over mice
        clc;
        disp(['Calculating tuning curve correlation between movie repeats for Neuropixels data:'])
        disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        current_mouse = current_area{mouse}; % subset a single mouse
        current_mouse_blockA =  current_mouse(:,:,1:30); % subset neuronal activity during the first block (repeats 1-30)
        current_mouse_blockB =  current_mouse(:,:,31:60); % subset neuronal activity during the second block (repeats 31-60)
        
        current_mouse_tuning_corr = []; % define an empty variable for average tuning curve across movie repeats
        for repeat1 = 1:30 % loop over movie repeats
            for repeat2 = 1:30 % loop over movie repeats
                current_tuning_blockA = corr(current_mouse_blockA(:,:,repeat1)',current_mouse_blockA(:,:,repeat2)'); % correlation between tuning curves of corresponding units across pair of movie repeats in block A
                current_mouse_tuning_corr(repeat1,repeat2,1) = median(diag(current_tuning_blockA),'omitnan'); % median tuning curve corr value across corresponding units
                imagesc(current_tuning_blockA)
                current_tuning_blockB = corr(current_mouse_blockB(:,:,repeat1)',current_mouse_blockB(:,:,repeat2)'); % correlation between tuning curves of corresponding units across pair of movie repeats in block B
                current_mouse_tuning_corr(repeat1,repeat2,2) = median(diag(current_tuning_blockB),'omitnan'); % median tuning curve corr value across corresponding units
            end
        end
        mean_tuning_corr_across_blocks = mean(current_mouse_tuning_corr,3,'omitnan'); % average across blocks
        
        % calculate tuning curve correlation as function of elapsed time
        for diagonal = 1:num_repeats-1 % loop over diagonals
            mean_tuning_elapsed_repeat(mouse,diagonal) = mean(diag(mean_tuning_corr_across_blocks,diagonal),'omitnan'); % average tuning curve corr values across diagonals
        end
    end
    mean_tuning_elapsed_repeat_per_area{area} = mean_tuning_elapsed_repeat; % store mean tuning curve values across mice for each area
end

% friedman test parameters
pval = []; % define empty variable for the pvalues
df = []; % define empty variable for the degrees of freedom
chi = []; % define empty variable for the chi square values

figure % visualization of mean tuning curve corr as function of time across areas and mice
for area = 1:6 % loop over areas
    current_area = mean_tuning_elapsed_repeat_per_area{area}; % subset values of specific visual area
    [pval(area),tbl,stats] = friedman(current_area,1,'off'); % perform friedman test for main effect of elapsed time
    df(area) = tbl{2,3};
    chi(area) = tbl{2,5};
    
    mean_across_mice = mean(current_area,'omitnan'); % calculate mean tuning curve corr across mice
    std_across_mice = std(current_area,'omitnan'); % calculate standard deviation across mice
    ste_across_mice = std_across_mice./sqrt(size(current_area,1)); % calculate standard error across mice
    
    hold on
    x = [1:length(mean_across_mice)]';
    y = mean_across_mice';
    dy = ste_across_mice';
    fill([x;flipud(x)],[y-dy;flipud(y+dy)],colors(area,:),'linestyle','none','facealpha',0.4);
    plt(area) = plot(y,'color',colors2(area,:),'linewidth',3);
end
lgd = legend(plt,brain_areas(1:6));
legend('boxoff')
lgd.Position = [0.15 0.225 0.15 0.15];
text(0.05,0.075,'1 repeat = 30 seconds','Units','normalized','FontSize',11)
ylabel('Tuning curve correlation')
xlabel('Elapsed time (# of movie repeats)')
title('Average across animals')
set(gca,'xtick',[1 10 20 29],'ytick',[0.1:0.1:0.5])

% correct for multiple comparisons using bonferroni-holm method
corrected_pval = bonf_holm(pval);

% define statistics summary table
VarNames = {'area','df','chi','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas(1:6)'),df(:),chi(:),pval(:),corrected_pval(:),'VariableNames',VarNames);

% display statistics summary table
clc;
disp(['Tuning curve correlation as a function of elapsed time'])
disp(['Friedman�s tests with Holm�Bonferroni correction:'])
disp(statistics)

%% Figure 3A - Between blocks stability - PV correlation for area LM across mice
% Calculating PV correlation between the 1st (repeats 1-2) and 2nd (repeats 3-5)
% halves of two different blocks of �Natural movie 3� in a single visual area (area LM).

nat_movie = 2; % natural movie 3
area = 2; % area LM
num_repeats = 5; % number of movie repeats in each block
cell_cutoff = 15; % cell count threshold of 15 units

% valid_mice - only mice from the 'Brain Observatory 1.1' group with at least 15 units
valid_mice = neuropixels_cell_count(:,area) >= cell_cutoff & movie_repeats(:,nat_movie) == num_repeats;
current_area = neuropixels_population_vectors(valid_mice,area,nat_movie); % subset all mice that met the requirements of 'valid mice'

pv_corr_across_mice = []; % define an empty variable for pvn corr matrices across mice (size of 120 by 120 by #mice)
for mouse = 1:length(current_area) % loop over mice
    current_mouse = current_area{mouse}; % subset a single mouse (size of #units by 30 time bins by 10 movie repeats)
    
    current_mouse_blockA =  current_mouse(:,:,1:5); % subset the neuronal activity during block A (repeats 1-5)
    current_mouse_blockA_half1 = mean(current_mouse_blockA(:,:,1:2),3,'omitnan'); % calculate the average activity during the first half of block A (repeats 1-2)
    current_mouse_blockA_half2 = mean(current_mouse_blockA(:,:,3:5),3,'omitnan'); % calculate the average activity during the second half of block A (repeats 3-5)
    
    current_mouse_blockB =  current_mouse(:,:,6:10); % subset the neuronal activity during block B (repeats 6-10)
    current_mouse_blockB_half1 = mean(current_mouse_blockB(:,:,1:2),3,'omitnan'); % calculate the average activity during the first half of block B (repeats 6-7)
    current_mouse_blockB_half2 = mean(current_mouse_blockB(:,:,3:5),3,'omitnan'); % calculate the average activity during the second half of block B (repeats 8-10)
    
    pv_corr_across_mice(:,:,mouse) = corr([current_mouse_blockA_half1,current_mouse_blockA_half2,...
        current_mouse_blockB_half1,current_mouse_blockB_half2]); % calculate the pv correlation across time bins of different halves and blocks
    
end

mean_pv_corr_across_mice = mean(pv_corr_across_mice,3,'omitnan'); % calculate average pv correlation across mice

figure % visualize average pv corr across mice for area LM
imagesc(mean_pv_corr_across_mice)
hold on
for line = 1:3
    plot(xlim,[30.5 30.5] +30*(line-1),'color',newmap3(1,:),'linewidth',1.25)
    plot([30.5 30.5] +30*(line-1),ylim,'color',newmap3(1,:),'linewidth',1.25)
end
set(gca,'xtick',15:30:120,'xticklabel',{'1st half','2nd half','1st half','2nd half'},...
    'ytick',15:30:120,'yticklabel',{'1st half','2nd half','1st half','2nd half'})
ytickangle(90)
xlabel('Block A                                           Block B')
ylabel('Block B                                  Block A')

title('Natural movie 3')
colormap(newmap3)
cb = colorbar;
cb.Label.String = 'PV correlation';
cb.FontSize = 12;

%% Figure 3B - Mean PV correlation within and between blocks of natural movie 3

nat_movie = 2; % natural movie 3
num_repeats = 5; % number of movie repeats in each block
cell_cutoff = 15; % cell count threshold of 15 units

pv_corr_between_blocks_across_areas = {}; % define empty variable that will store all pv corr values across blocks across mice and visual areas
for area = 1:6 % loop over visual areas
    % valid_mice - only mice from the 'Brain Observatory 1.1' group with at least 15 units
    valid_mice = neuropixels_cell_count(:,area) >= cell_cutoff & movie_repeats(:,nat_movie) == num_repeats;
    current_area = neuropixels_population_vectors(valid_mice,area,nat_movie); % subset all mice that met the requirements of 'valid mice'
    
    within_between_stability = []; % define empty variable that will store pv corr values for all mice
    for mouse = 1:length(current_area) % loop over mice
        clc;
        disp(['Calculating PV correlation between blocks for Neuropixels data:'])
        disp(['Stimulus: Natural movie 3 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse}; % subset a single mouse (size of #units by 30 time bins by 10 movie repeats)
        
        current_mouse_blockA =  current_mouse(:,:,1:5); % subset the neuronal activty for block A (repeats 1-5)
        current_mouse_blockB =  current_mouse(:,:,6:10); % subset the neuronal activty for block B (repeats 6-10)
        
        current_mouse_half_blocks = []; % define an empty variable that will store the mean neuronal activty for each half of the two blocks
        current_mouse_half_blocks(:,:,1) = mean(current_mouse_blockA(:,:,1:2),3,'omitnan'); % calculate the mean activity across repreats for the first half of block A
        current_mouse_half_blocks(:,:,2) = mean(current_mouse_blockA(:,:,3:5),3,'omitnan'); % calculate the mean activity across repreats for the second half of block A
        current_mouse_half_blocks(:,:,3) = mean(current_mouse_blockB(:,:,1:2),3,'omitnan'); % calculate the mean activity across repreats for the first half of block B
        current_mouse_half_blocks(:,:,4) = mean(current_mouse_blockB(:,:,3:5),3,'omitnan'); % calculate the mean activity across repreats for the second half of block B
        
        % calculate pv correlation between the four halves of the two blocks
        mean_pv_corr_across_mice = []; % define empty variable that will store the mean pv across halves for a given mouse
        for half1 = 1:4 % loop over blocks halves
            for half2 = 1:4 % loop over blocks halves
                current_pv_corr = corr(current_mouse_half_blocks(:,:,half1),current_mouse_half_blocks(:,:,half2)); % calculate the pv corr between time bins of different halves
                mean_pv_corr_across_mice(half1,half2) = mean(diag(current_pv_corr),'omitnan'); % calculate the mean pv corr across corresponding time bins
            end
        end
        
        % calculate and store the average pv corr between halves of the
        % same block ('within-block') and average pv corr between different
        % halves of different blocks ('between-blocks')
        within_block = [mean_pv_corr_across_mice(1,2),mean_pv_corr_across_mice(3,4)]; % subset pv corr between halves of the same block
        between_blocks = mean_pv_corr_across_mice(1:2,3:4); % subset pv corr between halves of different blocks
        
        within_between_stability(mouse,1) = mean(within_block(:),'omitnan'); % calculate the average pv corr within a block
        within_between_stability(mouse,2) = mean(between_blocks(:),'omitnan'); % calculate the average pv corr between blocks
        
    end
    pv_corr_between_blocks_across_areas{area} = within_between_stability; % store mean pv corr values of all mice of a given area
end

% define empty variables that will store the results of Wilcoxon signed-rank tests
pvalue = []; % define empty variable for the pvalues
zvalue = []; % define empty variable for the z values
figure('Units','Normalized','Position',[0.2 0.4 0.3 0.325]) % visualization of mean pv corr between blocks across areas and mice
for area = 1:6 % loop over areas
    current_area = pv_corr_between_blocks_across_areas{area}; % subset values of specific visual area
    
    mean_stability = mean(current_area,'omitnan'); % calculate mean pv corr across mice
    std_stability = std(current_area,'omitnan'); % calculate standard deviation across mice
    ste_stability = std_stability./sqrt(size(current_area,1)); % calculate standard error across mice
    
    [pvalue(area),~,stats] = signrank(current_area(:,1),current_area(:,2)); % perform two-sided wilcoxon signed-rank test for difference within-between blocks
    zvalue(area) = stats.zval; % store z values
    
    subplot(2,3,area)
    hold on
    errorbar(mean_stability,ste_stability,'o','color',[0.5 0.5 0.5],...
        'markerfacecolor',[0.5 0.5 0.5],'capsize',0,'linestyle','-','linewidth',3,'markersize',5)
    text(0.1,0.275,brain_areas{area},'Units','normalized','FontSize',12)
    text(0.05,0.15,['N=',num2str(size(current_area,1))],'Units','normalized','FontSize',10)
    
    xlim([0.5 2.5])
    ylim([0.75 1])
    
    if area >3
        set(gca,'xtick',[1,2],'xticklabel',{'Within block','Between blocks'})
        xtickangle(15)
    else
        set(gca,'xtick',[])
    end
    
    if area == 4 || area ==1
        ylabel('PV correlation')
    end
end

% correct for multiple comparisons using bonferroni-holm method
corrected_pvalue = bonf_holm(pvalue);

% define statistics summary table
VarNames = {'area','zvalue','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas(1:6)'),zvalue(:),pvalue(:),corrected_pvalue(:),'VariableNames',VarNames);

% display statistics summary table
clc;
disp(['PV correlation within blocks compared to between blocks'])
disp(['Two-tailed Wilcoxon signed-rank tests with Holm�Bonferroni correction:'])
disp(statistics)

%% Figure 3C - Ensemble rate correlation within and between blocks of natural movie 3

nat_movie = 2; % natural movie 3
num_repeats = 5; % number of movie repeats in each block
cell_cutoff = 15; % cell count threshold of 15 units

rate_corr_between_blocks_across_areas = {}; % define empty variable that will store all ensemble rate corr values across blocks across mice and visual areas
for area = 1:6 % loop over visual areas
    % valid_mice - only mice from the 'Brain Observatory 1.1' group with at least 15 units
    valid_mice = neuropixels_cell_count(:,area) >= cell_cutoff & movie_repeats(:,nat_movie) == num_repeats;
    current_area = neuropixels_population_vectors(valid_mice,area,nat_movie); % subset all mice that met the requirements of 'valid mice'
    
    within_between_stability = []; % define empty variable that will store pv corr values for all mice
    for mouse = 1:length(current_area) % loop over mice
        clc;
        disp(['Calculating ensemble rate correlation between blocks for Neuropixels data:'])
        disp(['Stimulus: Natural movie 3 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse}; % subset a single mouse (size of #units by 30 time bins by 10 movie repeats)
        
        current_mouse_blockA =  current_mouse(:,:,1:5); % subset the neuronal activty for block A (repeats 1-5)
        current_mouse_blockB =  current_mouse(:,:,6:10); % subset the neuronal activty for block B (repeats 6-10)
        
        % calculate ensemble rate vectors for each half of the two blocks
        current_mouse_half_blocks = []; % define an empty variable that will store the mean neuronal activty for each half of the two blocks
        current_mouse_half_blocks(:,1) = mean(mean(current_mouse_blockA(:,:,1:2),3,'omitnan'),2,'omitnan'); % calculate the mean activity across repreats and time bins for the first half of block A
        current_mouse_half_blocks(:,2) = mean(mean(current_mouse_blockA(:,:,3:5),3,'omitnan'),2,'omitnan'); % calculate the mean activity across repreats and time bins for the second half of block A
        current_mouse_half_blocks(:,3) = mean(mean(current_mouse_blockB(:,:,1:2),3,'omitnan'),2,'omitnan'); % calculate the mean activity across repreats and time bins for the first half of block B
        current_mouse_half_blocks(:,4) = mean(mean(current_mouse_blockB(:,:,3:5),3,'omitnan'),2,'omitnan'); % calculate the mean activity across repreats and time bins for the second half of block B
        
        rate_corr = corr(current_mouse_half_blocks); % calculate ensemble rate correlation between four halves of the two blocks
        
        % calculate and store the average pv corr between halves of the
        % same block ('within-block') and average pv corr between different
        % halves of different blocks ('between-blocks')
        within_block = [rate_corr(1,2),rate_corr(3,4)]; % subset ensemble rate corr between halves of the same block
        between_blocks = rate_corr(1:2,3:4); % subset ensemble rate corr between halves of different blocks
        
        within_between_stability(mouse,1) = mean(within_block(:),'omitnan'); % calculate the average ensemble rate corr within a block
        within_between_stability(mouse,2) = mean(between_blocks(:),'omitnan'); % calculate the average ensemble rate corr between blocks
        
    end
    rate_corr_between_blocks_across_areas{area} = within_between_stability; % store ensemble rate corr values of all mice of a given area
end

% define empty variables that will store the results of Wilcoxon signed-rank tests
pvalue = []; % define empty variable for the pvalues
zvalue = []; % define empty variable for the z values
plt = []; % store plot information for each area
figure('Units','Normalized','Position',[0.2 0.4 0.2 0.325]) % visualization of ensemble rate corr between blocks across areas and mice
for area = 1:6 % loop over areas
    current_area = rate_corr_between_blocks_across_areas{area}; % subset values of specific visual area
    
    mean_stability = mean(current_area,'omitnan'); % calculate mean ensemble rate corr across mice
    std_stability = std(current_area,'omitnan'); % calculate standard deviation across mice
    ste_stability = std_stability./sqrt(size(current_area,1)); % calculate standard error across mice
    
    [pvalue(area),~,stats] = signrank(current_area(:,1),current_area(:,2)); % perform two-sided wilcoxon signed-rank test for difference within-between blocks
    zvalue(area) = stats.zval; % store z values
    
    hold on
    plt(area) = errorbar(mean_stability,ste_stability,'o','color',colors(area,:),...
        'markerfacecolor',colors(area,:),'capsize',0,'linestyle','-','linewidth',3,'markersize',5);
    
end
xlim([0.5 2.5])
ylim([0.8255 1])
lgd = legend(plt,brain_areas(1:6));
legend('boxoff')
lgd.Position = [0.2 0.225 0.15 0.15];
set(gca,'xtick',[1,2],'xticklabel',{'Within block','Between blocks'},'ytick',0.85:0.05:1)
ylabel('Ensemble rate correlation')

% correct for multiple comparisons using bonferroni-holm method
corrected_pvalue = bonf_holm(pvalue);

% define statistics summary table
VarNames = {'area','zvalue','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas(1:6)'),zvalue(:),pvalue(:),corrected_pvalue(:),'VariableNames',VarNames);

% display statistics summary table
clc;
disp(['Ensemble rate correlation within blocks compared to between blocks'])
disp(['Two-tailed Wilcoxon signed-rank tests with Holm�Bonferroni correction:'])
disp(statistics)

%% Figure 3D - Tuning curve correlation within and between blocks of natural movie 3

nat_movie = 2; % natural movie 3
num_repeats = 5; % number of movie repeats in each block
cell_cutoff = 15; % cell count threshold of 15 units

tuning_corr_between_blocks_across_areas = {};  % define empty variable that will store all tuning curve corr values across blocks across mice and visual areas
for area = 1:6 % loop over visual areas
    % valid_mice - only mice from the 'Brain Observatory 1.1' group with at least 15 units
    valid_mice = neuropixels_cell_count(:,area) >= cell_cutoff & movie_repeats(:,nat_movie) == num_repeats;
    current_area = neuropixels_population_vectors(valid_mice,area,nat_movie); % subset all mice that met the requirements of 'valid mice'
    
    within_between_stability = []; % define empty variable that will store pv corr values for all mice
    for mouse = 1:length(current_area) % loop over mice
        clc;
        clc;
        disp(['Calculating tuning curve correlation between blocks for Neuropixels data:'])
        disp(['Stimulus: Natural movie 3 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse}; % subset a single mouse (size of #units by 30 time bins by 10 movie repeats)
        
        current_mouse_blockA =  current_mouse(:,:,1:5); % subset the neuronal activty for block A (repeats 1-5)
        current_mouse_blockB =  current_mouse(:,:,6:10); % subset the neuronal activty for block B (repeats 6-10)
        
        current_mouse_half_blocks = []; % define an empty variable that will store the mean neuronal activty for each half of the two blocks
        current_mouse_half_blocks(:,:,1) = mean(current_mouse_blockA(:,:,1:2),3,'omitnan'); % calculate the mean activity across repreats for the first half of block A
        current_mouse_half_blocks(:,:,2) = mean(current_mouse_blockA(:,:,3:5),3,'omitnan'); % calculate the mean activity across repreats for the second half of block A
        current_mouse_half_blocks(:,:,3) = mean(current_mouse_blockB(:,:,1:2),3,'omitnan'); % calculate the mean activity across repreats for the first half of block B
        current_mouse_half_blocks(:,:,4) = mean(current_mouse_blockB(:,:,3:5),3,'omitnan'); % calculate the mean activity across repreats for the second half of block B
        
        
        % calculate tuning curve correlation between the four halves of the two blocks
        tuning_corr_across_blocks = []; % define empty variable that will store the median tuning curve across halves for a given mouse
        for half1 = 1:4 % loop over blocks halves
            for half2 = 1:4 % loop over blocks halves
                current_tuning_corr = corr(current_mouse_half_blocks(:,:,half1)',current_mouse_half_blocks(:,:,half2)'); % calculate the tuning curve correlation between units across blocks halves
                tuning_corr_across_blocks(half1,half2) = median(diag(current_tuning_corr),'omitnan'); % calculate the median tuning curve corr across corresponding units
            end
        end
        
        % calculate and store the average tuning curve corr between halves of the
        % same block ('within-block') and average pv corr between different
        % halves of different blocks ('between-blocks')
        
        within_block = [tuning_corr_across_blocks(1,2),tuning_corr_across_blocks(3,4)]; % subset tuning curve corr between halves of the same block
        between_blocks = tuning_corr_across_blocks(1:2,3:4); % subset tuning curve corr between halves of different blocks
        
        within_between_stability(mouse,1) = mean(within_block(:),'omitnan'); % calculate the average tuning curve corr within a block
        within_between_stability(mouse,2) = mean(between_blocks(:),'omitnan'); % calculate the average tuning curve corr between blocks
        
    end
    tuning_corr_between_blocks_across_areas{area} = within_between_stability; % store mean tuning curve corr values of all mice of a given area
end

% define empty variables that will store the results of Wilcoxon signed-rank tests
pvalue = []; % define empty variable for the pvalues
zvalue = []; % define empty variable for the z values
plt = []; % store plot information for each area
figure('Units','Normalized','Position',[0.2 0.4 0.2 0.325]) % visualization of tuning curve corr between blocks across areas and mice
for area = 1:6 % loop over areas
    current_area = tuning_corr_between_blocks_across_areas{area}; % calculate mean tuning curve corr across mice
    
    mean_stability = mean(current_area,'omitnan'); % calculate mean tuning curve corr across mice
    std_stability = std(current_area,'omitnan'); % calculate standard deviation across mice
    ste_stability = std_stability./sqrt(size(current_area,1)); % calculate standard error across mice
    
    [pvalue(area),~,stats] = signrank(current_area(:,1),current_area(:,2)); % perform two-sided wilcoxon signed-rank test for difference within-between blocks
    zvalue(area) = stats.zval; % store z values
    
    hold on
    plt(area) = errorbar(mean_stability,ste_stability,'o','color',colors(area,:),...
        'markerfacecolor',colors(area,:),'capsize',0,'linestyle','-','linewidth',3,'markersize',5);
    
end
xlim([0.5 2.5])
ylim([0.35 0.7])
lgd = legend(plt,brain_areas(1:6));
legend('boxoff')
lgd.Position = [0.75 0.7 0.15 0.15];
set(gca,'xtick',[1,2],'xticklabel',{'Within block','Between blocks'},'ytick',0.4:0.1:0.7)
ylabel('Tuning curve correlation')

% correct for multiple comparisons using bonferroni-holm method
corrected_pvalue = bonf_holm(pvalue);

% define statistics summary table
VarNames = {'area','zvalue','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas(1:6)'),zvalue(:),pvalue(:),corrected_pvalue(:),'VariableNames',VarNames);

% display statistics summary table
clc;
disp(['Tuning curve correlation within blocks compared to between blocks'])
disp(['Two-tailed Wilcoxon signed-rank tests with Holm�Bonferroni correction:'])
disp(statistics)

%% Figure 3E - PV correlation across sessions for a single V1 mouse

area = 1; % area V1
nat_movie = 1; % natural movie 1
mouse = 3; % example mouse #3

% subset neuronal activity of example mouse (237 cells by 30 time bins by 30 movie repeats)
current_mouse = calcium_excitatory_population_vectors{area}{mouse,nat_movie};

% calculate average neuronal activity for each session
mean_current_mouse_sessions = []; % define empty variable that will store average neuronal activity
mean_current_mouse_sessions(:,:,1) =  mean(current_mouse(:,:,1:10),3,'omitnan'); % calculate mean activity rate across movie repeats for session 1 (repeats 1-10)
mean_current_mouse_sessions(:,:,2)  =  mean(current_mouse(:,:,11:20),3,'omitnan'); % calculate mean activity rate across movie repeats for session 2 (repeats 11-20)
mean_current_mouse_sessions(:,:,3)  =  mean(current_mouse(:,:,21:30),3,'omitnan'); % calculate mean activity rate across movie repeats for session 3 (repeats 21-30)

% calculate pv correlation across time bin of different sessions including
% only cells that were active in both compared sessions
pv_corr_pairwise_strict = []; % define empty variable that will store pv correlation across sessions
for sessionA = 1:3 % loop over sessions
    sessionA_activity = mean_current_mouse_sessions(:,:,sessionA); % subset neuronal activity for a single session
    rows = [1:30] +30*(sessionA-1); % define row indices for current session
    for sessionB = 1:3 % loop over sessions
        cols = [1:30] +30*(sessionB-1); % define column indices for current session
        sessionB_activity = mean_current_mouse_sessions(:,:,sessionB); % subset neuronal activity for a single session
        
        % valid_cells - cells that were active in both of compared sessions
        valid_cells = [mean(sessionA_activity,2,'omitnan') ~= 0] & [mean(sessionB_activity,2,'omitnan') ~= 0]; % true if average activity is above zero in both sessions
        valid_sessionA_activity = sessionA_activity(valid_cells,:); % subset only the cells that met the requirments of 'valid_cells'
        valid_sessionB_activity = sessionB_activity(valid_cells,:); % subset only the cells that met the requirments of 'valid_cells'
        
        pv_corr_pairwise_strict(rows,cols) = corr(valid_sessionA_activity,valid_sessionB_activity); % calculate the pv correlation between time bin of different sessions
        
    end
end

figure('unit','normalized','position',[0.3 0.3 0.215 0.3]) % visualize the pv correlation across sessions for example mouse
imagesc(pv_corr_pairwise_strict)
hold on
for line = 1:3
    plot(xlim,[30.5 30.5] *(line-1),'linewidth',2,'color',newmap3(1,:))
    plot([30.5 30.5] *(line-1),ylim,'linewidth',2,'color',newmap3(1,:))
end
colormap(newmap3)
set(gca,'xtick',15:30:90,'xticklabel',{'Session 1','Session 2','Session 3'},...
    'ytick',15:30:90,'yticklabel',{'Session 1','Session 2','Session 3'})
ytickangle(90)
title('Natural movie 1')
colormap(newmap3)
cb = colorbar;
cb.Label.String = 'PV correlation';
cb.FontSize = 12;

%% Figure 3F - Mean PV correlation within session and between sessions

nat_movie = 1; % natural movie 1
cell_cutoff = 20; % cell count threshold of 20 cells

within_between_session_stability_areas = {}; % define empty variable that will store the pv corr values of within and between sessions across mice and visual areas
for area = 1:6 % loop over visual areas
    % valid_mice - true only if the number of cells recorded is
    % atleast 20 in each of the recording sessions
    valid_mice = min(calcium_excitatory_cell_count{area},[],2)>= cell_cutoff;
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie); % mice that passed the requirments of 'valid_mice'
    
    within_between_session_stability = []; % define empty variable that will store the pv correlation values for all mice
    for mouse = 1:length(current_area) % loop over mice
        clc;
        disp(['Calculating PV correlation between sessions for Calcium imaging data:'])
        disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse}; % subset neuronal activity for a single mouse (#cells by 30 time bins by 30 movie repeats)
        
        %calculate the average activity rate for session halves across movie repeats
        mean_activity_per_half = []; % define an empty variable that will store the averaged neuronal responses of each session half
        for half = 1:6 % loop over session halves
            current_half = [1:5] + 5*(half-1); % define range of movie repeats for current half
            mean_activity_per_half(:,:,half) = mean(current_mouse(:,:,current_half),3,'omitnan'); % average the neuronal activity across chosen movie repeats
        end
        
        % calculate mean pv correlation between sessions halves using only
        % cells that were active in both compared time points (either
        % within or between sessions)
        mean_pv_corr = []; % define empty variable that will store pv corr values for all mice
        for halfA = 1:size(mean_activity_per_half,3) % loop over session halves
            halfA_activity = mean_activity_per_half(:,:,halfA); % subset neuronal activity for a single half
            
            for halfB = 1:size(mean_activity_per_half,3) % loop over session halves
                halfB_activity = mean_activity_per_half(:,:,halfB); % subset neuronal activity for a single half
                
                % valid_cells - cells that were active in both of compared time points
                valid_cells = [mean(halfA_activity,2,'omitnan') ~= 0] & [mean(halfB_activity,2,'omitnan') ~= 0]; % true if average activity is above zero in both session halves
                valid_halfA_activity = halfA_activity(valid_cells,:); % subset only the cells that met the requirments of 'valid_cells'
                valid_halfB_activity = halfB_activity(valid_cells,:); % subset only the cells that met the requirments of 'valid_cells'
                
                pv_corr = corr(valid_halfA_activity,valid_halfB_activity); % calculate the pv correlation between time bins of different halves
                mean_pv_corr(halfA,halfB) = mean(diag(pv_corr),'omitnan'); % calculate the mean pv corr across corresponding time bins
            end
        end
        
        within_session_values = [mean_pv_corr(1,2),mean_pv_corr(3,4),mean_pv_corr(5,6)]; % pv corr values between halves of the same session (within session)
        proximal_sessions_values = [mean_pv_corr(1:2,3:4),mean_pv_corr(3:4,5:6)]; % pv corr values between halves of proximal sessions (sessions 1&2 and sessions 2&3)
        distal_session_values = mean_pv_corr(1:2,5:6); % pv corr values between halves of distal sessions (sessions 1&3)
        
        within_between_session_stability(mouse,1) = mean(within_session_values(:),'omitnan'); % average pv corr values for within session
        within_between_session_stability(mouse,2) = mean(proximal_sessions_values(:),'omitnan'); % average pv corr values for proximal sessions
        within_between_session_stability(mouse,3) = mean(distal_session_values(:),'omitnan'); % average pv corr values for distal sessions
    end
    within_between_session_stability_areas{area} = within_between_session_stability; % store mean pv corr values for all mice of a given area
    
end

% define empty variables that will store the results of Wilcoxon signed-rank tests
pvalue = []; % define empty variable for the pvalues
zvalue = []; % define empty variable for the z values
ylims = [0.375,0.6; 0.425,0.65; 0.4,0.6; 0.35 0.55;0.2,0.35;0.3,0.55]; % y axis value limits for each visual area
figure('Units','Normalized','Position',[0.2 0.4 0.3 0.325]) % visualization of pv corr between sessions across areas and mice
for area = 1:6 % loop over visual areas
    current_area = within_between_session_stability_areas{area}; % subset pv corr values for specific area
    within_vs_between_sess = [current_area(:,1),mean(current_area(:,2:3),2,'omitnan')]; % calculate mean pv corr within and between sessions
    
    mean_stability = mean(within_vs_between_sess,'omitnan'); % calculate mean pv across mice
    std_stability = std(within_vs_between_sess,'omitnan'); % calculate standard deviation across mice
    ste_stability = std_stability./sqrt(size(within_vs_between_sess,1)); % calculate standard error across mice
    
    [pvalue(area),~,stats] = signrank(within_vs_between_sess(:,1),within_vs_between_sess(:,2)); % perform two-sided wilcoxon signed-rank test for difference within-between sessions
    zvalue(area) = stats.zval; % store z values
    
    subplot(2,3,area)
    hold on
    errorbar(mean_stability,ste_stability,'o','color',[0.5 0.5 0.5],...
        'markerfacecolor',[0.5 0.5 0.5],'capsize',0,'linestyle','-','linewidth',3,'markersize',5)
    text(0.1,0.275,brain_areas{area},'Units','normalized','FontSize',12)
    text(0.05,0.15,['N=',num2str(size(current_area,1))],'Units','normalized','FontSize',10)
    
    xlim([0.5 2.5])
    ylim([ylims(area,:)])
    
    if area >3
        set(gca,'xtick',[1,2],'xticklabel',{'Within session','Between sessions'})
        xtickangle(15)
    else
        set(gca,'xtick',[])
    end
    
    if area == 4 || area ==1
        ylabel('PV correlation')
    end
end

% correct for multiple comparisons using bonferroni-holm method
corrected_pvalue = bonf_holm(pvalue);

% define statistics summary table
VarNames = {'area','zvalue','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas(1:6)'),zvalue(:),pvalue(:),corrected_pvalue(:),'VariableNames',VarNames);

% display statistics summary table
clc;
disp(['PV correlation within sessions compared to between sessions'])
disp(['Two-tailed Wilcoxon signed-rank tests with Holm�Bonferroni correction:'])
disp(statistics)


%% Figure 3G - Ensemble rate correlation within session and between sessions
nat_movie = 1; % natural movie 1
cell_cutoff = 20; % cell count threshold of 20 cells

within_between_session_stability_areas = {}; % define empty variable that will store the ensemble rate corr values of within and between sessions across mice and visual areas
for area = 1:6 % loop over visual areas
    % valid_mice - true only if the number of cells recorded is
    % atleast 20 in each of the recording sessions
    valid_mice = min(calcium_excitatory_cell_count{area},[],2)>= cell_cutoff;
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie); % mice that passed the requirments of 'valid_mice'
    
    within_between_session_stability = []; % define empty variable that will store the ensemble rate correlation values for all mice
    for mouse = 1:length(current_area) % loop over mice
        clc;
        disp(['Calculating ensemble rate correlation between sessions for Calcium imaging data:'])
        disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse}; % subset neuronal activity for a single mouse (#cells by 30 time bins by 30 movie repeats)
        
        %calculate the average activity rate for session halves across movie repeats
        mean_activity_per_half = []; % define an empty variable that will store the averaged neuronal responses of each session half
        for half = 1:6 % loop over session halves
            current_half = [1:5] + 5*(half-1); % define range of movie repeats for current half
            mean_activity_per_half(:,:,half) = mean(current_mouse(:,:,current_half),3,'omitnan'); % average the neuronal activity across chosen movie repeats
        end
        
        % calculate mean pv correlation between sessions halves using only
        % cells that were active in both compared time points (either
        % within or between sessions)
        rate_corr = []; % define empty variable that will store ensemble rate corr values for all mice
        for halfA = 1:size(mean_activity_per_half,3) % loop over session halves
            halfA_activity = mean_activity_per_half(:,:,halfA); % subset neuronal activity for a single half
            for halfB = 1:size(mean_activity_per_half,3) % loop over session halves
                halfB_activity = mean_activity_per_half(:,:,halfB); % subset neuronal activity for a single half
                
                % valid_cells - cells that were active in both of compared time points
                valid_cells = [mean(halfA_activity,2,'omitnan') ~= 0] & [mean(halfB_activity,2,'omitnan') ~= 0]; % true if average activity is above zero in both session halves
                valid_halfA_activity = mean(halfA_activity(valid_cells,:),2,'omitnan'); % calculate average activity across time bins only for the cells that met the requirments of 'valid_cells'
                valid_halfB_activity = mean(halfB_activity(valid_cells,:),2,'omitnan'); %  calculate average activity across time bins only for cells that met the requirments of 'valid_cells'
                
                rate_corr(halfA,halfB) = corr(valid_halfA_activity,valid_halfB_activity); % calculate the ensemble rate correlation between time bin of different halves
            end
        end
        within_session_values = [rate_corr(1,2),rate_corr(3,4),rate_corr(5,6)]; % ensemble rate corr values between halves of the same session (within session)
        proximal_sessions_values = [rate_corr(1:2,3:4),rate_corr(3:4,5:6)]; % ensemble rate corr values between halves of proximal sessions (sessions 1&2 and sessions 2&3)
        distal_session_values = rate_corr(1:2,5:6); % ensemble rate corr values between halves of distal sessions (sessions 1&3)
        
        within_between_session_stability(mouse,1) = mean(within_session_values(:),'omitnan'); % average ensemble rate corr values for within session
        within_between_session_stability(mouse,2) = mean(proximal_sessions_values(:),'omitnan'); % average ensemble rate corr values for proximal sessions
        within_between_session_stability(mouse,3) = mean(distal_session_values(:),'omitnan'); % average ensemble rate corr values for distal sessions
        
    end
    within_between_session_stability_areas{area} = within_between_session_stability; % store mean ensemble rate corr values for all mice of a given area
    
end

% define empty variables that will store the results of Wilcoxon signed-rank tests
pvalue = []; % define empty variable for the pvalues
zvalue = []; % define empty variable for the z values
plt = []; % will store plot information for each visual area
figure('Units','Normalized','Position',[0.2 0.4 0.2 0.325]) % visualization of ensemble rate corr between sessions across areas and mice
for area = 1:6
    current_area = within_between_session_stability_areas{area}; % subset ensemble rate corr values for specific area
    within_vs_between_sess = [current_area(:,1),mean(current_area(:,2:3),2,'omitnan')]; % calculate mean ensemble rate corr within and between sessions
    
    mean_stability = mean(within_vs_between_sess,'omitnan'); % calculate mean pv across mice
    std_stability = std(within_vs_between_sess,'omitnan'); % calculate standard deviation across mice
    ste_stability = std_stability./sqrt(size(within_vs_between_sess,1)); % calculate standard error across mice
    
    [pvalue(area),~,stats] = signrank(within_vs_between_sess(:,1),within_vs_between_sess(:,2)); % perform two-sided wilcoxon signed-rank test for difference within-between sessions
    zvalue(area) = stats.zval; % store z values
    
    hold on
    plt(area) = errorbar(mean_stability,ste_stability,'o','color',colors(area,:),...
        'markerfacecolor',colors(area,:),'capsize',0,'linestyle','-','linewidth',3,'markersize',5);
    
end
xlim([0.5 2.5])
ylim([0.4 0.8])
lgd = legend(plt,brain_areas(1:6));
legend('boxoff')
lgd.Position = [0.2 0.225 0.15 0.15];
set(gca,'xtick',[1,2],'xticklabel',{'Within session','Between sessions'},'ytick',0.4:0.1:8)
ylabel('Ensemble rate correlation')

% correct for multiple comparisons using bonferroni-holm method
corrected_pvalue = bonf_holm(pvalue);

% define statistics summary table
VarNames = {'area','zvalue','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas(1:6)'),zvalue(:),pvalue(:),corrected_pvalue(:),'VariableNames',VarNames);

% display statistics summary table
clc;
disp(['Ensemble rate correlation within sessions compared to between sessions'])
disp(['Two-tailed Wilcoxon signed-rank tests with Holm�Bonferroni correction:'])
disp(statistics)


%% Figure 3H - Tuning curve correlation within session and between sessions
nat_movie = 1; % natural movie 1
cell_cutoff = 20; % cell count threshold of 20 cells

within_between_session_stability_areas = {}; % define empty variable that will store the tuning curve corr values of within and between sessions across mice and visual areas
for area = 1:6 % loop over visual areas
    % valid_mice - true only if the number of cells recorded is
    % above 20 in each of the recording sessions
    valid_mice = min(calcium_excitatory_cell_count{area},[],2)>= cell_cutoff;
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie); % mice that passed the requirments of 'valid_mice'
    
    within_between_session_stability = [];% define empty variable that will store the pv correlation values for all mice
    for mouse = 1:length(current_area) % loop over mice
        clc;
        disp(['Calculating tuning curve correlation between sessions for Calcium imaging data:'])
        disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse};% subset neuronal activity for a single mouse (#cells by 30 time bins by 30 movie repeats)
        %calculate the average activity rate for session halves across movie repeats
        mean_activity_per_half = []; % define an empty variable that will store the averaged neuronal responses of each session half
        for half = 1:6 % loop over session halves
            current_half = [1:5] + 5*(half-1); % define range of movie repeats for current half
            mean_activity_per_half(:,:,half) = mean(current_mouse(:,:,current_half),3,'omitnan'); % average the neuronal activity across chosen movie repeats
        end
        
        % calculate mean pv correlation between sessions halves using only
        % cells that were active in both compared time points (either
        % within or between sessions)
        median_tuning_corr = []; % define empty variable that will store tuning curve corr values for all mice
        for halfA = 1:size(mean_activity_per_half,3) % loop over session halves
            halfA_activity = mean_activity_per_half(:,:,halfA); % subset neuronal activity for a single half
            
            for halfB = 1:size(mean_activity_per_half,3) % loop over session halves
                halfB_activity = mean_activity_per_half(:,:,halfB); % subset neuronal activity for a single half
                
                % valid_cells - cells that were active in both of compared time points
                valid_cells = [mean(halfA_activity,2,'omitnan') ~= 0] & [mean(halfB_activity,2,'omitnan') ~= 0]; % true if average activity is above zero in both session halves
                valid_halfA_activity = halfA_activity(valid_cells,:); % subset only the cells that met the requirments of 'valid_cells'
                valid_halfB_activity = halfB_activity(valid_cells,:); % subset only the cells that met the requirments of 'valid_cells'
                
                tuning_corr = corr(valid_halfA_activity',valid_halfB_activity'); % calculate the tuning curve correlation between cell across different halves
                median_tuning_corr(halfA,halfB) = median(diag(tuning_corr),'omitnan'); % calculate the median tuning curve corr across corresponding cells
            end
        end
        
        within_session_values = [median_tuning_corr(1,2),median_tuning_corr(3,4),median_tuning_corr(5,6)]; % tuning curve corr values between halves of the same session (within session)
        proximal_sessions_values = [median_tuning_corr(1:2,3:4),median_tuning_corr(3:4,5:6)]; % tuning curve corr values between halves of proximal sessions (sessions 1&2 and sessions 2&3)
        distal_session_values = median_tuning_corr(1:2,5:6); % tuning curve corr values between halves of distal sessions (sessions 1&3)
        
        within_between_session_stability(mouse,1) = mean(within_session_values(:),'omitnan'); % average pv corr values for within session
        within_between_session_stability(mouse,2) = mean(proximal_sessions_values(:),'omitnan'); % average pv corr values for proximal sessions
        within_between_session_stability(mouse,3) = mean(distal_session_values(:),'omitnan'); % average pv corr values for distal sessions
    end
    within_between_session_stability_areas{area} = within_between_session_stability; % store tuning curve corr values for all mice of a given area
    
end

% define empty variables that will store the results of Wilcoxon signed-rank tests
pvalue = []; % define empty variable for the pvalues
zvalue = []; % define empty variable for the z values
plt = []; % will store plot information for each visual area
figure('Units','Normalized','Position',[0.2 0.4 0.2 0.325]) % visualization of tuning curve corr between sessions across areas and mice
for area = 1:6
    current_area = within_between_session_stability_areas{area}; % subset tuning curve corr values for specific area
    within_vs_between_sess = [current_area(:,1),mean(current_area(:,2:3),2,'omitnan')]; % calculate mean ensemble rate corr within and between sessions
    
    mean_stability = mean(within_vs_between_sess,'omitnan'); % calculate mean pv across mice
    std_stability = std(within_vs_between_sess,'omitnan'); % calculate standard deviation across mice
    ste_stability = std_stability./sqrt(size(within_vs_between_sess,1)); % calculate standard error across mice
    
    [pvalue(area),~,stats] = signrank(within_vs_between_sess(:,1),within_vs_between_sess(:,2)); % perform two-sided wilcoxon signed-rank test for difference within-between sessions
    zvalue(area) = stats.zval; % store z values
    
    hold on
    plt(area) = errorbar(mean_stability,ste_stability,'o','color',colors(area,:),...
        'markerfacecolor',colors(area,:),'capsize',0,'linestyle','-','linewidth',3,'markersize',5);
    
end
xlim([0.5 2.5])
ylim([0 0.425])
lgd = legend(plt,brain_areas(1:6));
legend('boxoff')
lgd.Position = [0.75 0.7 0.15 0.15];
set(gca,'xtick',[1,2],'xticklabel',{'Within session','Between sessions'},'ytick',0:0.1:4)
ylabel('Tuning curve correlation')

% correct for multiple comparisons using bonferroni-holm method
corrected_pvalue = bonf_holm(pvalue);

% define statistics summary table
VarNames = {'area','zvalue','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas(1:6)'),zvalue(:),pvalue(:),corrected_pvalue(:),'VariableNames',VarNames);

% display statistics summary table
clc;
disp(['Tuning curve correlation within sessions compared to between sessions'])
disp(['Two-tailed Wilcoxon signed-rank tests with Holm�Bonferroni correction:'])
disp(statistics)


%% Figure 3I - PV correlation between proximal and distal sessions

nat_movie = 1; % natural movie 1
cell_cutoff = 20; % cell count threshold of 20 cells

proximal_distal_session_stability_areas = {}; % define empty variable that will store the pv corr values between sessions across mice and visual areas
for area = 1:6 % loop over visual areas
    % valid_mice - true only if the number of cells recorded is
    % atleast 20 in each of the recording sessions
    valid_mice = min(calcium_excitatory_cell_count{area},[],2)>= cell_cutoff;
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie); % mice that passed the requirments of 'valid_mice'
    
    proximal_distal_session_stability = []; % define empty variable that will store the pv correlation values for all mice
    for mouse = 1:length(current_area) % loop over mice
        clc;
        disp(['Calculating PV correlation between proximal and distal sessions for Calcium imaging data:'])
        disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse}; % subset neuronal activity for a single mouse (#cells by 30 time bins by 30 movie repeats)
        
        %calculate the average activity rate across movie repeats for each session
        mean_activity_per_sess = []; % define an empty variable that will store the averaged neuronal responses of each session
        for sess = 1:3 % loop over session halves
            current_sess = [1:10] + 10*(sess-1); % define range of movie repeats for current session
            mean_activity_per_sess(:,:,sess) = mean(current_mouse(:,:,current_sess),3,'omitnan'); % average the neuronal activity across chosen movie repeats
        end
        
        % calculate mean pv correlation between sessions using only
        % cells that were active in both compared time points
        mean_pv_corr = []; % define empty variable that will store pv corr values for all mice
        for sessA = 1:size(mean_activity_per_sess,3) % loop over sessions
            sessA_activity = mean_activity_per_sess(:,:,sessA); % subset neuronal activity for a single session
            
            for sessB = 1:size(mean_activity_per_sess,3) % loop over sessions
                sessB_activity = mean_activity_per_sess(:,:,sessB); % subset neuronal activity for a single session
                
                % valid_cells - cells that were active in both of compared time points
                valid_cells = [mean(sessA_activity,2,'omitnan') ~= 0] & [mean(sessB_activity,2,'omitnan') ~= 0]; % true if average activity is above zero in both sessions
                valid_sessA_activity = sessA_activity(valid_cells,:); % subset only the cells that met the requirments of 'valid_cells'
                valid_sessB_activity = sessB_activity(valid_cells,:); % subset only the cells that met the requirments of 'valid_cells'
                
                pv_corr = corr(valid_sessA_activity,valid_sessB_activity); % calculate the pv correlation between time bins of different sessions
                mean_pv_corr(sessA,sessB) = nanmean(diag(pv_corr)); % calculate the mean pv corr across corresponding time bins
            end
        end
        
        proximal_sessions = mean(diag(mean_pv_corr,1),'omitnan'); % pv corr values between proximal sessions (sessions 1&2 and sessions 2&3)
        distal_sessions = diag(mean_pv_corr,2); % pv corr values between distal sessions (sessions 1&3)
        proximal_distal_session_stability(mouse,:) = [proximal_sessions,distal_sessions]; % store pv corr values of both proximal and distal values for all mice
    end
    proximal_distal_session_stability_areas{area} = proximal_distal_session_stability; % store mean pv corr values for all mice of a given area
end


% define empty variables that will store the results of Wilcoxon signed-rank tests
pvalue = []; % define empty variable for the pvalues
zvalue = []; % define empty variable for the z values
plt = []; % will store plot information for each visual area
mean_stability = []; % define empty variable that will store the mean pv corr values across mice for each area
ste_stability = []; % define empty variable that will store the standard error for the pv corr values across mice for each area
for area = 1:6
    current_area = proximal_distal_session_stability_areas{area}; % subset pv corr values for specific area
    proximal_vs_distal = [current_area(:,1)-current_area(:,2)]; % calculate the difference in pv corr values between proximal and distal sessions
    mean_stability(area) = mean(proximal_vs_distal,'omitnan'); % calculate the mean difference in pv corr values across mice
    std_stability = std(proximal_vs_distal,'omitnan'); % calculate the standard deviation across mice
    ste_stability(area) = std_stability./sqrt(length(proximal_vs_distal)); % calculate standard error across mice
    
    [pvalue(area),~,stats] = signrank(current_area(:,1),current_area(:,2),'tail','right'); % perform one-sided wilcoxon signed-rank test for difference between proximal and distal sessions
    zvalue(area) = stats.zval; % store z values
end



figure('Units','Normalized','Position',[0.2 0.4 0.185 0.325]) % visualization of pv corr between sessions across areas and mice
bar(mean_stability,'facecolor',[0.8 0.8 0.8],'edgecolor','none')
hold on
errorbar(mean_stability,ste_stability,'capsize',0,'linewidth',2,'color',[0.6 0.6 0.6],'linestyle','none','markerfacecolor',[0.6 0.6 0.6])
xlim([0 7])
set(gca,'xtick',[1:6],'xticklabel',brain_areas(1:6))
ylabel('PV correlation difference')
set(gca, 'box','off')

% correct for multiple comparisons using bonferroni-holm method
corrected_pvalue = bonf_holm(pvalue);

% define statistics summary table
VarNames = {'area','zvalue','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas(1:6)'),zvalue(:),pvalue(:),corrected_pvalue(:),'VariableNames',VarNames);

% display statistics summary table
clc;
disp(['PV correlation between proximal sessions compared to distal sessions'])
disp(['One-tailed Wilcoxon signed-rank tests with Holm�Bonferroni correction:'])
disp(statistics)


%% Figure 3J - Ensemble rate correlation as a function of elapsed days
nat_movie = 1; % natural movie 1
cell_cutoff = 20; % cell count threshold of 20 cells

rate_corr_across_sessions_areas = {}; % define empty variable that will store the ensemble rate corr values between sessions across mice and visual areas
elapsed_days_areas = {}; % define empty variable that will store the elapsed time (in days) between sessions across mice and visual areas
for area = 1:6 % loop over visual areas
    % valid_mice - true only if the number of cells recorded is
    % atleast 20 in each of the recording sessions
    valid_mice = min(calcium_excitatory_cell_count{area},[],2)>= cell_cutoff;
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie); % subset the neuronal activity of mice that passed the requirments of 'valid_mice'
    current_area_age = calcium_excitatory_sorted_mouse_age{area}(valid_mice,:); % subset the age in each session of mice that passed the requirments of 'valid_mice'
    
    rate_corr_across_sessions = []; % define empty variable that will store the ensemble rate correlation values for all mice
    elapsed_days = []; % define empty variable that will store the elapsed time (in days) between sessions for all mice
    for mouse = 1:length(current_area)% loop over mice
        clc;
        disp(['Calculating ensemble rate correlation between proximal and distal sessions for Calcium imaging data:'])
        disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse}; % subset neuronal activity of specific mouse
        current_mouse_area = current_area_age(mouse,:); % subset the age in each session for specific mouse
        
        %calculate the average activity rate across movie repeats for each session
        mean_activity_per_sess = []; % define an empty variable that will store the averaged neuronal responses of each session
        for sess = 1:3 % loop over session halves
            current_sess = [1:10] + 10*(sess-1); % define range of movie repeats for current session
            mean_activity_per_sess(:,:,sess) = mean(current_mouse(:,:,current_sess),3,'omitnan'); % average the neuronal activity across chosen movie repeats
        end
        
        % calculate the ensemble rate correlation between sessions using only
        % cells that were active in both compared time points
        rate_corr = []; % define empty variable that will store ensemble rate values for all mice
        for sessA = 1:size(mean_activity_per_sess,3) % loop over sessions
            sessA_activity = mean_activity_per_sess(:,:,sessA); % subset neuronal activity for a single session
            
            for sessB = 1:size(mean_activity_per_sess,3) % loop over sessions
                sessB_activity = mean_activity_per_sess(:,:,sessB); % subset neuronal activity for a single session
                
                % valid_cells - cells that were active in both of compared time points
                valid_cells = [mean(sessA_activity,2,'omitnan') ~= 0] & [mean(sessB_activity,2,'omitnan') ~= 0]; % true if average activity is above zero in both sessions
                valid_sessA_activity = mean(sessA_activity(valid_cells,:),2,'omitnan'); % calculate the average activity across time bins only for the cells that met the requirments of 'valid_cells'
                valid_sessB_activity = mean(sessB_activity(valid_cells,:),2,'omitnan'); % calculate the average activity across time bins only for the cells that met the requirments of 'valid_cells'
                
                rate_corr(sessA,sessB) = corr(valid_sessA_activity,valid_sessB_activity); % calculate the ensemble rate correlation between different sessions
            end
        end
        
        proximal_sessions = diag(rate_corr,1)'; % ensemble rate corr values between proximal sessions (sessions 1&2 and sessions 2&3)
        distal_sessions = diag(rate_corr,2); % ensemble rate corr values between distal sessions (sessions 1&3)
        
        rate_corr_across_sessions(mouse,:) = [proximal_sessions,distal_sessions]; % store ensemble rate corr values across sessions for all mice
        elapsed_days(mouse,:) = [current_mouse_area(2)-current_mouse_area(1),...
            current_mouse_area(3)-current_mouse_area(2),...
            current_mouse_area(3)-current_mouse_area(1)]; % calculate and store the elapsed time (in days) between sessions
    end
    rate_corr_across_sessions_areas{area} = rate_corr_across_sessions; % store ensemble rate corr values for all mice of a given area
    elapsed_days_areas{area} = elapsed_days; % store he elapsed time (in days) values for all mice of a given area
end

clearvars line % clear variable 'line' if it already exists
% define empty variables that will store the results of Pearson's correlation tests
r = []; % define empty variable for the Pearson's correlation coefficients
pvals =[]; % define empty variable for the pvalues
figure % visualize the ensemble rate corr as function of elapsed days
for area = 1:6 % loop over areas
    current_area = rate_corr_across_sessions_areas{area}; % subset ensemble rate corr values for specific area
    current_elapsed_day = elapsed_days_areas{area}; % subset elapsed days values for specific area
    
    % find and remove cased of elapsed day of 0 (two imaging sessions in the same
    % day (only 3 cases in area V1):
    within_day_intervals = sum(current_elapsed_day==0,2)>0; % find cases of elapsed day of zero
    current_area = current_area(~within_day_intervals,:); % remove mice with elapsed day zero
    current_elapsed_day = current_elapsed_day(~within_day_intervals,:); % remove mice with elapsed day zero
    
    % for each mouse, average the ensemble rate correlation values with the same interval
    % between sessions resulting in minimum of two or maximum of three ensemble rate corr values
    % per mouse (e.g. if a mouse performed three sessions on successive
    % days then the gap between session 1&2 and that of sessions 2&3 will
    % both be equal two 1 and the and the average value of the ensemble rate corr
    % of the two comperisons will be used):
    mean_correlations = []; % define empty variable to store mean ensemble rate corr values
    mean_elapsed_day = []; % define empty variable to store unique elapsed time values
    for mouse = 1:size(current_area,1) % loop over mice
        current_mouse = current_area(mouse,:); % subset ensemble rate corr values of specific mouse
        current_days = current_elapsed_day(mouse,:); % subset elapsed days values of specific mouse
        if current_days(1) == current_days(2) % test if the pairs of proximal sessions had the same gap
            % if true - average the ensemble rate corr and elapsed days values
            mean_correlations = [mean_correlations,[mean(current_mouse(1:2),'omitnan'),current_mouse(3)]]; % average ensemble rate corr values of proximal sessions
            mean_elapsed_day  = [mean_elapsed_day,[mean(current_days(1:2),'omitnan'),current_days(3)]]; % average elapsed day values of proximal sessions
        else % if false - append raw values without any manipulation
            mean_correlations = [mean_correlations,current_mouse];
            mean_elapsed_day  = [mean_elapsed_day,current_days];
        end
    end
    
    
    [~,i] = sort(mean_elapsed_day); % sort values based on elapsed time
    x = mean_elapsed_day(i)'; % sort values based on elapsed time
    y = mean_correlations(i)'; % sort values based on elapsed time
    
    % calculate 95% polynomial confidence interval for data
    [p,s] = polyfit(x,y,1);
    [yfit,dy] = polyconf(p,x,s,'predopt','curve');
    
    [r(area),pvals(area)] = corr(mean_elapsed_day',mean_correlations','tail','left'); % perform one-sided Pearson's correlation test
    
    subplot(2,3,area)
    hold on
    scatter(mean_elapsed_day',mean_correlations',20,[0.7 0.7 0.7],'filled','markerfacealpha',0.4)
    h = fill([x;flipud(x)],[yfit-dy;flipud(yfit+dy)],[0.6 0.6 0.6],'linestyle','none');
    line(x,yfit,'color',[0.6 0.6 0.6],'linewidth',2)
    set(h,'facealpha',.25)
    
    max_interval = max(mean_elapsed_day); % find max interval to define x limits for visualization
    ylim([-0.2 1.2])
    xlim([0 max_interval+1])
    
    
    if area == 4
        ylabel('Ensemble rate correlation')
    elseif area == 5
        xlabel('Days between recording sessions')
    end
    
end

% correct for multiple comparisons using bonferroni-holm method
corrected_pvalue = bonf_holm(pvals);

% add statistics to the plot
for area = 1:6
    hold on
    subplot(2,3,area)
    text(0.1, 1,[brain_areas{area}],'Units','normalized','color',[0 0 0],'FontSize',12)
    text(0.6, 1,['r = ',num2str(round(r(area),2))],'Units','normalized','color',[0 0 0])
    text(0.6, 0.925,['p = ',num2str(corrected_pvalue(area))],'Units','normalized','color',[0 0 0])
end

%% Figure 3K - Tuning curve correlation as a function of elapsed days
nat_movie = 1; % natural movie 1
cell_cutoff = 20; % cell count threshold of 20 cells

tuning_corr_across_sessions_areas = {}; % define empty variable that will store the tuning curve corr values between sessions across mice and visual areas
elapsed_days_areas = {}; % define empty variable that will store the elapsed time (in days) between sessions across mice and visual areas
for area = 1:6 % loop over visual areas
    % valid_mice - true only if the number of cells recorded is
    % atleast 20 in each of the recording sessions
    valid_mice = min(calcium_excitatory_cell_count{area},[],2)>= cell_cutoff;
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie); % subset the neuronal activity of mice that passed the requirments of 'valid_mice'
    current_area_age = calcium_excitatory_sorted_mouse_age{area}(valid_mice,:); % subset the age in each session of mice that passed the requirments of 'valid_mice'
    
    tuning_corr_across_sessions = []; % define empty variable that will store the tuning curve correlation values for all mice
    elapsed_days = []; % define empty variable that will store the elapsed time (in days) between sessions for all mice
    for mouse = 1:length(current_area) % loop over mice
        clc;
        disp(['Calculating tuning curve correlation between proximal and distal sessions for Calcium imaging data:'])
        disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse}; % subset neuronal activity of specific mouse
        current_mouse_area = current_area_age(mouse,:); % subset the age in each session for specific mouse
        
        %calculate the average activity rate across movie repeats for each session
        mean_activity_per_sess = []; % define an empty variable that will store the averaged neuronal responses of each session
        for sess = 1:3 % loop over session halves
            current_sess = [1:10] + 10*(sess-1); % define range of movie repeats for current session
            mean_activity_per_sess(:,:,sess) = mean(current_mouse(:,:,current_sess),3,'omitnan'); % average the neuronal activity across chosen movie repeats
        end
        
        % calculate the tuning curve correlation between sessions using only
        % cells that were active in both compared time points
        median_tuning_corr = []; % define empty variable that will store tuning curve values for all mice
        for sessA = 1:size(mean_activity_per_sess,3) % loop over sessions
            sessA_activity = mean_activity_per_sess(:,:,sessA); % subset neuronal activity for a single session
            
            for sessB =1:size(mean_activity_per_sess,3) % loop over sessions
                sessB_activity = mean_activity_per_sess(:,:,sessB); % subset neuronal activity for a single session
                
                % valid_cells - cells that were active in both of compared time points
                valid_cells = [mean(sessA_activity,2,'omitnan') ~= 0] & [mean(sessB_activity,2,'omitnan') ~= 0]; % true if average activity is above zero in both sessions
                valid_sessA_activity = sessA_activity(valid_cells,:)'; % calculate the average activity across time bins only for the cells that met the requirments of 'valid_cells'
                valid_sessB_activity = sessB_activity(valid_cells,:)'; % calculate the average activity across time bins only for the cells that met the requirments of 'valid_cells'
                
                tuning_corr = corr(valid_sessA_activity,valid_sessB_activity); % calculate the tuning curve correlation between sessions
                median_tuning_corr(sessA,sessB) = median(diag(tuning_corr),'omitnan'); % calculate the median tuning curve across corresponding cells
                
            end
        end
        
        proximal_sessions = diag(median_tuning_corr,1)'; % tuning curve corr values between proximal sessions (sessions 1&2 and sessions 2&3)
        distal_sessions = diag(median_tuning_corr,2); % tuning curve corr values between distal sessions (sessions 1&3)
        
        tuning_corr_across_sessions(mouse,:) = [proximal_sessions,distal_sessions]; % store tuning curve corr values across sessions for all mice
        elapsed_days(mouse,:) = [current_mouse_area(2)-current_mouse_area(1),...
            current_mouse_area(3)-current_mouse_area(2),...
            current_mouse_area(3)-current_mouse_area(1)]; % calculate and store the elapsed time (in days) between sessions
    end
    tuning_corr_across_sessions_areas{area} = tuning_corr_across_sessions;
    elapsed_days_areas{area} = elapsed_days;
end


clearvars line % clear variable 'line' if it already exists
% define empty variables that will store the results of Pearson's correlation tests
r = []; % define empty variable for the Pearson's correlation coefficients
pvals =[]; % define empty variable for the pvalues
figure % visualize the ensemble rate corr as function of elapsed days
for area = 1:6 % loop over areas
    current_area = tuning_corr_across_sessions_areas{area}; % subset tuning curve corr values for specific area
    current_elapsed_day = elapsed_days_areas{area}; % subset elapsed days values for specific area
    
    % find and remove cased of elapsed day of 0 (two imaging sessions in the same
    % day (only 3 cases in area V1):
    within_day_intervals = sum(current_elapsed_day==0,2)>0; % find cases of elapsed day of zero
    current_area = current_area(~within_day_intervals,:); % remove mice with elapsed day zero
    current_elapsed_day = current_elapsed_day(~within_day_intervals,:); % remove mice with elapsed day zero
    
    % for each mouse, average the tuning curve correlation values with the same interval
    % between sessions resulting in minimum of two or maximum of three tuning curve corr values
    % per mouse (e.g. if a mouse performed three sessions on successive
    % days then the gap between session 1&2 and that of sessions 2&3 will
    % both be equal two 1 and the and the average value of the tuning curve corr
    % of the two comperisons will be used):
    mean_correlations = []; % define empty variable to store mean tuning curve corr values
    mean_elapsed_day = []; % define empty variable to store unique elapsed time values
    for mouse = 1:size(current_area,1) % loop over mice
        current_mouse = current_area(mouse,:); % subset tuning curve corr values of specific mouse
        current_days = current_elapsed_day(mouse,:); % subset elapsed days values of specific mouse
        if current_days(1) == current_days(2) % test if the pairs of proximal sessions had the same gap
            % if true - average the tuning curve corr and elapsed days values
            mean_correlations = [mean_correlations,[mean(current_mouse(1:2),'omitnan'),current_mouse(3)]]; % average ensemble rate corr values of proximal sessions
            mean_elapsed_day  = [mean_elapsed_day,[mean(current_days(1:2),'omitnan'),current_days(3)]]; % average elapsed day values of proximal sessions
        else % if false - append raw values without any manipulation
            mean_correlations = [mean_correlations,current_mouse];
            mean_elapsed_day  = [mean_elapsed_day,current_days];
        end
    end
    
    
    [~,i] = sort(mean_elapsed_day); % sort values based on elapsed time
    x = mean_elapsed_day(i)'; % sort values based on elapsed time
    y = mean_correlations(i)'; % sort values based on elapsed time
    
    % calculate 95% polynomial confidence interval for data
    [p,s] = polyfit(x,y,1);
    [yfit,dy] = polyconf(p,x,s,'predopt','curve');
    
    [r(area),pvals(area)] = corr(mean_elapsed_day',mean_correlations','tail','left'); % perform one-sided Pearson's correlation test
    
    subplot(2,3,area)
    hold on
    scatter(mean_elapsed_day',mean_correlations',20,[0.7 0.7 0.7],'filled','markerfacealpha',0.4)
    h = fill([x;flipud(x)],[yfit-dy;flipud(yfit+dy)],[0.6 0.6 0.6],'linestyle','none');
    line(x,yfit,'color',[0.6 0.6 0.6],'linewidth',2)
    set(h,'facealpha',.25)
    
    max_interval = max(mean_elapsed_day); % find max interval to define x limits for visualization
    ylim([-0.2 1.2])
    xlim([0 max_interval+1])
    
    if area == 4
        ylabel('Tuning curve correlation')
    elseif area == 5
        xlabel('Days between recording sessions')
    end
    
end

% correct for multiple comparisons using bonferroni-holm method
corrected_pvalue = bonf_holm(pvals);

% add statistics to the plot
for area = 1:6
    hold on
    subplot(2,3,area)
    text(0.1, 1,[brain_areas{area}],'Units','normalized','color',[0 0 0],'FontSize',12)
    text(0.6, 1,['r = ',num2str(round(r(area),2))],'Units','normalized','color',[0 0 0])
    text(0.6, 0.925,['p = ',num2str(round(corrected_pvalue(area),4))],'Units','normalized','color',[0 0 0])
end

[~,i] = sort(pvals','descend')
pvals(i)

%% Figure 4A - PV correlation across days grouped by layers and areas
% PV correlation between the two halves of the same session (within session),
% between halves of two temporally proximal sessions (proximal session)
% and between halves of two temporally distal sessions (distal session)
% grouped based on cortical layers for each of six visual areas

depth_list = [150 250; 251 350; 351 500; 501 700]; % define imaging depth ranges for each layer
layer_list = {'L2/3','L4','L5','L6'}; %  define layer names

nat_movie = 1; % natural movie 1
cell_cutoff = 20; % cell count threshold of 20 cells

valid_mice_num = []; % define empty variable that will store the number of FOV in each layer (for visualization)
elapsed_sess_pv_depth = {}; % define empty variable that will store the pv corr values of within and between sessions across mice and visual areas
for area = 1:6 % loop over areas
    
    for depth = 1:size(depth_list,1) % loop over layers
        
        % valid_mice_num - true only if the number of cells recorded is
        % atleast 20 in each of the recording sessions
        valid_cell_num = min(calcium_excitatory_cell_count{area},[],2) >= cell_cutoff;
        
        % valid_depth - true only if the imaging depth of the FOV is
        % within the specific imaging depth range
        valid_depth = calcium_excitatory_imaging_depth{area} >= depth_list(depth,1) & calcium_excitatory_imaging_depth{area} <= depth_list(depth,2);
        
        % valid_mice - mice that passed the reqirements of both
        % 'valid_mice_num' and 'valid_depth'
        valid_mice = valid_cell_num & valid_depth';
        valid_mice_num(area,depth) = sum(valid_mice); % store the number of FOV that passed the requirments of 'valid mice'
        
        elapsed_sess_pv = []; % define empty variable that will store the pv corr values across mice
        if valid_mice_num(area) >= 0 % test if there are any FOV that passed the criteria of 'valid_mice'
            
            current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie); % subsetthe neuronal activity of mice that paseed the the criteria of 'valid_mice'
            
            for mouse = 1:size(current_area,1) % loop over mice
                clc;
                disp(['Calculating PV correlation across layers for Calcium imaging data:'])
                disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Depth: ',layer_list{depth},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
                
                current_mouse = current_area{mouse}; % subset neuronal activity for a single mouse (#cells by 30 time bins by 30 movie repeats)
                
                %calculate the average activity rate for session halves across movie repeats
                mean_activity_per_half = []; % define an empty variable that will store the averaged neuronal responses of each session half
                for half = 1:6 % loop over session halves
                    current_half = [1:5] + 5*(half-1); % define range of movie repeats for current half
                    mean_activity_per_half(:,:,half) = mean(current_mouse(:,:,current_half),3,'omitnan'); % average the neuronal activity across chosen movie repeats
                end
                
                % calculate mean pv correlation between sessions halves using only
                % cells that were active in both compared time points (either
                % within or between sessions)
                mean_pv_corr = []; % define empty variable that will store pv corr values for all mice
                for halfA = 1:size(mean_activity_per_half,3) % loop over session halves
                    halfA_activity = mean_activity_per_half(:,:,halfA); % subset neuronal activity for a single half
                    
                    for halfB = 1:size(mean_activity_per_half,3) % loop over session halves
                        halfB_activity = mean_activity_per_half(:,:,halfB); % subset neuronal activity for a single half
                        
                        % valid_cells - cells that were active in both of compared time points
                        valid_cells = [mean(halfA_activity,2,'omitnan') ~= 0] & [mean(halfB_activity,2,'omitnan') ~= 0]; % true if average activity is above zero in both session halves
                        valid_halfA_activity = halfA_activity(valid_cells,:); % subset only the cells that met the requirments of 'valid_cells'
                        valid_halfB_activity = halfB_activity(valid_cells,:); % subset only the cells that met the requirments of 'valid_cells'
                        
                        pv_corr = corr(valid_halfA_activity,valid_halfB_activity); % calculate the pv correlation between time bins of different halves
                        mean_pv_corr(halfA,halfB) = mean(diag(pv_corr),'omitnan'); % calculate the mean pv corr across corresponding time bins
                    end
                end
                
                within_session_values = [mean_pv_corr(1,2),mean_pv_corr(3,4),mean_pv_corr(5,6)]; % pv corr values between halves of the same session (within session)
                proximal_sessions_values = [mean_pv_corr(1:2,3:4),mean_pv_corr(3:4,5:6)]; % pv corr values between halves of proximal sessions (sessions 1&2 and sessions 2&3)
                distal_session_values = mean_pv_corr(1:2,5:6); % pv corr values between halves of distal sessions (sessions 1&3)
                
                elapsed_sess_pv(mouse,1)=  mean(within_session_values(:),'omitnan'); % average pv corr values for within session;
                elapsed_sess_pv(mouse,2)= mean(proximal_sessions_values(:),'omitnan'); % average pv corr values for proximal sessions
                elapsed_sess_pv(mouse,3)= mean(distal_session_values(:),'omitnan'); % average pv corr values for distal sessions
                
                
            end
            
        end
        elapsed_sess_pv_depth{area,depth} = elapsed_sess_pv; % store mean pv corr values for all mice of a given area
    end
end


ylims = [0 0.7; 0.3 0.7; 0.3 0.7; 0.175 0.7; 0 0.5; 0.2 0.7]; % y axis value limits for each visual area
figure('units','normalized','position',[0.3 0.3 0.4 0.4]) % visualization of pv corr between sessions across areas and imaging depth
for area = 1:6 % loop over visual areas
    for depth = 1:size(depth_list,1) % loop over layers
        curent_area_depth = elapsed_sess_pv_depth{area,depth};  % subset pv corr values for specific area and cortical layer
        
        if ~isempty(curent_area_depth) % test if there are mice recorded from specific combination of cortical layer and visual area
            
            mean_stability = mean(curent_area_depth,'omitnan'); % calculate mean pv across mice
            std_stability = std(curent_area_depth,'omitnan'); % calculate standard deviation across mice
            ste_stability=std_stability./sqrt(size(curent_area_depth,1)); % calculate standard error across mice
            
            subplot(2,3,area)
            hold on
            errorbar(mean_stability ,ste_stability ,'o','color',new_jet_colormap(depth+20*(depth-1),:),...
                'markerfacecolor',new_jet_colormap(depth+20*(depth-1),:),'capsize',0,'linestyle','-','linewidth',3,'markersize',3)
            
            if area == 4 || area == 1
                ylabel('PV correlation')
            end
            
            if area > 3
                set(gca,'xtick',[1:3],'xticklabel',{'Within session','Proximal sessions','Distal sessions'})
                xtickangle(15)
            else
                set(gca,'xtick',[])
            end
            xlim([0.5 3.5])
        end
    end
    if area == 1
        lgd = legend({'L2/3','L4','L5','L6'});
        legend('boxoff');
        lgd.Position = [0.18 0.65 0.025 0.025];
    end
    ylim(ylims(area,:))
    text(0.7,0.9,brain_areas{area},'Units','normalized','FontSize',15)
    text(0.67,0.8,['N=',num2str(sum(valid_mice_num(area,:),2))],'Units','normalized','FontSize',12)
    
end

%% Figure 4B - PV correlation across days grouped by layers
% PV correlation between the two halves of the same session (within session),
% between halves of two temporally proximal sessions (proximal session)
% and between halves of two temporally distal sessions (distal session)
% grouped based on cortical layers

depth_list = [150 250; 251 350; 351 500; 501 700]; % define imaging depth ranges for each layer
layer_list = {'L2/3','L4','L5','L6'}; %  define layer names

nat_movie = 1; % natural movie 1
cell_cutoff = 20; % cell count threshold of 20 cells

nat_movie = 1; % natural movie 1
cell_cutoff = 20; % cell count threshold of 20 cells

valid_mice_num = []; % define empty variable that will store the number of FOV in each layer (for visualization)
elapsed_sess_pv_depth = {}; % define empty variable that will store the pv corr values of within and between sessions across mice and visual areas
for area = 1:6 % loop over areas
    
    for depth = 1:size(depth_list,1) % loop over layers
        
        % valid_mice_num - true only if the number of cells recorded is
        % atleast 20 in each of the recording sessions
        valid_cell_num = min(calcium_excitatory_cell_count{area},[],2) >= cell_cutoff;
        
        % valid_depth - true only if the imaging depth of the FOV is
        % within the specific imaging depth range
        valid_depth = calcium_excitatory_imaging_depth{area} >= depth_list(depth,1) & calcium_excitatory_imaging_depth{area} <= depth_list(depth,2);
        
        % valid_mice - mice that passed the reqirements of both
        % 'valid_mice_num' and 'valid_depth'
        valid_mice = valid_cell_num & valid_depth';
        valid_mice_num(area,depth) = sum(valid_mice); % store the number of FOV that passed the requirments of 'valid mice'
        
        elapsed_sess_pv = []; % define empty variable that will store the pv corr values across mice
        if valid_mice_num(area) >= 0 % test if there are any FOV that passed the criteria of 'valid_mice'
            
            current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie); % subsetthe neuronal activity of mice that paseed the the criteria of 'valid_mice'
            
            for mouse = 1:size(current_area,1) % loop over mice
                clc;
                disp(['Calculating PV correlation across layers for Calcium imaging data:'])
                disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Depth: ',layer_list{depth},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
                
                current_mouse = current_area{mouse}; % subset neuronal activity for a single mouse (#cells by 30 time bins by 30 movie repeats)
                
                %calculate the average activity rate for session halves across movie repeats
                mean_activity_per_half = []; % define an empty variable that will store the averaged neuronal responses of each session half
                for half = 1:6 % loop over session halves
                    current_half = [1:5] + 5*(half-1); % define range of movie repeats for current half
                    mean_activity_per_half(:,:,half) = mean(current_mouse(:,:,current_half),3,'omitnan'); % average the neuronal activity across chosen movie repeats
                end
                
                % calculate mean pv correlation between sessions halves using only
                % cells that were active in both compared time points (either
                % within or between sessions)
                mean_pv_corr = []; % define empty variable that will store pv corr values for all mice
                for halfA = 1:size(mean_activity_per_half,3) % loop over session halves
                    halfA_activity = mean_activity_per_half(:,:,halfA); % subset neuronal activity for a single half
                    
                    for halfB = 1:size(mean_activity_per_half,3) % loop over session halves
                        halfB_activity = mean_activity_per_half(:,:,halfB); % subset neuronal activity for a single half
                        
                        % valid_cells - cells that were active in both of compared time points
                        valid_cells = [mean(halfA_activity,2,'omitnan') ~= 0] & [mean(halfB_activity,2,'omitnan') ~= 0]; % true if average activity is above zero in both session halves
                        valid_halfA_activity = halfA_activity(valid_cells,:); % subset only the cells that met the requirments of 'valid_cells'
                        valid_halfB_activity = halfB_activity(valid_cells,:); % subset only the cells that met the requirments of 'valid_cells'
                        
                        pv_corr = corr(valid_halfA_activity,valid_halfB_activity); % calculate the pv correlation between time bins of different halves
                        mean_pv_corr(halfA,halfB) = mean(diag(pv_corr),'omitnan'); % calculate the mean pv corr across corresponding time bins
                    end
                end
                
                within_session_values = [mean_pv_corr(1,2),mean_pv_corr(3,4),mean_pv_corr(5,6)]; % pv corr values between halves of the same session (within session)
                proximal_sessions_values = [mean_pv_corr(1:2,3:4),mean_pv_corr(3:4,5:6)]; % pv corr values between halves of proximal sessions (sessions 1&2 and sessions 2&3)
                distal_session_values = mean_pv_corr(1:2,5:6); % pv corr values between halves of distal sessions (sessions 1&3)
                
                elapsed_sess_pv(mouse,1)=  mean(within_session_values(:),'omitnan'); % average pv corr values for within session;
                elapsed_sess_pv(mouse,2)= mean(proximal_sessions_values(:),'omitnan'); % average pv corr values for proximal sessions
                elapsed_sess_pv(mouse,3)= mean(distal_session_values(:),'omitnan'); % average pv corr values for distal sessions
                
                
            end
            
        end
        elapsed_sess_pv_depth{area,depth} = elapsed_sess_pv; % store mean pv corr values for all mice of a given area
    end
end

figure % visualization of pv corr between sessions across imaging depth
for depth = 1:size(depth_list,1) % loop over layers
    curent_depth = cell2mat(elapsed_sess_pv_depth(:,depth)); % subset pv corr values for specific cortical layer
    
    mean_stability = mean(curent_depth,'omitnan'); % calculate mean pv across mice
    std_stability = std(curent_depth,'omitnan'); % calculate standard deviation across mice
    ste_stability=std_stability./sqrt(size(curent_depth,1)); % calculate standard error across mice
    
    hold on
    errorbar(mean_stability ,ste_stability ,'o','color',new_jet_colormap(depth+20*(depth-1),:),...
        'markerfacecolor',new_jet_colormap(depth+20*(depth-1),:),'capsize',0,'linestyle','-','linewidth',3,'markersize',3)
    
end
ylabel('PV correlation')
set(gca,'xtick',[1:3],'xticklabel',{'Within session','Proximal sessions','Distal sessions'},...
    'ytick',[0.3:0.1:0.6])
xtickangle(15)
xlim([0.5 3.5])
lgd = legend({'L2/3','L4','L5','L6'});
legend('boxoff');
lgd.Position = [0.2 0.2 0.2 0.2];
ylim([0.25 0.65])


%% Figure 4C - PV similarity index across layers

depth_list = [150 250; 251 350; 351 500; 501 700]; % define imaging depth ranges for each layer
layer_list = {'L2/3','L4  ','L5  ','L6  '}; %  define layer names

nat_movie = 1; % natural movie 1
cell_cutoff = 20; % cell count threshold of 20 cells

valid_mice_num = []; % define empty variable that will store the number of FOV in each layer (for visualization)
elapsed_sess_pv_depth = {}; % define empty variable that will store the pv corr values between sessions across mice and visual areas
for area = 1:6 % loop over areas
    
    for depth = 1:size(depth_list,1) % loop over layers
        
        % valid_mice_num - true only if the number of cells recorded is
        % atleast 20 in each of the recording sessions
        valid_cell_num = min(calcium_excitatory_cell_count{area},[],2) >= cell_cutoff;
        
        % valid_depth - true only if the imaging depth of the FOV is
        % within the specific imaging depth range
        valid_depth = calcium_excitatory_imaging_depth{area} >= depth_list(depth,1) & calcium_excitatory_imaging_depth{area} <= depth_list(depth,2);
        
        % valid_mice - mice that passed the reqirements of both
        % 'valid_mice_num' and 'valid_depth'
        valid_mice = valid_cell_num & valid_depth';
        valid_mice_num(area,depth) = sum(valid_mice); % store the number of FOV that passed the requirments of 'valid mice'
        
        elapsed_sess_pv = []; % define empty variable that will store the pv corr values across mice
        if valid_mice_num(area) >= 0 % test if there are any FOV that passed the criteria of 'valid_mice'
            
            current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie); % subsetthe neuronal activity of mice that paseed the the criteria of 'valid_mice'
            
            for mouse = 1:size(current_area,1) % loop over mice
                clc;
                disp(['Calculating PV correlation across layers for Calcium imaging data:'])
                disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Depth: ',layer_list{depth},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
                
                current_mouse = current_area{mouse}; % subset neuronal activity for a single mouse (#cells by 30 time bins by 30 movie repeats)
                
                %calculate the average activity rate across movie repeats for each session
                mean_activity_per_sess = []; % define an empty variable that will store the averaged neuronal responses of each session
                for sess = 1:3 % loop over session halves
                    current_sess = [1:10] + 10*(sess-1); % define range of movie repeats for current session
                    mean_activity_per_sess(:,:,sess) = mean(current_mouse(:,:,current_sess),3,'omitnan'); % average the neuronal activity across chosen movie repeats
                end
                
                % calculate mean pv correlation between sessions using only
                % cells that were active in both compared time points
                mean_pv_corr = []; % define empty variable that will store pv corr values for all mice
                for sessA = 1:size(mean_activity_per_sess,3) % loop over sessions
                    sessA_activity = mean_activity_per_sess(:,:,sessA); % subset neuronal activity for a single session
                    
                    for sessB = 1:size(mean_activity_per_sess,3) % loop over sessions
                        sessB_activity = mean_activity_per_sess(:,:,sessB); % subset neuronal activity for a single session
                        
                        % valid_cells - cells that were active in both of compared time points
                        valid_cells = [mean(sessA_activity,2,'omitnan') ~= 0] & [mean(sessB_activity,2,'omitnan') ~= 0]; % true if average activity is above zero in both sessions
                        valid_sessA_activity = sessA_activity(valid_cells,:); % subset only the cells that met the requirments of 'valid_cells'
                        valid_sessB_activity = sessB_activity(valid_cells,:); % subset only the cells that met the requirments of 'valid_cells'
                        
                        pv_corr = corr(valid_sessA_activity,valid_sessB_activity); % calculate the pv correlation between time bins of different sessions
                        mean_pv_corr(sessA,sessB) = nanmean(diag(pv_corr)); % calculate the mean pv corr across corresponding time bins
                    end
                end
                
                proximal_sessions = mean(diag(mean_pv_corr,1),'omitnan'); % pv corr values between proximal sessions (sessions 1&2 and sessions 2&3)
                distal_sessions = diag(mean_pv_corr,2); % pv corr values between distal sessions (sessions 1&3)
                
                elapsed_sess_pv(mouse,:) = [proximal_sessions,distal_sessions]; % store pv corr values of both proximal and distal values for all mice
                
            end
            
        end
        elapsed_sess_pv_depth{area,depth} = elapsed_sess_pv; % store mean pv corr values for all mice of a given area
    end
end

% define empty variables that will store the results of Wilcoxon signed-rank tests
pvalues = []; % define empty variable for the pvalues

mean_stability = []; % define empty variable that will store the mean pv similarity score values across mice for each layer
ste_stability = []; % define empty variable that will store the standard error for the pv similarity score values across mice for e
for depth = 1:4 % loop over cortical layers
    current_depth = cell2mat(elapsed_sess_pv_depth(:,depth)); % subset pv corr values for specific layer
    [pvalues(depth),~,stats] = signrank(current_depth(:,1),current_depth(:,2),'tail','right'); % perform one-sided wilcoxon signed-rank test for difference between proximal and distal sessions
    
    current_depth(current_depth<0) = 0; % rectify negative pv corr values into zero
    
    % calculate population vector similarity score for each mouse defined
    % as: [distal - proximal] / [distal + proximal] (ranges between -1 and 1)
    pv_similarity_score = [current_depth(:,2)-current_depth(:,1)]./[current_depth(:,2)+current_depth(:,1)];
    
    mean_stability(depth) = mean(pv_similarity_score,'omitnan'); % calculate mean pv similarity score across mice
    std_stability = std(pv_similarity_score,'omitnan'); % calculate the standard deviation across mice
    ste_stability(depth) = std_stability./sqrt(length(pv_similarity_score)); % calculate standard error across mice
    
end

figure % visualize pv similarity score between proximal as distal sessions across cortical layers
hold on
xlim([0.5 4.5])
for depth = 1:4
    errorbar(depth,mean_stability(depth),ste_stability(depth),'o','markerfacecolor',new_jet_colormap(depth+20*(depth-1),:),'color',new_jet_colormap(depth+20*(depth-1),:),'capsize',0,'linewidth',2,'linestyle','none')
end
plot(xlim,[0 0],'--','color',[0.4 0.4 0.4])
ylim([-0.2 0.05])
ylabel('Population vector similarity index')
set(gca,'xtick',[1:4],'xticklabel',{'L2/3','L4','L5','L6'})

% correct for multiple comparisons using bonferroni-holm method
corrected_pvalues = bonf_holm(pvalues);

% define statistics summary table
VarNames = {'depth','pvalue','bonf_holm'};
statistics = table(cell2mat(layer_list'),pvalues(:),corrected_pvalues(:),'VariableNames',VarNames);

% display statistics summary table
clc;
disp(['PV correlation between proximal sessions compared to distal sessions'])
disp(['One-tailed Wilcoxon signed-rank tests with Holm�Bonferroni correction:'])
disp(statistics)

%% Figure 4D - SST cells visual fields across movie repeats
% visualize the neuronal responses of three exmaple calcium imaging SST cells
% during the presentation of 'Natural movie 1'.

nat_movie = 1; % natural movie 1
area = 1; % area V1
mouse = 9; % exmaple mouse #9

% subset the neuronal activity of a single exmaple mouse (sst cre line)
example_mouse = calcium_inhibitory_population_vectors{area}{mouse,nat_movie}*30;
cell_list = [5,4,7];
figure('units','normalized','position',[0.3 0.3 0.3 0.5])
for current_cell = 1:length(cell_list)
    % example_cell - binned activity rate (calcium events) for example cell
    % (30 movie repeats x 30 time bins)
    example_cell = squeeze(example_mouse(cell_list(current_cell),:,:))';
    
    % gaussian smoothing of neuronal activity for each movie repeat
    smooth_example_cell = [];
    for repeat = 1:size(example_cell,1) % loop over movie repeats
        smooth_example_cell(repeat,:) = imgaussfilt(example_cell(repeat,:),2); % sigma = 2;
    end
    
    % normalize neuronal activity for each movie repeat based on peak firing rate
    norm_example_cell = smooth_example_cell ./ max(smooth_example_cell,[],2);
    
    % visualize normalized activity patterns for each cell
    subplot(2,3,current_cell)
    imagesc(norm_example_cell)
    hold on
    plot(xlim, [10 10]+0.5,'--','linewidth',2,'color','w')
    plot(xlim, [20 20]+0.5,'--','linewidth',2,'color','w')
    if current_cell ==1
        ylabel('Movie repeat')
        text(0.575, 0.925,['Sess 1'],'Units','normalized','color','w')
        text(0.575, 0.6,['Sess 2'],'Units','normalized','color','w')
        text(0.575, 0.275,['Sess 3'],'Units','normalized','color','w')
    elseif current_cell == 3
        cb = colorbar;
        set(cb,'position',[0.925 0.55 0.03 0.32])
        
    end
    colormap(newmap3)
    title(['Cell #',num2str(cell_list(current_cell))])
    
    % calculate mean activity rate across movie repeats for each session
    mean_example_cell = [];
    mean_example_cell(1,:) = mean(example_cell(1:10,:),1,'omitnan'); % average activity for session 1 (repeats 1-10)
    mean_example_cell(2,:) = mean(example_cell(11:20,:),1,'omitnan'); % average activity for session 2 (repeats 11-20)
    mean_example_cell(3,:) = mean(example_cell(21:30,:),1,'omitnan'); % average activity for session 3 (repeats 21-30)
    
    % visualize average activity patterns (tuning curve) for session of each cell
    subplot(2,3,current_cell+3)
    hold on
    plot(mean_example_cell(1,:),'-','markersize',20,'linewidth',1.5,'color',[0.3 0.3 0.3])
    plot(mean_example_cell(2,:)','-','markersize',20,'linewidth',1.5,'color',[0.5 0.5 0.5])
    plot(mean_example_cell(3,:)','-','markersize',20,'linewidth',1.5,'color',[0.7 0.7 0.7])
    
    xlim([1 30])
    if current_cell ==1
        ylabel('Mean activity rate (events/sec)')
        legend({'Sess 1','Sess 2','Sess 3'})
        legend('boxoff')
    elseif current_cell ==2
        xlabel('Time in movie (sec)')
    end
    
end
suptitle('SST interneurons:')

%% Figure 4E - VIP cells visual fields across movie repeats
% visualize the neuronal responses of three exmaple calcium imaging VIP cells
% during the presentation of 'Natural movie 1'.

nat_movie = 1; % natural movie 1
area = 1; % area V1
mouse = 23; % exmaple mouse #23

% subset the neuronal activity of a single exmaple mouse (vip cre line)
example_mouse = calcium_inhibitory_population_vectors{area}{mouse,nat_movie}*30;

cell_list = [4,10,16];
figure('units','normalized','position',[0.3 0.3 0.3 0.5])
for current_cell = 1:length(cell_list)
    % example_cell - binned activity rate (calcium events) for example cell
    % (30 movie repeats x 30 time bins)
    example_cell = squeeze(example_mouse(cell_list(current_cell),:,:))';
    
    % gaussian smoothing of neuronal activity for each movie repeat
    smooth_example_cell = [];
    for repeat = 1:size(example_cell,1) % loop over movie repeats
        smooth_example_cell(repeat,:) = imgaussfilt(example_cell(repeat,:),2); % sigma = 2;
    end
    
    % normalize neuronal activity for each movie repeat based on peak firing rate
    norm_example_cell = smooth_example_cell ./ max(smooth_example_cell,[],2);
    
    % visualize normalized activity patterns for each cell
    subplot(2,3,current_cell)
    imagesc(norm_example_cell)
    hold on
    plot(xlim, [10 10]+0.5,'--','linewidth',2,'color','w')
    plot(xlim, [20 20]+0.5,'--','linewidth',2,'color','w')
    if current_cell ==1
        ylabel('Movie repeat')
        text(0.575, 0.925,['Sess 1'],'Units','normalized','color','w')
        text(0.575, 0.6,['Sess 2'],'Units','normalized','color','w')
        text(0.575, 0.275,['Sess 3'],'Units','normalized','color','w')
    elseif current_cell == 3
        cb = colorbar;
        set(cb,'position',[0.925 0.55 0.03 0.32])
        
    end
    colormap(newmap3)
    title(['Cell #',num2str(cell_list(current_cell))])
    
    % calculate mean activity rate across movie repeats for each session
    mean_example_cell = [];
    mean_example_cell(1,:) = mean(example_cell(1:10,:),1,'omitnan'); % average activity for session 1 (repeats 1-10)
    mean_example_cell(2,:) = mean(example_cell(11:20,:),1,'omitnan'); % average activity for session 2 (repeats 11-20)
    mean_example_cell(3,:) = mean(example_cell(21:30,:),1,'omitnan'); % average activity for session 3 (repeats 21-30)
    
    % visualize average activity patterns (tuning curve) for session of each cell
    subplot(2,3,current_cell+3)
    hold on
    plot(mean_example_cell(1,:),'-','markersize',20,'linewidth',1.5,'color',[0.3 0.3 0.3])
    plot(mean_example_cell(2,:)','-','markersize',20,'linewidth',1.5,'color',[0.5 0.5 0.5])
    plot(mean_example_cell(3,:)','-','markersize',20,'linewidth',1.5,'color',[0.7 0.7 0.7])
    
    xlim([1 30])
    if current_cell ==1
        ylabel('Mean activity rate (events/sec)')
        legend({'Sess 1','Sess 2','Sess 3'})
        legend('boxoff')
    elseif current_cell ==2
        xlabel('Time in movie (sec)')
    end
    
end
suptitle('VIP interneurons:')

%% Figure 4F - Pvalb cells visual fields across movie repeats
% visualize the neuronal responses of three exmaple calcium imaging PVALB cells
% during the presentation of 'Natural movie 1'.

nat_movie = 1; % natural movie 1
area = 1; % area V1
mouse = 38; % exmaple mouse #38

% subset the neuronal activity of a single exmaple mouse (sst cre line)
example_mouse = calcium_inhibitory_population_vectors{area}{mouse,nat_movie}*30;
cell_list = [1,15,16];
figure('units','normalized','position',[0.3 0.3 0.3 0.5])
for current_cell = 1:length(cell_list)
    % example_cell - binned activity rate (calcium events) for example cell
    % (30 movie repeats x 30 time bins)
    example_cell = squeeze(example_mouse(cell_list(current_cell),:,:))';
    
    % gaussian smoothing of neuronal activity for each movie repeat
    smooth_example_cell = [];
    for repeat = 1:size(example_cell,1) % loop over movie repeats
        smooth_example_cell(repeat,:) = imgaussfilt(example_cell(repeat,:),2); % sigma = 2;
    end
    
    % normalize neuronal activity for each movie repeat based on peak firing rate
    norm_example_cell = smooth_example_cell ./ max(smooth_example_cell,[],2);
    
    % visualize normalized activity patterns for each cell
    subplot(2,3,current_cell)
    imagesc(norm_example_cell)
    hold on
    plot(xlim, [10 10]+0.5,'--','linewidth',2,'color','w')
    plot(xlim, [20 20]+0.5,'--','linewidth',2,'color','w')
    if current_cell ==1
        ylabel('Movie repeat')
        text(0.575, 0.925,['Sess 1'],'Units','normalized','color','w')
        text(0.575, 0.6,['Sess 2'],'Units','normalized','color','w')
        text(0.575, 0.275,['Sess 3'],'Units','normalized','color','w')
    elseif current_cell == 3
        cb = colorbar;
        set(cb,'position',[0.925 0.55 0.03 0.32])
        
    end
    colormap(newmap3)
    title(['Cell #',num2str(cell_list(current_cell))])
    
    % calculate mean activity rate across movie repeats for each session
    mean_example_cell = [];
    mean_example_cell(1,:) = mean(example_cell(1:10,:),1,'omitnan'); % average activity for session 1 (repeats 1-10)
    mean_example_cell(2,:) = mean(example_cell(11:20,:),1,'omitnan'); % average activity for session 2 (repeats 11-20)
    mean_example_cell(3,:) = mean(example_cell(21:30,:),1,'omitnan'); % average activity for session 3 (repeats 21-30)
    
    % visualize average activity patterns (tuning curve) for session of each cell
    subplot(2,3,current_cell+3)
    hold on
    plot(mean_example_cell(1,:),'-','markersize',20,'linewidth',1.5,'color',[0.3 0.3 0.3])
    plot(mean_example_cell(2,:)','-','markersize',20,'linewidth',1.5,'color',[0.5 0.5 0.5])
    plot(mean_example_cell(3,:)','-','markersize',20,'linewidth',1.5,'color',[0.7 0.7 0.7])
    
    xlim([1 30])
    if current_cell ==1
        ylabel('Mean activity rate (events/sec)')
        legend({'Sess 1','Sess 2','Sess 3'})
        legend('boxoff')
    elseif current_cell ==2
        xlabel('Time in movie (sec)')
    end
    
end
suptitle('Pvalb interneurons:')

%% Figure 4G - Excitatory VS inhibitory - PV corr across movie repeats
% Calculate the difference in PV correlation as a function of time for the
% inhibitory and excitatory Cre lines imaged from areas V1, LM and PM

nat_movie = 1; % natural movie 1
num_repeats = 10; % number of movie repeats in each imaging session
cell_cutoff = 20; % cell count threshold of 20 cells
cell_type = {'Excitatory','Inhibitory'}; % define cell subtypes names
area_list = [1,2,4]; % areas included are V1,LM and PM

elapse_repeat_pv_corr_areas = {};
for subtype = 1:2 % loop over cell types
    for area = 1:3 % loop over visual areas
        if subtype == 1 % excitatory cre lines
            % valid_mice - true only if the number of cells recorded is
            % atleast 20 in each of the recording sessions
            valid_mice = min(calcium_excitatory_cell_count{area_list(area)},[],2) >= cell_cutoff;
            current_area = calcium_excitatory_population_vectors{area_list(area)}(valid_mice,nat_movie); % subset all mice that met the requirements of 'valid mice'
        elseif subtype == 2 % inhibitory cre lines
            current_area = calcium_inhibitory_population_vectors{area_list(area)}(:,nat_movie); % subset all mice from a specific area
        end
        
        elapse_repeat_pv_corr = []; % define an empty variable for average pv corr as a function of elapsed time across mice (size of #mice by 9)
        for mouse = 1:length(current_area) % loop over mice
            clc;
            disp(['Calculating PV correlation between movie repeats for calcium imaging data:'])
            disp(['Stimulus: Natural movie 1 | Cell type: ',cell_type{subtype},' | Area: ',brain_areas{area_list(area)},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
            
            current_mouse = current_area{mouse}; % subset the neuronal activity of a single mouse (size of #cells by 30 time bins by 30 movie repeats)
            
            current_mouse_sess1 =  current_mouse(:,:,1:10); % subset neuronal activity during the first session (repeats 1-10)
            valid_cells = mean(squeeze(mean(current_mouse_sess1,2,'omitnan')),2,'omitnan')>0; % include only cells that were active in at least one movie repeat
            valid_current_mouse_sess1 = current_mouse_sess1(valid_cells,:,:); % subset only the cells that passed the requirements of 'valid_cells'
            
            current_mouse_sess2 =  current_mouse(:,:,11:20); % subset neuronal activity during the second session (repeats 11-20)
            valid_cells = mean(squeeze(mean(current_mouse_sess2,2,'omitnan')),2,'omitnan')>0; % include only cells that were active in at least one movie repeat
            valid_current_mouse_sess2 = current_mouse_sess2(valid_cells,:,:); % subset only the cells that passed the requirements of 'valid_cells'
            
            current_mouse_sess3 =  current_mouse(:,:,21:30); % subset neuronal activity during the third session (repeats 21-30)
            valid_cells = mean(squeeze(mean(current_mouse_sess3,2,'omitnan')),2,'omitnan')>0; % include only cells that were active in at least one movie repeat
            valid_current_mouse_sess3 = current_mouse_sess3(valid_cells,:,:); % subset only the cells that passed the requirements of 'valid_cells'
            
            pv_corr_across_session = []; % define an empty variable for average pv corr across movie repeats
            for repeat1 = 1:10 % loop over movie repeats
                for repeat2 = 1:10 % loop over movie repeats
                    pv_corr = corr(valid_current_mouse_sess1(:,:,repeat1),valid_current_mouse_sess1(:,:,repeat2)); % PV corr between time bins of pair of movie repeats in session 1
                    pv_corr_across_session(repeat1,repeat2,1) = mean(diag(pv_corr),'omitnan'); % average PV corr across corresponding time bins
                    
                    pv_corr = corr(valid_current_mouse_sess2(:,:,repeat1),valid_current_mouse_sess2(:,:,repeat2)); % PV corr between time bins of pair of movie repeats in session 2
                    pv_corr_across_session(repeat1,repeat2,2) = mean(diag(pv_corr),'omitnan'); % average PV corr across corresponding time bins
                    
                    pv_corr = corr(valid_current_mouse_sess3(:,:,repeat1),valid_current_mouse_sess3(:,:,repeat2)); % PV corr between time bins of pair of movie repeats in session 3
                    pv_corr_across_session(repeat1,repeat2,3) = mean(diag(pv_corr),'omitnan'); % average PV corr across corresponding time bins
                    
                end
            end
            mean_pv_corr_across_mice = mean(pv_corr_across_session,3,'omitnan'); % average across sessions
            
            % calculate pv correlation as function of elapsed time
            for diagonal = 1:num_repeats-1 % loop over diagonals
                elapse_repeat_pv_corr(mouse,diagonal) = mean(diag(mean_pv_corr_across_mice,diagonal),'omitnan'); % average pv corr values across diagonals
            end
        end
        elapse_repeat_pv_corr_areas{subtype,area} = elapse_repeat_pv_corr; % store mean pv corr values across mice for each area
    end
end

% friedman test parameters
pvalues = []; % define empty variable for the pvalues
df = []; % define empty variable for the degrees of freedom
chi = []; % define empty variable for the chi square values

figure('units','normalized','position',[0.3 0.3 0.4 0.4]) % visualization of mean pv corr as function of time across areas and mice
for area = 1:3 % loop over areas
    plt = []; % store plot information of each cell type
    for subtype = 1:2 % loop over cell types
        current_area = elapse_repeat_pv_corr_areas{subtype,area}; % subset the mean pv corr values of specific combination of cell type and visual area
        
        [pvalues(subtype,area),tbl,stats] = friedman(current_area,[1],'off'); % perform friedman test for main effect of elapsed time
        df(subtype,area) = tbl{2,3}; % store degrees of freedom
        chi(subtype,area) = tbl{2,5}; % store chi square values
        
        norm_elapse_repeat_pv_corr = current_area-current_area(:,1); % scale the pv corr values of each mouse based on the minimal interval
        
        mean_norm_vals = mean(norm_elapse_repeat_pv_corr,'omitnan'); % calculate mean pv corr across mice
        std_norm_vals = std(norm_elapse_repeat_pv_corr,'omitnan'); % calculate standard deviation across mice
        ste_norm_vals = std_norm_vals./sqrt(size(current_area,1)); % calculate standard error across mice
        
        subplot(1,3,area)
        hold on
        x = [1:length(mean_norm_vals)]';
        y = mean_norm_vals';
        dy = ste_norm_vals';
        hold on
        if subtype == 1
            fill([x;flipud(x)],[y-dy;flipud(y+dy)],[0.8 0.8 0.8],'linestyle','none','facealpha',0.4);
            plt(1) = plot(mean_norm_vals,'color',[0.6 0.6 0.6],'linewidth',2);
        elseif subtype == 2
            fill([x;flipud(x)],[y-dy;flipud(y+dy)],colors(area_list(area),:),'linestyle','none','facealpha',0.4);
            plt(2) = plot(mean_norm_vals,'color',colors2(area_list(area),:),'linewidth',2);
        end
        xlim([0.5 9.5])
        ylim([-0.07 0])
        if area == 1
            ylabel('PV correlation difference')
        elseif area == 2
            xlabel('Elapsed time (# of movie repeats)')
            set(gca,'ytick',[])
        else
            set(gca,'ytick',[])
        end
    end
    text(0.8, 0.95,brain_areas{area_list(area)},'Units','normalized','color',[0 0 0 ],'fontsize',15)
    lgd = legend(plt,{'Excitatory','Inhibitory'},'location','southwest');
    legend('boxoff')
    %lgd.Position = [0.15 0.225 0.15 0.15];
    
end
suptitle('Seconds to minutes:                                         1 repeats = 30 seconds')

% correct for multiple comparisons using bonferroni-holm method
corrected_pval = bonf_holm(pvalues(2,:));

% define statistics summary table
VarNames = {'area','df','chi','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas([1,2,4])'),df(2,:)',chi(2,:)',pvalues(2,:)',corrected_pval(:),'VariableNames',VarNames);

% display statistics summary table
clc;
disp(['PV correlation as a function of elapsed time for inhibitory Cre lines'])
disp(['Friedman�s tests with Holm�Bonferroni correction:'])
disp(statistics)

%% Figure 4H - Excitatory VS inhibitory - PV corr across blocks

nat_movie = 2; % natural movie 3
cell_cutoff = 20; % cell count threshold of 20 cells
cell_type = {'Excitatory','Inhibitory'}; % define cell subtypes names
area_list = [1,2,4]; % areas included are V1,LM and PM

within_between_stability_area = {}; % define empty variable that will store all pv corr values across blocks across mice and visual areas
for subtype = 1:2 % loop over cell types
    for area = 1:3 % loop over visual areas
        if subtype == 1 % excitatory cre lines
            % valid_mice - true only if the number of cells recorded is
            % atleast 20 in session A
            valid_mice = calcium_excitatory_cell_count{area_list(area)}(:,1) >= cell_cutoff;
            current_area = calcium_excitatory_population_vectors{area_list(area)}(valid_mice,nat_movie); % subset all mice that met the requirements of 'valid mice'
        elseif subtype == 2 % inhibitory cre lines
            current_area = calcium_inhibitory_population_vectors{area_list(area)}(:,nat_movie); % subset all mice
        end
        
        within_between_stability = []; % define empty variable that will store pv corr values for all mice
        for mouse = 1:length(current_area) % loop over mice
            clc;
            disp(['Calculating PV correlation between blocks for calcium imaging data:'])
            disp(['Stimulus: Natural movie 3 | Cell type: ',cell_type{subtype},' | Area: ',brain_areas{area_list(area)},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
            
            current_mouse = current_area{mouse}; % subset a single mouse (size of #cells by 30 time bins by 10 movie repeats)
            
            current_mouse_blocks = []; % define an empty variable that will store the mean neuronal activty for each half of the two blocks
            current_mouse_blocks(:,:,1) = mean(current_mouse(:,:,1:2),3,'omitnan'); % calculate the mean activity across repreats for the first half of block A (repeats 1-2)
            current_mouse_blocks(:,:,2) = mean(current_mouse(:,:,3:5),3,'omitnan'); % calculate the mean activity across repreats for the second half of block A (repeats 3-5)
            current_mouse_blocks(:,:,3) = mean(current_mouse(:,:,6:7),3,'omitnan'); % calculate the mean activity across repreats for the first half of block B (repeats 6-7)
            current_mouse_blocks(:,:,4) = mean(current_mouse(:,:,8:10),3,'omitnan'); % calculate the mean activity across repreats for the second half of block B (repeats 8-10)
            
            % calculate pv correlation between the four halves of the two blocks
            mean_pv_corr_across_mice = []; % define empty variable that will store the mean pv across halves for a given mouse
            for half1 = 1:4 % loop over blocks halves
                for half2 = 1:4 % loop over blocks halves
                    current_pv = corr(current_mouse_blocks(:,:,half1),current_mouse_blocks(:,:,half2)); % calculate the pv corr between time bins of different halves
                    mean_pv_corr_across_mice(half1,half2) = mean(diag(current_pv),'omitnan'); % calculate the mean pv corr across corresponding time bins
                    
                end
            end
            
            % calculate and store the average pv corr between halves of the
            % same block ('within-block') and average pv corr between different
            % halves of different blocks ('between-blocks')
            within_block = [mean_pv_corr_across_mice(1,2),mean_pv_corr_across_mice(3,4)]; % subset pv corr between halves of the same block
            between_blocks = mean_pv_corr_across_mice(1:2,3:4); % subset pv corr between halves of different blocks
            
            within_between_stability(mouse,1) = mean(within_block(:),'omitnan'); % calculate the average pv corr within a block
            within_between_stability(mouse,2) = mean(between_blocks(:),'omitnan'); % calculate the average pv corr between blocks
            
        end
        within_between_stability_area{subtype,area} = within_between_stability; % store mean pv corr across mice for each area
    end
end

% define empty variables that will store the results of Wilcoxon signed-rank tests
pvalue = []; % define empty variable for the pvalues
zvalue = []; % define empty variable for the z values
figure('units','normalized','position',[0.3 0.3 0.4 0.4]) % visualization of pv corr between blocks across areas and mice
for area = 1:3 % loop over visual areas
    plt = []; % store plot information for each cell type
    for subtype = 1:2 % loop over cell types
        current_area = within_between_stability_area{subtype,area}; % calculate mean pv corr across mice
        
        [pvalues(subtype,area),~,stats] = signrank(current_area(:,1),current_area(:,2)); % perform two-sided wilcoxon signed-rank test for difference within-between blocks
        zvalues(subtype,area) = stats.zval; % store z values
        
        mean_stability = mean(current_area,'omitnan'); % calculate mean pv corr across mice
        std_stability = std(current_area,'omitnan'); % calculate standard deviation across mice
        ste_stability = std_stability./sqrt(size(current_area,1)); % calculate standard error across mice
        
        subplot(1,3,area)
        hold on
        if subtype == 1
            plt(1) = errorbar(mean_stability,ste_stability,'o','color',[0.7 0.7 0.7],...
                'markerfacecolor',[0.7 0.7 0.7],'capsize',0,'linestyle','-','linewidth',3,'markersize',5);
        elseif subtype ==2
            plt(2) = errorbar(mean_stability,ste_stability,'o','color',colors(area_list(area),:),...
                'markerfacecolor',colors(area_list(area),:),'capsize',0,'linestyle','-','linewidth',3,'markersize',5);
        end
        
        xlim([0.5 2.5])
        ylim([0.3 0.6])
        
        set(gca,'xtick',1:2,'xticklabel',{'Within block','Between blocks'})
        xtickangle(15)
    end
    
    if area == 1
        ylabel('PV correlation')
    else
        set(gca,'ytick',[])
    end
    
    text(0.8, 0.95,brain_areas{area_list(area)},'Units','normalized','color',[0 0 0 ],'fontsize',15)
    lgd = legend(plt,{'Excitatory','Inhibitory'},'location','southwest');
    legend('boxoff')
end
suptitle('Minutes to hours:')

% correct for multiple comparisons using bonferroni-holm method
corrected_pval = bonf_holm(pvalues(2,:));

% define statistics summary table
VarNames = {'area','zvalue','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas([1,2,4])'),zvalues(2,:)',pvalues(2,:)',corrected_pval(:),'VariableNames',VarNames);

% display statistics summary table
clc;
disp(['PV correlation within blocks compared to between blocks for inhibitory Cre lines'])
disp(['Two-tailed Wilcoxon signed-rank tests with Holm�Bonferroni correction:'])
disp(statistics)

%% Figure 4I - Excitatory VS inhibitory - PV corr across sessions
% Calculate the PV correlation between the two halves of the same session,
% between halves of two temporally proximal sessions and between halves of
% two temporally distal sessions for the inhibitory and excitatory Cre lines

nat_movie = 1; % natural movie 1
cell_cutoff = 20; % cell count threshold of 20 cells
cell_type = {'Excitatory','Inhibitory'}; % define cell subtypes names
area_list = [1,2,4]; % areas included are V1,LM and PM

within_between_session_stability_areas = {}; % define empty variable that will store the pv corr values of within and between sessions across mice and visual areas
for subtype = 1:2 % loop over cell types
    for area = 1:3 % loop over visual areas
        if subtype == 1 % excitatory cre lines
            % valid_mice - true only if the number of cells recorded is
            % atleast 20 in each of the recording sessions
            valid_mice = min(calcium_excitatory_cell_count{area_list(area)},[],2) >= cell_cutoff;
            current_area = calcium_excitatory_population_vectors{area_list(area)}(valid_mice,nat_movie); % mice that passed the requirments of 'valid_mice'
        elseif subtype == 2 % inhibitory cre lines
            current_area = calcium_inhibitory_population_vectors{area_list(area)}(:,nat_movie); % subset all mice
        end
        
        within_between_session_stability = []; % define empty variable that will store the pv correlation values for all mice
        for mouse = 1:length(current_area) % loop over mice
            clc;
            disp(['Calculating PV correlation between sessions for calcium imaging data:'])
            disp(['Stimulus: Natural movie 1 | Cell type: ',cell_type{subtype},' | Area: ',brain_areas{area_list(area)},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
            
            current_mouse = current_area{mouse}; % subset neuronal activity for a single mouse (#cells by 30 time bins by 30 movie repeats)
            
            %calculate the average activity rate for session halves across movie repeats
            mean_activity_per_half = []; % define an empty variable that will store the averaged neuronal responses of each session half
            for half = 1:6 % loop over session halves
                current_half = [1:5] + 5*(half-1); % define range of movie repeats for current half
                mean_activity_per_half(:,:,half) = mean(current_mouse(:,:,current_half),3,'omitnan');  % average the neuronal activity across chosen movie repeats
            end
            
            
            % calculate mean pv correlation between sessions halves using only
            % cells that were active in both compared time points (either
            % within or between sessions)
            mean_pv_corr = []; % define empty variable that will store pv corr values for all mice
            for halfA = 1:size(mean_activity_per_half,3) % loop over session halves
                halfA_activity = mean_activity_per_half(:,:,halfA); % subset neuronal activity for a single half
                
                for halfB = 1:size(mean_activity_per_half,3) % loop over session halves
                    halfB_activity = mean_activity_per_half(:,:,halfB); % subset neuronal activity for a single half
                    
                    % valid_cells - cells that were active in at least one of the compared time points
                    valid_cells = [mean(halfA_activity,2,'omitnan') ~= 0] | [mean(halfB_activity,2,'omitnan') ~= 0]; % true if average activity is above zero in at least one of the session halves
                    valid_halfA_activity = halfA_activity(valid_cells,:); % subset only the cells that met the requirments of 'valid_cells'
                    valid_halfB_activity = halfB_activity(valid_cells,:); % subset only the cells that met the requirments of 'valid_cells'
                    
                    pv_corr = corr(valid_halfA_activity,valid_halfB_activity); % calculate the pv correlation between time bins of different halves
                    mean_pv_corr(halfA,halfB) = mean(diag(pv_corr),'omitnan'); % calculate the mean pv corr across corresponding time bins
                    
                end
            end
            
            within_session_values = [mean_pv_corr(1,2),mean_pv_corr(3,4),mean_pv_corr(5,6)]; % pv corr values between halves of the same session (within session)
            proximal_sessions_values = [mean_pv_corr(1:2,3:4),mean_pv_corr(3:4,5:6)]; % pv corr values between halves of proximal sessions (sessions 1&2 and sessions 2&3)
            distal_session_values = mean_pv_corr(1:2,5:6); % pv corr values between halves of distal sessions (sessions 1&3)
            
            within_between_session_stability(mouse,1) = mean(within_session_values(:),'omitnan'); % average pv corr values for within session
            within_between_session_stability(mouse,2) = mean(proximal_sessions_values(:),'omitnan'); % average pv corr values for proximal sessions
            within_between_session_stability(mouse,3) = mean(distal_session_values(:),'omitnan'); % average pv corr values for distal sessionsend
        end
        within_between_session_stability_areas{subtype,area} = within_between_session_stability; % store mean pv corr values for all mice of a given area and cell type
    end
end

% define empty variables that will store the results of Wilcoxon signed-rank tests
pvalue = []; % define empty variable for the pvalues
zvalue = []; % define empty variable for the z values
figure('units','normalized','position',[0.3 0.3 0.4 0.4]) % visualization of pv corr between sessions across areas and mice
for area = 1:3 % loop over visual areas
    plt = []; % store plot information for each cell type
    for subtype = 1:2
        current_area = within_between_session_stability_areas{subtype,area}; % subset pv corr values for specific area and cell type
        
        [pvalues(subtype,area),~,stats] = signrank(current_area(:,2),current_area(:,3),'tail','right'); % perform one-sided wilcoxon signed-rank test for difference between proximal and distal sessions
        zvalues(subtype,area) = stats.zval; % store z values
        
        mean_stability = mean(current_area,'omitnan'); % calculate mean pv across mice
        std_stability = std(current_area,'omitnan'); % calculate standard deviation across mice
        ste_stability = std_stability./sqrt(size(current_area,1)); % calculate standard error across mice
        
        subplot(1,3,area)
        hold on
        if subtype == 1
            plt(1) = errorbar(mean_stability,ste_stability,'o','color',[0.7 0.7 0.7],...
                'markerfacecolor',[0.7 0.7 0.7],'capsize',0,'linestyle','-','linewidth',3,'markersize',5);
        elseif subtype ==2
            plt(2) = errorbar(mean_stability,ste_stability,'o','color',colors(area_list(area),:),...
                'markerfacecolor',colors(area_list(area),:),'capsize',0,'linestyle','-','linewidth',3,'markersize',5);
        end
        
        xlim([0.5 3.5])
        ylim([0.175 0.625])
        
        set(gca,'xtick',1:3,'xticklabel',{'Within block','Proximal sessions','Distal sessions'})
        xtickangle(15)
    end
    
    if area == 1
        ylabel('PV correlation')
    else
        set(gca,'ytick',[])
    end
    
    text(0.8, 0.95,brain_areas{area_list(area)},'Units','normalized','color',[0 0 0 ],'fontsize',15)
    lgd = legend(plt,{'Excitatory','Inhibitory'},'location','southwest');
    legend('boxoff')
end
suptitle('Minutes to days:')

% correct for multiple comparisons using bonferroni-holm method
corrected_pval = bonf_holm(pvalues(2,:));

% define statistics summary table
VarNames = {'area','zvalue','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas([1,2,4])'),zvalues(2,:)',pvalues(2,:)',corrected_pval(:),'VariableNames',VarNames);

% display statistics summary table
clc;
disp(['PV correlation  of proximal sessions compared to distal sessions for inhibitory Cre lines'])
disp(['One-tailed Wilcoxon signed-rank tests with Holm�Bonferroni correction:'])
disp(statistics)

%% Figure 5A - Visual hierarchy - ensemble rate correlation between movie repeats - Neuropixels

nat_movie = 1; % natural movie 1
cell_cutoff = 15; % cell count threshold of 15 units
num_repeats = 30; % number of movie repeats in each block

area_list = [1,2,7,8]; % deifine a list of visual areas to be compared (area V1, LM dLGN and LP)

elapse_repeat_rate_corr_area = {}; % define an empty variable for ensemble rate corr across areas

% calculate for each visual area of each mouse the ensemble rate correlation between
% pairs of movie repeats as a function of elapsed time
for area = 1:4 % loop over brain areas
    % valid_mice - only mice from the 'Functional connectivity' group with at least 15 units
    valid_mice = neuropixels_cell_count(:,area_list(area),nat_movie) >= cell_cutoff & movie_repeats(:,nat_movie) == num_repeats;
    
    current_area = neuropixels_population_vectors(valid_mice,area_list(area),nat_movie); % subset all mice that met the requirements of 'valid mice'
    elapse_repeat_rate_corr = []; % define an empty variable for ensemble rate corr matrices across mice (size of #mice by 29)
    
    for mouse = 1:size(current_area,1) % loop over mice
        clc;
        disp(['Calculating ensemble rate correlation between repeats for Neuropixels data:'])
        disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area_list(area)},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse};  % subset the neuronal activity of a single mouse (size of #units by 30 time bins by 60 movie repeats)
        current_mouse_blockA = squeeze(mean(current_mouse(:,:,1:30),2,'omitnan')); % subset and calculate mean neuronal activity in each movie repeat during the first block (repeats 1-30)
        current_mouse_blockB = squeeze(mean(current_mouse(:,:,31:60),2,'omitnan')); % subset and calculate mean neuronal activity in each movie repeat during the second block (repeats 31-60)
        
        % calculate ensemble rate correlation between pairs of movie repeats in each block
        mean_rate_corr = [];
        mean_rate_corr(:,:,1) = corr(current_mouse_blockA); % ensemble rate correlation for block A
        mean_rate_corr(:,:,2) = corr(current_mouse_blockB); % ensemble rate correlation for block B
        
        mean_rate_corr(mean_rate_corr<0) = 0; % rectify negative ensemble rate corr values to zero
        mean_rate_corr = mean(mean_rate_corr,3,'omitnan');  % average ensemble rate correlation across blocks
        
        % calculate ensemble rate correlation as function of elapsed time
        for diagonal = 1:num_repeats-1 % loop over diagonals
            elapse_repeat_rate_corr(mouse,diagonal) = mean(diag(mean_rate_corr,diagonal),'omitnan'); % average ensemble rate corr values across diagonals
        end
    end
    elapse_repeat_rate_corr_area{area} = elapse_repeat_rate_corr; % store ensemble rate corr values across mice for each area
end

% define color schemes for visualization
new_colors = [0 0.7 0.8;0.6 0.6 0.6;0.6 0.2 0.6;0.4 0.4 0.4];
new_colors2 = [0 0.5 0.6;0.5 0.5 0.5;0.5 0.1 0.4;0.3 0.3 0.3];
plt = []; % define empty variable that will store plot information for each area
rate_similarity_index = {}; % define empty variable that will store ensemble rate similarity index values of indevidual mice across areas
figure('units','normalized','position',[0.3 0.3 0.25 0.4]) % visualization of ensemble rate similarity index across areas
for area = 1:4 % loop over areas
    current_area = elapse_repeat_rate_corr_area{area}; % subset ensemble rate corr values of specific visual area
    rate_similarity_index{area} = [current_area-current_area(:,1)]./[current_area+current_area(:,1)]; % calculate for each mouse its ensemble rate similarity index
    
    mean_elapse_repeat = mean(rate_similarity_index{area},'omitnan'); % calculate mean ensemble rate similarity index across mice
    std_elapse_repeat = std(rate_similarity_index{area},'omitnan'); % calculate standard deviation across mice
    ste_elapse_repeat = std_elapse_repeat./sqrt(size(current_area,1)); % calculate standard error across mice
    
    x = [1:length(mean_elapse_repeat)]';
    y = mean_elapse_repeat';
    dy = ste_elapse_repeat';
    
    hold on
    fill([x;flipud(x)],[y-dy;flipud(y+dy)],new_colors(area,:),'linestyle','none','facealpha',0.4);
    plt(area) = plot(y,'color',new_colors2(area,:),'linewidth',3);
end
lgd = legend(plt,brain_areas(area_list));
legend('boxoff')
lgd.Position = [0.2 0.225 0.1 0.1];
text(0.05,0.075,'1 repeat = 30 seconds','Units','normalized','FontSize',11)
set(gca,'xtick',[1,10,20,29])
xlabel('Elapsed time (# of movie repeats)')
ylabel('Ensemble rate similarity index')
title('Neuropixels - seconds to minutes')


% statistical analysis for differences across pairs of visual areas
V1_rate_corr = rate_similarity_index{1};  % subset ensemble rate similarity indices of area V1
V2_rate_corr = rate_similarity_index{2}; % subset ensemble rate similarity indices of area LM

LGN_rate_corr = rate_similarity_index{3}; % subset ensemble rate similarity indices of area dLGN
LP_rate_corr = rate_similarity_index{4}; % subset ensemble rate similarity indices of area LP

% perform multiple Mann-Whitney rank-sum tests btween pairs of visual
% areas and visualize the results on the existing figure
for trial = 1:29 % loop over time intervals
    p = ranksum(V1_rate_corr(:,trial),V2_rate_corr(:,trial));
    if p < 0.05
        hold on
        plot(trial,0,'*','color',new_colors(1,:),'markersize',8)
    end
    
    p = ranksum(LGN_rate_corr(:,trial),LP_rate_corr(:,trial)); % perform two-sided Mann-Whitney rank-sum between area dLGN and area LP
    if p < 0.05
        hold on
        plot(trial,0.005,'*','color',new_colors(3,:),'markersize',8)
    end
end


%% Figure 5B - Visual hierarchy - tuning curve correlation between movie repeats - Neuropixels

nat_movie = 1; % natural movie 1
cell_cutoff = 15; % cell count threshold of 15 units
num_repeats = 30; % number of movie repeats in each block

area_list = [1,2,7,8]; % deifine a list of visual areas to be compared (area V1, LM dLGN and LP)

elapse_repeat_tuning_corr_area = {}; % define an empty variable for tuning curve corr across areas

% calculate for each visual area of each mouse the tuning curve correlation between
% pairs of movie repeats as a function of elapsed time
for area = 1:4 % loop over brain areas
    % valid_mice - only mice from the 'Functional connectivity' group with at least 15 units
    valid_mice = neuropixels_cell_count(:,area_list(area),nat_movie) >= cell_cutoff & movie_repeats(:,nat_movie) == num_repeats;
    
    current_area = neuropixels_population_vectors(valid_mice,area_list(area),nat_movie); % subset all mice that met the requirements of 'valid mice'
    elapse_repeat_tuning_corr = []; % define an empty variable for ensemble rate corr matrices across mice (size of #mice by 29)
    
    for mouse = 1:size(current_area,1) % loop over mice
        clc;
        disp(['Calculating tuning curve correlation between repeats for Neuropixels data:'])
        disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area_list(area)},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse}; % subset the neuronal activity of a single mouse (size of #units by 30 time bins by 60 movie repeats)
        current_mouse_blockA = current_mouse(:,:,1:30); % subset the neuronal activity in each movie repeat during the first block (repeats 1-30)
        current_mouse_blockB = current_mouse(:,:,31:60); % subset the neuronal activity in each movie repeat during the second block (repeats 31-60)
        
        % calculate the medina tuning curve correlation between pairs of movie repeats in each block
        mean_tuning_corr = []; % define an empty variale to store the median tuning curve corr for each block
        for repeat1 = 1:num_repeats % loop over movie repeats
            for repeat2 = 1:num_repeats % loop over movie repeats
                tuning_corr_blockA = corr(current_mouse_blockA(:,:,repeat1)',current_mouse_blockA(:,:,repeat2)'); % calculate the tuning curve correlation between pair of movie repeats in block A
                mean_tuning_corr(repeat1,repeat2,1) = median(diag(tuning_corr_blockA),'omitnan');  % calculate the median tuning curve correlation across corresponding units
                
                tuning_corr_blockB = corr(current_mouse_blockB(:,:,repeat1)',current_mouse_blockB(:,:,repeat2)'); % calculate the tuning curve correlation between pair of movie repeats in block B
                mean_tuning_corr(repeat1,repeat2,2) = median(diag(tuning_corr_blockB),'omitnan'); % calculate the median tuning curve correlation across corresponding units
            end
        end
        mean_tuning_corr(mean_tuning_corr<0) = 0; % rectify negative tuning curve corr values to zero
        mean_tuning_corr = mean(mean_tuning_corr,3,'omitnan');  % average tuning curve correlation across blocks
        
        % calculate tuning curve correlation as function of elapsed time
        for diagonal = 1:num_repeats-1 % loop over diagonals
            elapse_repeat_tuning_corr(mouse,diagonal) = mean(diag(mean_tuning_corr,diagonal),'omitnan'); % average tuning curve corr values across diagonals
        end
    end
    elapse_repeat_tuning_corr_area{area} = elapse_repeat_tuning_corr; % store tuning curve corr values across mice for each area
end

% define color schemes for visualization
new_colors = [0 0.7 0.8;0.6 0.6 0.6;0.6 0.2 0.6;0.4 0.4 0.4];
new_colors2 = [0 0.5 0.6;0.5 0.5 0.5;0.5 0.1 0.4;0.3 0.3 0.3];
plt = []; % define empty variable that will store plot information for each area
tuning_similarity_index = {}; % define empty variable that will store tuning curve similarity index values of indevidual mice across areas
figure('units','normalized','position',[0.3 0.3 0.25 0.4]) % visualization of tuning curve similarity index across areas
for area = 1:4
    current_area = elapse_repeat_tuning_corr_area{area}; % subset tuning curve corr values of specific visual area
    tuning_similarity_index{area} = [current_area-current_area(:,1)]./[current_area+current_area(:,1)]; % calculate for each mouse its tuning curve similarity index
    
    mean_elapse_repeat = mean(tuning_similarity_index{area},'omitnan'); % calculate mean tuning curve similarity index across mice
    std_elapse_repeat = std(tuning_similarity_index{area},'omitnan'); % calculate standard deviation across mice
    ste_elapse_repeat = std_elapse_repeat./sqrt(size(current_area,1)); % calculate standard error across mice
    
    x = [1:length(mean_elapse_repeat)]';
    y = mean_elapse_repeat';
    dy = ste_elapse_repeat';
    
    hold on
    fill([x;flipud(x)],[y-dy;flipud(y+dy)],new_colors(area,:),'linestyle','none','facealpha',0.4);
    plt(area) = plot(y,'color',new_colors2(area,:),'linewidth',3);
end
text(0.05,0.075,'1 repeat = 30 seconds','Units','normalized','FontSize',11)
set(gca,'xtick',[1,10,20,29],'ytick',[-0.4:0.1:0])
xlabel('Elapsed time (# of movie repeats)')
ylabel('Tuning curve similarity index')
title('Neuropixels - seconds to minutes')

% statistical analysis for differences across pairs of visual areas
V1_tuning_corr = tuning_similarity_index{1}; % subset tuning curve similarity indices of area V1
V2_tuning_corr = tuning_similarity_index{2}; % subset tuning curve similarity indices of area LM

LGN_tuning_corr = tuning_similarity_index{3}; % subset tuning curve similarity indices of area dLGN
LP_tuning_corr = tuning_similarity_index{4}; % subset tuning curve similarity indices of area LP

% perform multiple Mann-Whitney rank-sum tests btween pairs of visual
% areas and visualize the results on the existing figure
for trial = 1:29 % loop over time intervals
    p = ranksum(V1_tuning_corr(:,trial),V2_tuning_corr(:,trial)); % perform two-sided Mann-Whitney rank-sum between area V1 and area LM
    if p < 0.05
        hold on
        plot(trial,0,'*','color',new_colors(1,:),'markersize',8)
    end
    
    p = ranksum(LGN_tuning_corr(:,trial),LP_tuning_corr(:,trial)); % perform two-sided Mann-Whitney rank-sum between area dLGN and area LP
    if p < 0.05
        hold on
        plot(trial,-0.4,'*','color',new_colors(3,:),'markersize',8)
    end
end
lgd = legend(plt,brain_areas(area_list));
legend('boxoff')
lgd.Position = [0.2 0.225 0.1 0.1];

%% Figure 5C - Visual hierarchy - ensemble rate correlation between movie repeats - calcium imaging

nat_movie = 1; % natural movie 1
cell_cutoff = 20; % cell count threshold of 20 cells

elapse_repeat_rate_corr_area = {}; % define an empty variable for ensemble rate corr across areas and mice
for area = 1:2 % loop over visual areas
    % valid_mice - true only if the number of cells recorded is
    % atleast 20 in each of the recording sessions
    valid_mice = min(calcium_excitatory_cell_count{area},[],2) >= cell_cutoff;
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie); % subset all mice that met the requirements of 'valid mice'
    
    elapse_repeat_rate_corr = [];  % define an empty variable for ensemble rate corr as a function of elapsed time across mice (size of #mice by 9)
    for mouse = 1:length(current_area) % loop over mice
        clc;
        disp(['Calculating ensemble rate correlation between repeats for calcium imaging data:'])
        disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse}; % subset the neuronal activity of a single mouse (size of #cells by 30 time bins by 30 movie repeats)
        
        current_mouse_sess1 =  current_mouse(:,:,1:10); % subset neuronal activity during the first session (repeats 1-10)
        valid_cells = mean(squeeze(mean(current_mouse_sess1,2,'omitnan')),2,'omitnan')>0; % include only cells that were active in at least one movie repeat
        valid_current_mouse_sess1 = squeeze(mean(current_mouse_sess1(valid_cells,:,:),2,'omitnan')); % calculate the average activity over time bins only for the cells that passed the requirements of 'valid_cells'
        
        current_mouse_sess2 =  current_mouse(:,:,11:20); % subset neuronal activity during the second session (repeats 11-20)
        valid_cells = mean(squeeze(mean(current_mouse_sess2,2,'omitnan')),2,'omitnan')>0; % include only cells that were active in at least one movie repeat
        valid_current_mouse_sess2 = squeeze(mean(current_mouse_sess2(valid_cells,:,:),2,'omitnan')); % calculate the average activity over time bins only for the cells that passed the requirements of 'valid_cells'
        
        current_mouse_sess3 =  current_mouse(:,:,21:30); % subset neuronal activity during the third session (repeats 21-30)
        valid_cells = mean(squeeze(mean(current_mouse_sess3,2,'omitnan')),2,'omitnan')>0; % include only cells that were active in at least one movie repeat
        valid_current_mouse_sess3 = squeeze(mean(current_mouse_sess3(valid_cells,:,:),2,'omitnan')); % calculate the average activity over time bins only for the cells that passed the requirements of 'valid_cells'
        
        rate_corr_across_session = []; % define an empty variable for ensemble rate corr across movie repeats
        for repeat1 = 1:10 % loop over movie repeats
            for repeat2 = 1:10 % loop over movie repeats
                
                rate_corr_across_session(repeat1,repeat2,1) = corr(valid_current_mouse_sess1(:,repeat1),valid_current_mouse_sess1(:,repeat2)); % calculate the ensemble rate corr between pair of movie repeats in session 1
                rate_corr_across_session(repeat1,repeat2,2) = corr(valid_current_mouse_sess2(:,repeat1),valid_current_mouse_sess2(:,repeat2)); % calculate the ensemble rate corr between pair of movie repeats in session 2
                rate_corr_across_session(repeat1,repeat2,3) = corr(valid_current_mouse_sess3(:,repeat1),valid_current_mouse_sess3(:,repeat2)); % calculate the ensemble rate corr between pair of movie repeats in session 3
                
            end
        end
        rate_corr_across_session(rate_corr_across_session<0) = 0; % rectify negative ensemble rate corr values to zero
        mean_rate_corr_across_mice = mean(rate_corr_across_session,3,'omitnan'); % average ensemble rate corr across sessions
        
        % calculate ensemble rate correlation as function of elapsed time
        for diagonal = 1:9 % loop over diagonals
            elapse_repeat_rate_corr(mouse,diagonal) = mean(diag(mean_rate_corr_across_mice,diagonal),'omitnan'); % average ensemble rate corr across diagonals
        end
        
    end
    elapse_repeat_rate_corr_area{area} = elapse_repeat_rate_corr; % store mean ensemble rate corr values across mice for each area
end

% define color schemes for visualization
new_colors = [0 0.7 0.8;0.7 0.7 0.7];
new_colors2 = [0 0.5 0.6;0.5 0.5 0.5];
plt = []; % define empty variable that will store plot information for each area
rate_similarity_index = {}; % define empty variable that will store ensemble rate similarity index values of indevidual mice across areas
figure('units','normalized','position',[0.3 0.3 0.25 0.4]) % visualization of ensemble rate similarity index across areas
for area = 1:2 % loop over areas
    current_area = elapse_repeat_rate_corr_area{area}; % subset ensemble rate corr values of specific visual area
    rate_similarity_index{area} = [current_area-current_area(:,1)]./[current_area+current_area(:,1)]; % calculate for each mouse its ensemble rate similarity index
    
    mean_elapse_repeat = mean(rate_similarity_index{area},'omitnan'); % calculate mean ensemble rate similarity index across mice
    std_elapse_repeat = std(rate_similarity_index{area},'omitnan'); % calculate standard deviation across mice
    ste_elapse_repeat = std_elapse_repeat./sqrt(size(current_area,1)); % calculate standard error across mice
    
    x = [1:length(mean_elapse_repeat)]';
    y = mean_elapse_repeat';
    dy = ste_elapse_repeat';
    
    hold on
    fill([x;flipud(x)],[y-dy;flipud(y+dy)],new_colors(area,:),'linestyle','none','facealpha',0.4);
    plt(area) = plot(y,'color',new_colors2(area,:),'linewidth',3);
end

text(0.05,0.075,'1 repeat = 30 seconds','Units','normalized','FontSize',11)
set(gca,'xtick',[1:9])
xlabel('Elapsed time (# of movie repeats)')
ylabel('Ensemble rate similarity index')
title('Calcium imaging - seconds to minutes')

% statistical analysis for differences across pairs of visual areas
V1_rate_corr = rate_similarity_index{1}; % subset ensemble rate similarity indices of area V1
V2_rate_corr = rate_similarity_index{2}; % subset ensemble rate similarity indices of area LM

% perform multiple Mann-Whitney rank-sum tests btween pairs of visual
% areas and visualize the results on the existing figure
for trial = 1:9 % loop over time intervals
    p = ranksum(V1_rate_corr(:,trial),V2_rate_corr(:,trial)); % perform two-sided Mann-Whitney rank-sum between area V1 and area LM
    if p < 0.05
        hold on
        plot(trial,-0.01,'*','color',[0.4 0.4 0.4],'markersize',8)
    end
    
end
lgd = legend(plt,brain_areas(1:2));
legend('boxoff')
lgd.Position = [0.175 0.19 0.1 0.1];

%% Figure 5D - Visual hierarchy - tuning curve correlation between movie repeats - calcium imaging

nat_movie = 1; % natural movie 1
cell_cutoff = 20; % cell count threshold of 20 cells

elapse_repeat_tuning_corr_area = {}; % define an empty variable for tuning curve corr across areas and mice
for area = 1:2% loop over visual areas
    % valid_mice - true only if the number of cells recorded is
    % atleast 20 in each of the recording sessions
    valid_mice = min(calcium_excitatory_cell_count{area},[],2) >= cell_cutoff;
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie); % subset all mice that met the requirements of 'valid mice'
    
    elapse_repeat_tuning_corr = []; % define an empty variable for ensemble rate corr as a function of elapsed time across mice (size of #mice by 9)
    for mouse = 1:length(current_area) % loop over mice
        clc;
        disp(['Calculating tuning curve correlation between repeats for calcium imaging data:'])
        disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse}; % subset the neuronal activity of a single mouse (size of #cells by 30 time bins by 30 movie repeats)
        
        current_mouse_sess1 =  current_mouse(:,:,1:10); % subset neuronal activity during the first session (repeats 1-10)
        valid_cells = mean(squeeze(mean(current_mouse_sess1,2,'omitnan')),2,'omitnan')>0; % include only cells that were active in at least one movie repeat
        valid_current_mouse_sess1 = current_mouse_sess1(valid_cells,:,:); % subset only the cells that passed the requirements of 'valid_cells'
        
        current_mouse_sess2 =  current_mouse(:,:,11:20); % subset neuronal activity during the second session (repeats 11-20)
        valid_cells = mean(squeeze(mean(current_mouse_sess2,2,'omitnan')),2,'omitnan')>0; % include only cells that were active in at least one movie repeat
        valid_current_mouse_sess2 = current_mouse_sess2(valid_cells,:,:); % subset only the cells that passed the requirements of 'valid_cells'
        
        current_mouse_sess3 =  current_mouse(:,:,21:30); % subset neuronal activity during the third session (repeats 21-30)
        valid_cells = mean(squeeze(mean(current_mouse_sess3,2,'omitnan')),2,'omitnan')>0; % include only cells that were active in at least one movie repeat
        valid_current_mouse_sess3 = current_mouse_sess3(valid_cells,:,:); % subset only the cells that passed the requirements of 'valid_cells'
        
        tuning_corr_across_session = [];% define an empty variable for mean tuning curve across movie repeats
        for repeat1 = 1:10 % loop over movie repeats
            for repeat2 = 1:10 % loop over movie repeats
                tuning_corr = corr(valid_current_mouse_sess1(:,:,repeat1)',valid_current_mouse_sess1(:,:,repeat2)'); % tuning curve corr between pair of movie repeats in session 1
                tuning_corr_across_session(repeat1,repeat2,1) = mean(diag(tuning_corr),'omitnan'); % calculate the average tuning curve corr across values of corresponding cells
                
                tuning_corr = corr(valid_current_mouse_sess2(:,:,repeat1)',valid_current_mouse_sess2(:,:,repeat2)'); % tuning curve corr between pair of movie repeats in session 2
                tuning_corr_across_session(repeat1,repeat2,2) = mean(diag(tuning_corr),'omitnan'); % calculate the average tuning curve corr across values of corresponding cells
                
                tuning_corr = corr(valid_current_mouse_sess3(:,:,repeat1)',valid_current_mouse_sess3(:,:,repeat2)'); % tuning curve corr between pair of movie repeats in session 3
                tuning_corr_across_session(repeat1,repeat2,3) = mean(diag(tuning_corr),'omitnan'); % calculate the average tuning curve corr across values of corresponding cells
                
            end
        end
        tuning_corr_across_session(tuning_corr_across_session<0) = 0; % rectify negative tuning curve corr values to zero
        mean_tuning_corr_across_mice = mean(tuning_corr_across_session,3,'omitnan'); % average tuning curve corr across sessions
        
        % calculate tuning curve correlation as function of elapsed time
        for diagonal = 1:9 % loop over diagonals
            elapse_repeat_tuning_corr(mouse,diagonal) = mean(diag(mean_tuning_corr_across_mice,diagonal),'omitnan'); % average tuning curve corr across diagonals
        end

    end
    elapse_repeat_tuning_corr_area{area} = elapse_repeat_tuning_corr; % store mean tuning curve corr values across mice for each area
end

% define color schemes for visualization
new_colors = [0 0.7 0.8;0.7 0.7 0.7];
new_colors2 = [0 0.5 0.6;0.5 0.5 0.5];
plt = []; % define empty variable that will store plot information for each area
tuning_similarity_index = {}; % define empty variable that will store tuning curve similarity index values of indevidual mice across areas
figure('units','normalized','position',[0.3 0.3 0.25 0.4]) % visualization oftuning curve similarity index across areas
for area = 1:2
    current_area = elapse_repeat_tuning_corr_area{area}; % subset tuning curve corr values of specific visual area
    tuning_similarity_index{area} = [current_area-current_area(:,1)]./[current_area+current_area(:,1)]; % calculate for each mouse its tuning curve similarity index
    
    mean_elapse_repeat = mean(tuning_similarity_index{area},'omitnan'); % calculate mean tuning curve similarity index across mice
    std_elapse_repeat = std(tuning_similarity_index{area},'omitnan'); % calculate standard deviation across mice
    ste_elapse_repeat = std_elapse_repeat./sqrt(size(current_area,1));  % calculate standard error across mice
    
    x = [1:length(mean_elapse_repeat)]';
    y = mean_elapse_repeat';
    dy = ste_elapse_repeat';
    
    hold on
    fill([x;flipud(x)],[y-dy;flipud(y+dy)],new_colors(area,:),'linestyle','none','facealpha',0.4);
    plt(area) = plot(y,'color',new_colors2(area,:),'linewidth',3);
end

text(0.05,0.075,'1 repeat = 30 seconds','Units','normalized','FontSize',11)
set(gca,'xtick',[1:9],'ytick',[-0.08:0.02:0])
xlim([0.5 9.5])
xlabel('Elapsed time (# of movie repeats)')
ylabel('Tuning curve similarity index')
title('Calcium imaging - seconds to minutes')

% statistical analysis for differences across pairs of visual areas
V1_tuning_corr = tuning_similarity_index{1}; % subset tuning curve similarity indices of area V1
V2_tuning_corr = tuning_similarity_index{2}; % subset tuning curve similarity indices of area LM

% perform multiple Mann-Whitney rank-sum tests btween pairs of visual
% areas and visualize the results on the existing figure
for trial = 1:9 % loop over time intervals
    p = ranksum(V1_tuning_corr(:,trial),V2_tuning_corr(:,trial)); % perform two-sided Mann-Whitney rank-sum between area V1 and area LM
    if p < 0.05
        hold on
        plot(trial,-0.005,'*','color',[0.4 0.4 0.4],'markersize',8)
    end
    
end
lgd = legend(plt,brain_areas(1:2));
legend('boxoff')
lgd.Position = [0.175 0.19 0.1 0.1];


%% Figure 5E - Visual hierarchy - Ensemble rate correlation between blocks - calcium imaging

nat_movie = 2; % natural movie 3
cell_cutoff = 20; % cell count threshold of 20 cells

within_between_stability_area = {}; % define empty variable that will store all ensemble rate corr values across blocks across mice and visual areas
for area = 1:2 % loop over visual areas
    % valid_mice - true only if the number of cells recorded is
            % atleast 20 in session A
    valid_mice = calcium_excitatory_cell_count{area}(:,1) >= cell_cutoff;
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie); % subset all mice that met the requirements of 'valid mice'
    
    within_between_stability = []; % define empty variable that will store ensemble rate corr values for all mice
    for mouse = 1:length(current_area) % loop over mice
        clc;
        disp(['Calculating ensemble rate correlation between blocks for calcium imaging data:'])
        disp(['Stimulus: Natural movie 3 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse}; % subset a single mouse (size of #cells by 30 time bins by 10 movie repeats)
        
        mean_activity = [];
        mean_activity(:,1) = mean(mean(current_mouse(:,:,1:2),2,'omitnan'),3,'omitnan'); % calculate the mean activity across repreats and time bins for the first half of block A (repeats 1-2)
        mean_activity(:,2) = mean(mean(current_mouse(:,:,3:5),2,'omitnan'),3,'omitnan'); % calculate the mean activity across repreats and time bins for the second half of block A (repeats 3-5)
        mean_activity(:,3) = mean(mean(current_mouse(:,:,6:7),2,'omitnan'),3,'omitnan'); % calculate the mean activity across repreats and time bins for the first half of block B (repeats 6-7)
        mean_activity(:,4) = mean(mean(current_mouse(:,:,8:10),2,'omitnan'),3,'omitnan'); % calculate the mean activity across repreats and time bins for the second half of block B (repeats 8-10)
        
        rate_corr = corr(mean_activity); % calculate the ensemble rate correlation between the four halves of the two blocks
        rate_corr(rate_corr<0) = 0; % rectify negative ensemble rate corr values to zero
        
        % calculate and store the average ensemble rate corr between halves of the
            % same block ('within-block') and average ensemble rate corr between different
            % halves of different blocks ('between-blocks')
            within_block = [rate_corr(1,2),rate_corr(3,4)]; % subset ensemble rate corr between halves of the same block
            between_blocks = rate_corr(1:2,3:4); % subset ensemble rate corr between halves of different blocks
            
            within_between_stability(mouse,1) = mean(within_block(:),'omitnan'); % calculate the average ensemble rate corr within a block
            within_between_stability(mouse,2) = mean(between_blocks(:),'omitnan'); % calculate the average ensemble rate corr between blocks

    end
    within_between_stability_area{area} = within_between_stability; % store mean ensemble rate corr across mice for each area
end

rate_similarity_index = []; % define empty variable that will store ensemble rate similarity index values of indevidual mice across areas
plt = []; % define empty variable that will store plot information for each area
figure('units','normalized','position',[0.3 0.3 0.2 0.3]) % visualization of ensemble rate corr across areas
for area = [2,1] % loop over areas
    current_area = within_between_stability_area{area}; % subset ensemble rate corr values of specific visual area
    
    rate_similarity_index{area} = [current_area(:,2)-current_area(:,1)]./[current_area(:,2) + current_area(:,1)]; % calculate for each mouse its ensemble rate similarity index
    
    mean_stability = mean(current_area,'omitnan'); % calculate mean ensemble rate corr across mice
    std_stability = std(current_area,'omitnan'); % calculate standard deviation across mice
    ste_stability = std_stability./sqrt(size(current_area,1)); % calculate standard error across mice
    
    hold on
    if area == 1
        plt(1) = errorbar(mean_stability,ste_stability,'-o','color',colors(area,:),'linewidth',3,'capsize',0);
    elseif area == 2
        plt(2) = errorbar(mean_stability,ste_stability,'-o','color',[0.6 0.6 0.6],'linewidth',3,'capsize',0);
    end
end
xlim([0.5 2.5])
set(gca,'xtick',1:2,'xticklabel',{'Within block','Between blocks'})
ylabel('Ensemble rate correlation')
legend(plt,brain_areas(1:2),'location','southwest')
legend('boxoff')


figure('units','normalized','position',[0.5 0.4 0.1 0.2]) % visualize distribution of ensemble rate similarity index across areas
xlim([0 3])
hold on
plot(xlim,[0 0 ],'--','color',[0.2 0.2 0.2])
stability_diff_areas = nan(size(rate_similarity_index{1},1),2);
stability_diff_areas(:,1) = rate_similarity_index{1};
stability_diff_areas(1:size(rate_similarity_index{2},1),2) = rate_similarity_index{2};
figure_boxplot(stability_diff_areas)
ylabel('Similarity index')
set(gca,'xtick',1:2,'xticklabel',{'V1','LM'})
ylim([-0.225 0.125])

% statistical analysis for differences across pairs of visual areas
[pval,~,stat] = ranksum(rate_similarity_index{1},rate_similarity_index{2}); % perform two-sided Mann-Whitney rank-sum between area V1 and area LM
zval = stat.zval;

clc;
disp(['Ensemble rate similarity index of V1 compared to LM'])
disp(['Two-tailed Mann-Whitney rank-sum test: Z=',num2str(zval),' p=',num2str(pval)])


%% Figure 5F - Visual hierarchy - Tuning curve correlation between blocks - calcium imaging

nat_movie = 2; % natural movie 3
cell_cutoff = 20; % cell count threshold of 20 cells

within_between_stability_area = {}; % define empty variable that will store all tuning curve corr values across blocks across mice and visual areas
for area = 1:2 % loop over visual areas
        % valid_mice - true only if the number of cells recorded is
            % atleast 20 in session A
    valid_mice = calcium_excitatory_cell_count{area}(:,1) >= cell_cutoff;
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie); % subset all mice that met the requirements of 'valid mice'
    
    within_between_stability = []; % define empty variable that will store ensemble rate corr values for all mice
    for mouse = 1:length(current_area) % loop over mice
        clc;
        disp(['Calculating tuning curve correlation between blocks for calcium imaging data:'])
        disp(['Stimulus: Natural movie 3 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse}; % subset a single mouse (size of #cells by 30 time bins by 10 movie repeats)
        
        tuning_corr_across_mice = [];
        tuning_corr_across_mice(:,:,1) = mean(current_mouse(:,:,1:2),3,'omitnan'); % calculate the mean activity across repreats for the first half of block A (repeats 1-2)
        tuning_corr_across_mice(:,:,2) = mean(current_mouse(:,:,3:5),3,'omitnan'); % calculate the mean activity across repreats for the second half of block A (repeats 3-5)
        tuning_corr_across_mice(:,:,3) = mean(current_mouse(:,:,6:7),3,'omitnan'); % calculate the mean activity across repreats for the first half of block B (repeats 6-7)
        tuning_corr_across_mice(:,:,4) = mean(current_mouse(:,:,8:10),3,'omitnan'); % calculate the mean activity across repreats for the second half of block B (repeats 8-10)
        
        % calculate the tuning curve correlation between the four halves of the two blocks
        mean_tuning_corr_across_mice = []; % define an empty variable to store the tuning curve correlation across block halves
        for half1 = 1:4 % loop over blocks halves
            for half2 = 1:4 % loop over blocks halves
                tuning_corr =  corr(tuning_corr_across_mice(:,:,half1)',tuning_corr_across_mice(:,:,half2)');% calculate the tuning curve correlation between block halves
                mean_tuning_corr_across_mice(half1,half2) =  median(diag(tuning_corr),'omitnan'); % calculate the median tuning curve corr across corresponding cells
            end
        end
        mean_tuning_corr_across_mice(mean_tuning_corr_across_mice<0) = 0; % rectify negative tuning curve corr values to zero
        
        % calculate and store the average tuning curve corr between halves of the
            % same block ('within-block') and average tuning curve corr between different
            % halves of different blocks ('between-blocks')
            within_block = [mean_tuning_corr_across_mice(1,2), mean_tuning_corr_across_mice(3,4)]; % subset tuning curve corr between halves of the same block
            between_blocks = mean_tuning_corr_across_mice(1:2,3:4); % subset tuning curve corr between halves of different blocks
            
            within_between_stability(mouse,1) = mean(within_block(:),'omitnan'); % calculate the average tuning curve corr within a block
            within_between_stability(mouse,2) = mean(between_blocks(:),'omitnan'); % calculate the average tuning curve corr between blocks

         end
    within_between_stability_area{area} = within_between_stability; % store mean tuning curve corr across mice for each area
end

tuning_similarity_index = [];  % define empty variable that will store tuning curve similarity index values of indevidual mice across areas
plt = []; % define empty variable that will store plot information for each area
figure('units','normalized','position',[0.3 0.3 0.2 0.3]) % visualization of tuning curve corr across areas
for area = [2,1]
    current_area = within_between_stability_area{area}; % subset tuning curve corr values of specific visual area
    
    tuning_similarity_index{area} = [current_area(:,2)-current_area(:,1)]./[current_area(:,2) + current_area(:,1)]; % calculate for each mouse its tuning curve similarity index
    
    mean_stability = mean(current_area,'omitnan'); % calculate mean tuning curve corr across mice
    std_stability = std(current_area,'omitnan'); % calculate standard deviation across mice
    ste_stability = std_stability./sqrt(size(current_area,1)); % calculate standard error across mice
    
    hold on
    if area == 1
        plt(1) = errorbar(mean_stability,ste_stability,'-o','color',colors(area,:),'linewidth',3,'capsize',0);
    elseif area == 2
        plt(2) = errorbar(mean_stability,ste_stability,'-o','color',[0.6 0.6 0.6],'linewidth',3,'capsize',0);
    end
end
xlim([0.5 2.5])
ylim([0.3375 0.43])
set(gca,'xtick',1:2,'xticklabel',{'Within block','Between blocks'})
ylabel('Tuning curve correlation')
legend(plt,brain_areas(1:2),'location','southwest')
legend('boxoff')


figure('units','normalized','position',[0.5 0.4 0.1 0.2]) % visualize distribution of tuning curve similarity index across areas
xlim([0 3])
hold on
plot(xlim,[0 0 ],'--','color',[0.2 0.2 0.2])
stability_diff_areas = nan(size(tuning_similarity_index{1},1),2);
stability_diff_areas(:,1) = tuning_similarity_index{1};
stability_diff_areas(1:size(tuning_similarity_index{2},1),2) = tuning_similarity_index{2};
figure_boxplot(stability_diff_areas)
ylabel('Similarity index')
set(gca,'xtick',1:2,'xticklabel',{'V1','LM'})
ylim([-0.175 0.125])

% statistical analysis for differences across pairs of visual areas
[pval,~,stat] = ranksum(tuning_similarity_index{1},tuning_similarity_index{2}); % perform two-sided Mann-Whitney rank-sum between area V1 and area
zval = stat.zval;

clc;
disp(['Tuning curve similarity index of V1 compared to LM'])
disp(['Two-tailed Mann-Whitney rank-sum test: Z=',num2str(zval),' p=',num2str(pval)])


%% Figure 5G - Visual hierarchy - Ensemble rate correlation between sessions - calcium imaging

nat_movie = 1; % natural movie 1
cell_cutoff = 20; % cell count threshold of 20 cells

within_between_session_stability_area = {}; % define empty variable that will store the ensemble rate corr values of within and between sessions across mice and visual areas
for area = 1:2 % loop over visual areas
    % valid_mice - true only if the number of cells recorded is atleast 20 in each of the recording sessions
    valid_mice = min(calcium_excitatory_cell_count{area},[],2)>= cell_cutoff;
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie); % mice that passed the requirments of 'valid_mice'
    
    within_between_session_stability = []; % define empty variable that will store the ensemble rate correlation values for all mice
    for mouse = 1:length(current_area) % loop over visual areas
        clc;
        disp(['Calculating ensemble rate correlation between sessions for calcium imaging data:'])
        disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse}; % subset neuronal activity for a single mouse (#cells by 30 time bins by 30 movie repeats)
        
       %calculate the average activity rate for session halves across movie repeats
            mean_activity_per_half = []; % define an empty variable that will store the averaged neuronal responses of each session half
            for half = 1:6 % loop over session halves
                current_half = [1:5] + 5*(half-1); % define range of movie repeats for current half
                mean_activity_per_half(:,:,half) = mean(current_mouse(:,:,current_half),3,'omitnan');  % average the neuronal activity across chosen movie repeats
            end
        
            % calculate ensemble rate correlation between sessions halves using only
            % cells that were active in both compared time points (either within or between sessions)
        rate_corr = []; % define empty variable that will store ensemble rate corr values for all mice
        for halfA = 1:size(mean_activity_per_half,3) % loop over session halves
            halfA_activity = mean_activity_per_half(:,:,halfA); % subset neuronal activity for a single half
            
            for halfB = 1:size(mean_activity_per_half,3) % loop over session halves
                halfB_activity = mean_activity_per_half(:,:,halfB); % subset neuronal activity for a single half
                
                % valid_cells - cells that were active in both of compared time points
                valid_cells = [mean(halfA_activity,2,'omitnan') ~= 0] & [mean(halfB_activity,2,'omitnan') ~= 0]; % true if average activity is above zero in both session halves
                valid_halfA_activity = mean(halfA_activity(valid_cells,:),2,'omitnan'); % calculate average activity across time bins only for the cells that met the requirments of 'valid_cells'
                valid_halfB_activity = mean(halfB_activity(valid_cells,:),2,'omitnan'); % calculate average activity across time bins only for the cells that met the requirments of 'valid_cells'
                
                rate_corr(halfA,halfB) = corr(valid_halfA_activity,valid_halfB_activity); % calculate ensemble rate corr between sessions halves
                
            end
        end
        rate_corr(rate_corr<0) = 0; % rectify negative ensemble rate corr values to zero
        
        within_session_values = [rate_corr(1,2),rate_corr(3,4),rate_corr(5,6)]; % ensemble rate corr values between halves of the same session (within session)
            proximal_sessions_values = [rate_corr(1:2,3:4),rate_corr(3:4,5:6)]; % ensemble rate corr values between halves of proximal sessions (sessions 1&2 and sessions 2&3)
            distal_session_values = rate_corr(1:2,5:6); % ensemble rate corr values between halves of distal sessions (sessions 1&3)
            
            within_between_session_stability(mouse,1) = mean(within_session_values(:),'omitnan'); % average ensemble rate corr values for within session
            within_between_session_stability(mouse,2) = mean(proximal_sessions_values(:),'omitnan'); % average ensemble rate corr values for proximal sessions
            within_between_session_stability(mouse,3) = mean(distal_session_values(:),'omitnan'); % average ensemble rate corr values for distal sessionsend

    end
    within_between_session_stability_area{area} = within_between_session_stability; % store ensemble rate corr values for all mice of a given area 
end

rate_similarity_index = {}; % define empty variable that will store ensemble rate similarity index values of indevidual mice across areas
figure('units','normalized','position',[0.3 0.3 0.2 0.3]) % visualization of ensemble rate corr across areas
for area = [2,1] % loop over areas
    current_area = within_between_session_stability_area{area}; % subset ensemble rate corr values of specific visual area
    
    rate_similarity_index{area} = [current_area-current_area(:,1)]./[current_area+current_area(:,1)]; % calculate for each mouse its ensemble rate similarity index
    
    mean_stability = mean(rate_similarity_index{area},'omitnan'); % calculate mean ensemble rate corr across mice
    std_stability = std(rate_similarity_index{area},'omitnan'); % calculate standard deviation across mice
    ste_stability = std_stability ./sqrt(size(rate_similarity_index{area},1)); % calculate standard error across mice
    
    hold on
    if area == 2
        errorbar(mean_stability,ste_stability,'o','color',[0.7 0.7 0.7],...
            'markerfacecolor',[0.7 0.7 0.7],'capsize',0,'linestyle','-','linewidth',2)
    elseif area == 1
        errorbar(mean_stability,ste_stability,'o','color',colors(area,:),...
            'markerfacecolor',colors(area,:),'capsize',0,'linestyle','-','linewidth',2)
    end
end
xlim([0.5 3.5])
ylim([-0.225 0])
set(gca,'xtick',1:3,'xticklabel',{'Within session','Proximal sessions','Distal sessions'})
xtickangle(15)
ylabel('Ensemble rate similarity index')
title('Calcium imaging - minutes to days:')

% perform multiple Mann-Whitney rank-sum tests btween pairs of visual
% areas and visualize the results on the existing figure
for sess = 1:3 % loop over time intervals
    p = ranksum(rate_similarity_index{1}(:,sess),rate_similarity_index{2}(:,sess)); % perform two-sided Mann-Whitney rank-sum between area V1 and area LM
    if p < 0.05
        hold on
        plot(sess,0,'*','color',[0.4 0.4 0.4],'markersize',8)
    end
end

%% Figure 5H - Visual hierarchy - Tuning curve correlation between sessions - calcium imaging

nat_movie = 1; % natural movie 1
cell_cutoff = 20; % cell count threshold of 20 cells

within_between_session_stability_area = {}; % define empty variable that will store the tuning curve corr values of within and between sessions across mice and visual areas
for area = 1:2 % loop over areas
    % valid_mice - true only if the number of cells recorded is atleast 20 in each of the recording sessions
    valid_mice = min(calcium_excitatory_cell_count{area},[],2)>= cell_cutoff;
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie); % mice that passed the requirments of 'valid_mice'
    
    within_between_session_stability = []; % define empty variable that will store the ensemble rate correlation values for all mice
    for mouse = 1:length(current_area) % loop over visual areas
        clc;
        disp(['Calculating tuning curve correlation between sessions for calcium imaging data:'])
        disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse}; % subset neuronal activity for a single mouse (#cells by 30 time bins by 30 movie repeats)
        
        %calculate the average activity rate for session halves across movie repeats
            mean_activity_per_half = []; % define an empty variable that will store the averaged neuronal responses of each session half
            for half = 1:6 % loop over session halves
                current_half = [1:5] + 5*(half-1); % define range of movie repeats for current half
                mean_activity_per_half(:,:,half) = mean(current_mouse(:,:,current_half),3,'omitnan');  % average the neuronal activity across chosen movie repeats
            end
        
            % calculate tuning curve correlation between sessions halves using only
            % cells that were active in both compared time points (either within or between sessions)
        median_tuning_corr = []; % define empty variable that will store tuning curve corr values for all mice
        for halfA = 1:size(mean_activity_per_half,3) % loop over session halves
            halfA_activity = mean_activity_per_half(:,:,halfA); % subset neuronal activity for a single half
            
            for halfB = 1:size(mean_activity_per_half,3) % loop over session halves
                halfB_activity = mean_activity_per_half(:,:,halfB); % subset neuronal activity for a single half
                
                % valid_cells - cells that were active in both of compared time points
                valid_cells = [nanmean(halfA_activity,2) ~= 0] & [nanmean(halfB_activity,2) ~= 0]; % true if average activity is above zero in both session halves
                valid_halfA_activity = halfA_activity(valid_cells,:); % subset only the cells that met the requirments of 'valid_cells'
                valid_halfB_activity = halfB_activity(valid_cells,:); % subset only the cells that met the requirments of 'valid_cells'
                
                tuning_corr = corr(valid_halfA_activity',valid_halfB_activity'); % calculate tuning curve corr between sessions halves
                median_tuning_corr(halfA,halfB) = median(diag(tuning_corr),'omitnan'); % calculate the median tuning curve corr across corresponding cells
            end
        end
        median_tuning_corr(median_tuning_corr<0) = 0; % rectify negative tuning curve corr values to zero
        
        within_session_values = [median_tuning_corr(1,2),median_tuning_corr(3,4),median_tuning_corr(5,6)]; % tuning curve corr values between halves of the same session (within session)
            proximal_sessions_values = [median_tuning_corr(1:2,3:4),median_tuning_corr(3:4,5:6)]; % tuning curve corr values between halves of proximal sessions (sessions 1&2 and sessions 2&3)
            distal_session_values = median_tuning_corr(1:2,5:6); % tuning curve corr values between halves of distal sessions (sessions 1&3)
            
            within_between_session_stability(mouse,1) = mean(within_session_values(:),'omitnan'); % average tuning curve corr values for within session
            within_between_session_stability(mouse,2) = mean(proximal_sessions_values(:),'omitnan'); % average tuning curve corr values for proximal sessions
            within_between_session_stability(mouse,3) = mean(distal_session_values(:),'omitnan'); % average tuning curve corr values for distal sessionsend

    end
    within_between_session_stability_area{area} = within_between_session_stability; % store tuning curve corr values for all mice of a given area 
end

tuning_similarity_index = {}; % define empty variable that will store tuning curve similarity index values of indevidual mice across areas
figure('units','normalized','position',[0.3 0.3 0.2 0.3]) % visualization of ensemble rate corr across areas
for area = [2,1] % loop over areas
    current_area = within_between_session_stability_area{area}; % subset ensemble rate corr values of specific visual area
    
    tuning_similarity_index{area} = [current_area-current_area(:,1)]./[current_area+current_area(:,1)]; % calculate for each mouse its tuning curve similarity index
    
    mean_stability = mean(tuning_similarity_index{area},'omitnan'); % calculate mean tuning curve corr across mice
    std_stability = std(tuning_similarity_index{area},'omitnan'); % calculate standard deviation across mice
    ste_stability = std_stability ./sqrt(size(tuning_similarity_index{area},1)); % calculate standard error across mice
    
    hold on
    if area == 2
        errorbar(mean_stability,ste_stability,'o','color',[0.7 0.7 0.7],...
            'markerfacecolor',[0.7 0.7 0.7],'capsize',0,'linestyle','-','linewidth',2)
    elseif area == 1
        errorbar(mean_stability,ste_stability,'o','color',colors(area,:),...
            'markerfacecolor',colors(area,:),'capsize',0,'linestyle','-','linewidth',2)
    end
end
xlim([0.5 3.5])
ylim([-0.225 0])
set(gca,'xtick',1:3,'xticklabel',{'Within session','Proximal sessions','Distal sessions'})
xtickangle(15)
ylabel('Tuning curve similarity index')
title('Calcium imaging - minutes to days:')

% perform multiple Mann-Whitney rank-sum tests btween pairs of visual
% areas and visualize the results on the existing figure
for sess = 1:3 % loop over time intervals
    p = ranksum(tuning_similarity_index{1}(:,sess),tuning_similarity_index{2}(:,sess)); % perform two-sided Mann-Whitney rank-sum between area V1 and area LM
    if p < 0.05
        hold on
        plot(sess,-0.1,'*','color',[0.4 0.4 0.4],'markersize',8)
    end
end

%% Figure 6A - Internal structure in reduced space - single animal v1 example

area = 1; % area V1
mouse = 91; % example mouse #91
current_area = calcium_excitatory_population_vectors{area}{mouse,3}*30; % subset neuronal activity of a single example mouse
valid_cells = mean(mean(current_area(:,:,1:10),3,'omitnan'),2,'omitnan')>0; % find the cells that were active during session 1
valid_cells_sess1 = current_area(valid_cells,:,1:10); % subset only the cells the passed the requirments of 'valid_cells'

% sort cells based on their time bin with peak activiy in movie repeat #4
[~,i] = max(valid_cells_sess1(:,:,4),[],2);
[~,sorted_ind1] = sort(i); % store sorted indices for repeat 4

% sort cells based on their time bin with peak activiy in movie repeat #8
[~,i] = max(valid_cells_sess1(:,:,8),[],2);
[~,sorted_ind3] = sort(i); % store sorted indices for repeat 8

repeat_list =[4,6,8]; %define a list with example movie repeats

figure('units','normalized','position',[0.3 0.3 0.275 0.4]) % visualize sorted neuronal activity during repeats 4, 6 and 8
for trial = 1:3 % loop over movie repeats
    
    subplot(2,3,trial) % sorted based on repeat #4
    current_trial = valid_cells_sess1(:,:,repeat_list(trial));
    
    current_trial_sorted = current_trial(sorted_ind1,:);
    imagesc(current_trial_sorted,[0 30])
    title(['Repeat #',num2str(repeat_list(trial))])
    if trial == 1
        ylabel('Sorted by repeat 4')
    else
        set(gca,'ytick',[])
    end
    
    subplot(2,3,trial+3) % sorted based on repeat #8
    current_trial = valid_cells_sess1(:,:,repeat_list(trial));
    current_trial_sorted = current_trial(sorted_ind3,:);
    imagesc(current_trial_sorted,[0 30])
    if trial == 1
        ylabel('Sorted by repeat 8')
    elseif trial ==2
        xlabel('Time in movie')
    end
    if trial ~=1
        set(gca,'ytick',[])
    end
end
colormap(flipud(gray))
cb = colorbar('Ticks',[0 30],'TickLabels',[0 30]);
set(cb,'position',[0.925 0.55 0.02 0.31])
cb.Label.String = 'Activity rate (events/sec)';
cb.FontSize = 8;
suptitle(['Mouse #91 - area V1:'])


% visualize the internal structure of the same example mouse
load('Figure6A.mat','state') % load seed for tSNE visualization
rng(state) % set seed

% reshape the neuronal activity from 3D (475 cells x 90 time bins x 10 repeats) into 2D (475 cells x 900 time bins across repeats)
reshped_example_mouse = reshape(valid_cells_sess1,[size(valid_cells_sess1,1),size(valid_cells_sess1,2)*size(valid_cells_sess1,3)]);
timebin_colors = repmat([1:90],[1,10]); % define a colorscheme vector for time bins

dim_reduce = tsne(reshped_example_mouse','Algorithm','exact','Distance','cosine',...
    'NumDimensions',2,'NumPCAComponents',20,'Perplexity',100); % tSNE dim reduction on neuronal activity of example mouse

figure % visualization of internal structure
scatter(-dim_reduce(:,1),dim_reduce(:,2),15,timebin_colors,'filled')
set(gca, 'YAxisLocation','right')
colormap(new_jet_colormap)
xlabel('Component 1')
ylabel('Component 2')
text(0.65, 0.075,['Natural movie 1'],'Units','normalized','color',[0.6 0.6 0.6],'fontsize',15)
xlim([-15 15])
ylim([-15 15])
cb = colorbar('Ticks',[1 90],'TickLabels',{'Start','End'});
set(cb,'position',[0.2 0.55 0.04 0.35])
cb.Label.String = 'Time in movie';
cb.FontSize = 12;


%% Figure 6C - internal structures of a single neuropixels mouse

nat_movie = 1; % natural movie 1
num_repeats = 30; % number of repeats per block ('Functional connectivity' group)
mouse = 35; % example mouse #35

structures_area_mats = []; % define an empty variable that will store the zectorized internal structures of all visual areas (will be in the size of 435 x 360)
structures_area_labels = []; % define an empty variable that will store the label of each vectorize internal structure to its corresponding visual area
triu_ind = boolean(triu(ones(30),1)); % define a boolian matrix with true values in the upper half of the matrix
for area = 1:6 % loop over visual areas
    current_area = neuropixels_population_vectors{mouse,area,nat_movie}; % subset the neuronal activity of a single area from the exmaple mouse
    for repeat = 1:size(current_area,3) % loop over movie repeats
        current_repeat = current_area(:,:,repeat); % subset the neuronal activity of a single movie repeat
        current_structure = corr(current_repeat); % calculate the internal struture for a single movie repeat (correlation between time bins)
        structures_area_mats = [structures_area_mats,current_structure(triu_ind)]; % vectorize and store the internal structure
        structures_area_labels = [structures_area_labels,area]; % store the area label of the current internal structure
    end
end

load('Figure6C.mat','state') % load seed for tSNE visualization
rng(state) % set seed

dim_reduce = tsne(structures_area_mats','Algorithm','exact','Distance','Cosine',...
    'NumDimensions',3,'NumPCAComponents',10,'Perplexity',30); % perform tSNE dim reduction on the internal structures of a single mouse

figure('units','normalized','position',[0.4 0.55 0.225 0.35]) % visualize the relationship between internal structures of a single mouse
for area = 1:6
    scatter3(-dim_reduce(structures_area_labels==area,2),dim_reduce(structures_area_labels==area,1),dim_reduce(structures_area_labels==area,3),50,colors(area,:),'filled')
    hold on
end
grid on
ax = gca;
ax.GridAlpha = 0.1;
ylim([-25 25])
xlabel('Component 2')
ylabel('Component 1')
zlabel('Component 3')
legend(brain_areas,'Location','best')
legend('boxoff')
view([15 10 2.5])

% visualize two exmaple internal structure matrices
figure('units','normalized','position',[0.4 0.25 0.1 0.175]) 
imagesc(corr(neuropixels_population_vectors{mouse,1,nat_movie}(:,:,49))) % repeat #49 from area V1
xlabel('Time in movie (sec)')
ylabel('Time in movie (sec)')
title('V1 - repeat 49')
colormap(newmap3)

figure('units','normalized','position',[0.525 0.25 0.1 0.175])
imagesc(corr(neuropixels_population_vectors{mouse,3,nat_movie}(:,:,37))) % repeat #37 from area AL
xlabel('Time in movie (sec)')
ylabel('Time in movie (sec)')
title('AL - repeat 37')
colormap(newmap3)

%% Figure 6D - Internal structures of example neuropixels pseudo-mouse

% Creation of two example pseudo-mice using the "Functional connectivity" group.
% Since this procedure involves the random subsampling of mice and cells it requires the loading of saved seed
load('Figure6D.mat','state') % load seed
rng(state) % set seed

nat_movie = 1; % natural movie 1
subset_population_vectors = neuropixels_population_vectors(movie_repeats(:,nat_movie) == 30,1:6,nat_movie); % subset only mice from the "Functional connectivity" group

% split the dataset into two independent group of mice
pseudo_mouseA_ind = sort(randperm(size(subset_population_vectors,1),size(subset_population_vectors,1)./2)); % indices for group A (pseudo-mouse A)
pseudo_mouseB_ind = find(~ismember([1:size(subset_population_vectors,1)],pseudo_mouseA_ind)); % indices for group B (pseudo-mouse B)

% for each visual area in each pseudo mouse, pool all units across mice
% to create 12 pseudo areas (6 areas x 2 pseudo-mice)
pseudo_area_cell_num = []; % define an empty variable that will store the number of units in each pseudo-area
pseudo_mouseA = {}; % define an empty variable that will store the pooled units across mice for each pseudo-area of pseudo-mouse A
pseudo_mouseB = {}; % define an empty variable that will store the pooled units across mice for each pseudo-area of pseudo-mouse B
for area = 1:6 % loop over areas
    pseudo_mouseA{area} = cell2mat(subset_population_vectors(pseudo_mouseA_ind,area,nat_movie)); % pool units across mice for a single visual area (pseudo-mouse A)  
    pseudo_mouseB{area} = cell2mat(subset_population_vectors(pseudo_mouseB_ind,area,nat_movie)); % pool units across mice for a single visual area (pseudo-mouse B) 
    pseudo_area_cell_num(1,area) = size(pseudo_mouseA{area},1); % store number of cells in pseudo-area of pseudo-mouse A
    pseudo_area_cell_num(2,area) = size(pseudo_mouseB{area},1); % store number of cells in pseudo-area of pseudo-mouse B
end

% min_cell_num - the minimal number of cells across pseudo-areas and pseudo-mice.
% this number will be used to randomly subsample the same number of cells for all pseudo-ares 
min_cell_num = min(pseudo_area_cell_num(:)); 

% subsampling randomly the same number of units for all pseudo-ares  
pseudo_mouseA_subset = {}; % define an empty variable that will store the subsampled units for each pseudo-area of pseudo-mouse A
pseudo_mouseB_subset = {}; % define an empty variable that will store the subsampled units for each pseudo-area of pseudo-mouse B
for area = 1:6 % loop over areas
    % random sampling #min_cell_num of units from each pseudo-area of pseudo-mouse A
    subset_cell_ids_mouseA = sort(randperm(pseudo_area_cell_num(1,area),min_cell_num)); 
    pseudo_mouseA_subset{area} = pseudo_mouseA{area}(subset_cell_ids_mouseA,:,:);
    
    % random sampling #min_cell_num of units from each pseudo-area of pseudo-mouse B
    subset_cell_ids_mouseB = sort(randperm(pseudo_area_cell_num(2,area),min_cell_num));
    pseudo_mouseB_subset{area} = pseudo_mouseB{area}(subset_cell_ids_mouseB,:,:);
end

% calculating the internal structures for each movie repeat for all
% pseudo-areas of both example pseudo-mice
internal_structures_pseudoA = []; % define an empty variable the will store the internal structures of pseudo-mouse A
internal_structures_pseudoB = []; % define an empty variable the will store the internal structures of pseudo-mouse B
internal_structures_labels = []; % define an empty variable that will store the label of each internal structure

triu_ind = boolean(triu(ones(30),1)); % define a boolian matrix with true values in the upper half of the matrix
for area = 1:6 % loop over areas
    for repeat = 1:60 % loop over movie repeats
        current_structure_mouseA = corr(pseudo_mouseA_subset{area}(:,:,repeat)); % calculate the internal struture for a single movie repeat (correlation between time bins) for pseudo-mouse A
        internal_structures_pseudoA = [internal_structures_pseudoA,current_structure_mouseA(triu_ind)]; % vectorize and store the internal structure
        
        current_structure_mouseB = corr(pseudo_mouseB_subset{area}(:,:,repeat)); % calculate the internal struture for a single movie repeat (correlation between time bins) for pseudo-mouse B
        internal_structures_pseudoB = [internal_structures_pseudoB,current_structure_mouseB(triu_ind)]; % vectorize and store the internal structure
        
        internal_structures_labels = [internal_structures_labels,area]; % store the area label of the current internal structure
    end
end

all_internal_structures = [internal_structures_pseudoA,internal_structures_pseudoB]; % combine the internal structures of both pseud-mice into a single matrix (will be used for dim reduction)
all_internal_structures_labels = [internal_structures_labels,internal_structures_labels+6]; % create a vector containing the labels for all internal structures of both pseudo-mice

dim_reduce = tsne(all_internal_structures','Algorithm','exact','Distance','Cosine',...
    'NumDimensions',3,'NumPCAComponents',20,'Perplexity',30); % perform tSNE dim reduction on the internal structures of both pseudo-mice

plt = []; % define an empty variable that will store the plot information for each internal structure of each pseudo-area
figure('units','normalized','position',[0.3 0.3 0.215 0.355]) % visulaize the relationship between the internal structures of emaple pseudo-mice
for area = 1:6 % loop over areas
    current_structuresA = [all_internal_structures_labels ==area];
    current_structuresB = [all_internal_structures_labels ==area+6];
    hold on
    plt(area) = scatter3(dim_reduce(current_structuresA,3),dim_reduce(current_structuresA,2),dim_reduce(current_structuresA,1),60,colors(area,:),'filled');
    scatter3(dim_reduce(current_structuresB,3),dim_reduce(current_structuresB,2),dim_reduce(current_structuresB,1),60,colors2(area,:),'filled')
end

grid on
ax = gca;
ax.GridAlpha = 0.1;
xlabel('Component 1')
ylabel('Component 2')
zlabel('Component 3')
view([2 24.4])
legend(plt,brain_areas(1:6),'Location','best')
legend('boxoff')

%% Figure 6E - Internal structure similarity across technologies

% pool over all cells of each visual area of each dataset in order to
% create two pseudo mice (one for each of the two recordings techniques)
pseudo_mice = {}; % define an empty varialbe that will store the pooled cells of each area for each pseudo-mouse (1st row - calcium imaging, 2nd row - neuropixels)
for area = 1:6 % loop over areas
    pseudo_mice{1,area} = cell2mat(calcium_excitatory_population_vectors{area}(:,1)); % pool all cells of a single visual area from the calcium-imaging dataset
    
    valid_mice = neuropixels_cell_count(:,area,1) >0; % include only mice that were recorded from the currect visual area
    current_area = neuropixels_population_vectors(valid_mice,area); % subset all mice that passed the requirments of 'valid_mice'
    pseudo_mice{2,area} = cell2mat(cellfun(@(x) x(:,:,1:20),current_area,'UniformOutput',false)); % pool all units from a single visual area from the neuropixels dataset including neuronal activity during the first 20 movie repeats
end

calcium_internal_structure = []; % create an empty variable that will store the internal structures of the calcium imaging pseudo-mouse
calcium_area_labels = []; % create an empty variable that will store the labels for calcium imaging internal structures

neuropixels_internal_structure = []; % create an empty variable that will store the internal structures of the calcium imaging pseudo-mouse
neuropixels_area_labels = []; % create an empty variable that will store the labels for neuropixels internal structures

triu_ind = boolean(triu(ones(30),1)); % define a boolian matrix with true values in the upper half of the matrix
for area = 1:6 % loop over areas
    calcium_current_area = pseudo_mice{1,area}; % subset the neuronal activity of a single pseudo-area of the calcium imaging pseudo-mouse
    internal_structure_area = []; % define an empty variable that will store the internal structures of the current calcium imaging pseudo-area
    for repeat = 1:30 % loop over movie repeats
        current_structure = corr(calcium_current_area(:,:,repeat)); % calculate the internal structure of a single movie repeat (correlation across time bins)
        internal_structure_area(:,repeat) = current_structure(triu_ind); % vectorize and store the internal structure
    end
    
    calcium_internal_structure = [calcium_internal_structure,internal_structure_area]; % store the internal structures of the current calcium imaging pseudo-area
    calcium_area_labels = [calcium_area_labels,ones(1,30)*area]; % create and store the labels for the current calcium imaging pseud-area
    
    neuropixels_current_area = pseudo_mice{2,area};  % subset the neuronal activity of a single pseudo-area of the neuropixels pseudo-mouse
    internal_structure_area = []; % define an empty variable that will store the internal structures of the current neuropixels pseudo-area
    for repeat = 1:20 % loop over movie repeats
        current_structure = corr(neuropixels_current_area(:,:,repeat)); % calculate the internal structure of a single movie repeat (correlation across time bins)
        internal_structure_area(:,repeat) = current_structure(triu_ind); % vectorize and store the internal structure
    end
    neuropixels_internal_structure = [neuropixels_internal_structure,internal_structure_area]; % store the internal structures of the current neuropixels pseudo-area
    neuropixels_area_labels = [neuropixels_area_labels,ones(1,20)*area]; % create and store the labels for the current calcium imaging pseud-area
end

% normalization of internal structures within each dataset.
% without this normalization the tSNE dim reduction will seperate the
% internal structures into two cluster corresponding to the different
% technologies. in order to examine the similarities between technologies
% the it is required to normalize the range of values within each dataset.
neuropixels_internal_structure_zscore = zscore(neuropixels_internal_structure,[],2); % normalize (z-score) the internal structures of the neuropixels pseudo-mouse
calcium_internal_structure_zscore = zscore(calcium_internal_structure,[],2); % normalize (z-score) the internal structures of the calcium imaging pseudo-mouse

% calculate the median internal structure for each pseudo-area of each pseudo-mouse
calcium_median_internal_structure= []; % define an empty variable that will store the median internal structure of each pseudo-area of the calcium imaging pseudo-mouse
neuropixels_median_internal_structure = []; % define an empty variable that will store the median internal structure of each pseudo-area of the neuropixels pseudo-mouse
for area = 1:6
    current_area_calcium = calcium_internal_structure_zscore(:,calcium_area_labels == area); % subset the normalized internal structure of a single pseudo-area of the calcium imaging pseudo-mouse
    calcium_median_internal_structure(:,area)=  median(current_area_calcium,2,'omitnan'); % calculate the median internal structure for a single pseudo-area (calcium imaging)
    
    current_area_neuropixels = neuropixels_internal_structure_zscore(:,neuropixels_area_labels == area); % subset the normalized internal structure of a single pseudo-area of the neuropixels pseudo-mouse
    neuropixels_median_internal_structure(:,area) = median(current_area_neuropixels,2,'omitnan'); % calculate the median internal structure for a single pseudo-area (neruopixels)
end


normalized_internal_structures_both_tech = [neuropixels_internal_structure_zscore,calcium_internal_structure_zscore]; % combine the normalized median internal structures of both neuropixels and calcium imaging pseudo-mice into a single matrix
area_labels_both_tech =[neuropixels_area_labels,calcium_area_labels+6]; % create a vector with labels for the internal structures of both pseudo-mice


load('figure6E.mat','state') % load seed for tSNE visualization
rng(state) % set seed
dim_reduce = tsne(normalized_internal_structures_both_tech','Algorithm','exact','Distance','correlation',...
    'NumDimensions',3,'NumPCAComponents',20,'Perplexity',40); % perform tSNE dim reduction on the internal structures of both pseudo-mice
plt = []; % define empty variable to store plot information for each visual area of each pseudo-mouse
figure('units','normalized','position',[0.3 0.3 0.25 0.4]) % visualize the relationship between internal structures of both pseudo-mice
for area =  1:6 % loop over visual areas
    hold on
    current_dimA = dim_reduce(area_labels_both_tech==area,1:3);
    current_dimB = dim_reduce(area_labels_both_tech==area+6,1:3);
    mean_metaA = nanmedian(current_dimA);
    mean_metaB = nanmedian(current_dimB);
    
    
    scatter3(current_dimA(:,2),current_dimA(:,1),current_dimA(:,3),20,colors(area,:),'filled','MarkerfaceAlpha',0.3)
    scatter3(current_dimB(:,2),current_dimB(:,1),current_dimB(:,3),20,colors2(area,:),'filled','MarkerfaceAlpha',0.3)
    
    plt(area) =  scatter3(mean_metaA(:,2),mean_metaA(:,1),mean_metaA(:,3),200,colors(area,:),'filled');
    scatter3(mean_metaB(:,2),mean_metaB(:,1),mean_metaB(:,3),200,colors2(area,:),'filled');
    plot3([mean_metaA(:,2),mean_metaB(:,2)],[mean_metaA(:,1),mean_metaB(:,1)],[mean_metaA(:,3),mean_metaB(:,3)],'color',colors(area,:),'linewidth',3)
    
end
legend(plt,brain_areas(1:6))
legend('boxoff')
grid on
ax = gca;
ax.GridAlpha = 0.1;
xlabel('Component 2')
ylabel('Component 1')
zlabel('Component 3')
set(gca, 'Zdir', 'reverse')
view([15 10 3])

% calculate the correlation distance between the median internal structures
% of the two pseudo-mice
similarity_across_technologies = pdist2(calcium_median_internal_structure',neuropixels_median_internal_structure','correlation');
figure('units','normalized','position',[0.55 0.4 0.2 0.275])
imagesc(similarity_across_technologies) % visualize the correlation distance matrix between internal structures of the two pseudo-mice
colormap(flipud(newmap3))
set(gca,'xtick',1:6,'xticklabel',brain_areas(1:6),'ytick',1:6,'yticklabel',brain_areas(1:6))
xlabel('Calcium imaging')
ylabel('Neuropixels')
cb = colorbar;
cb.Label.String = 'Correlation distance';
cb.FontSize = 10;


%% Figure 6F - Internal structures of two pseudo-mice (reduced and matrices)

nat_movie = 1; % natural movie 1
subset_population_vectors = {}; % define an empty variable that will variable that will store the neuronal activity of mice from the 'Function connectivity' group
for area = 1:6 % loop over areas
    for mouse = 1:size(neuropixels_population_vectors_tsne,1) % loop over mice
        if ~isempty(neuropixels_population_vectors_tsne{mouse,area,nat_movie}) && movie_repeats(mouse,nat_movie) == 30 % check if mouse was recorded from the currvent area and from the "Functional connectivity" group
            subset_population_vectors{mouse,area} = neuropixels_population_vectors_tsne{mouse,area,nat_movie}(:,:,31:60); % subset the neuronal activity of a single mouse during block B
        end
        
    end
end
% split the dataset into two independent group of mice
load('neuropixels_structure_ind.mat','pseudo_mouseA_ind','pseudo_mouseB_ind') % load pseudo-mice indices

% for each visual area in each pseudo mouse, pool all units across mice
% to create 12 pseudo areas (6 areas x 2 pseudo-mice)
pseudo_mouseA = {}; % define an empty variable that will store the pooled units across mice for each pseudo-area of pseudo-mouse A
pseudo_mouseB = {}; % define an empty variable that will store the pooled units across mice for each pseudo-area of pseudo-mouse B
for area = 1:6 % loop over areas
    pseudo_mouseA{area} = cell2mat(subset_population_vectors(pseudo_mouseA_ind,area)); % pool units across mice for a single visual area (pseudo-mouse A)  
    pseudo_mouseB{area} = cell2mat(subset_population_vectors(pseudo_mouseB_ind,area)); % pool units across mice for a single visual area (pseudo-mouse B) 
end

% perform tsne dim reduction for each pseudo-area of each pseudo-mouse
dim_reduce_areas = {}; % define an empty variable that will store the dim reduction components of each pseudo-area of each pseud-mouse
for area = 1:6 % loop over visual areas
    load(['neuropixels_structure_STATE',num2str(area),'.mat'],'state') % load seed for tSNE analysis of a single area
    rng(state) % set seed
    
    % reshape the neuronal activity of pseudo-mouse A from 3D (#cells x 90 time bins x 30 repeats) into 2D (#cells x 27000 time bins across repeats)
    pseudo_mouseA_strcut = reshape(pseudo_mouseA{area},[size(pseudo_mouseA{area},1),90*30]); 
    
    dim_reduceA = tsne(pseudo_mouseA_strcut','Algorithm','exact','Distance','cosine',...
        'NumDimensions',3,'NumPCAComponents',20,'Perplexity',200); % tSNE dim reduction on neuronal activity of a single visual area of pseudo-mouse A
    dim_reduce_areas{1,area} = dim_reduceA; % store the tSNE components for the current pseudo-area of pseudo-mouse A
    
    % reshape the neuronal activity of pseudo-mouse B from 3D (#cells x 90 time bins x 30 repeats) into 2D (#cells x 27000 time bins across repeats)
    pseudo_mouseB_strcut = reshape(pseudo_mouseB{area},[size(pseudo_mouseB{area},1),90*30]);
    
    dim_reduceB = tsne(pseudo_mouseB_strcut','Algorithm','exact','Distance','cosine',...
        'NumDimensions',3,'NumPCAComponents',20,'Perplexity',200); % tSNE dim reduction on neuronal activity of a single visual area of pseudo-mouse B
    dim_reduce_areas{2,area} = dim_reduceB; % store the tSNE components for the current pseudo-area of pseudo-mouse B
    
end

% visualization of reduced space internal structures of all six visual
% areas of both neuropixels pseudo-mice
view_list_pseudoA = [-127.1 24.4;-174.7 63.6;4.9 19.6;-96.3 57.2;-22.7 6;34.5 52.4]; % defining the 3D rotation of each plot of pseudo-mouse A
view_list_pseudoB = [55.3 -32.4;79.7 70.8;-3.5 -25.2;133.2 26;31.6 11.2;61.2 13.2]; % defining the 3D rotation of each plot of pseudo-mouse B
figure('units','normalized','position',[0.1 0.5 0.8 0.4]) % visualize the internal structures
for area = 1:6 % loop over areas
    subplot(2,6,area)
    scatter3(dim_reduce_areas{1,area}(:,1),dim_reduce_areas{1,area}(:,2),dim_reduce_areas{1,area}(:,3),5,repmat([1:90],[1 30]),'filled')
    view(view_list_pseudoA(area,:))
    set(gca,'xtick',[],'ytick',[],'ztick',[])
    grid off
    subplot(2,6,area+6)
    scatter3(dim_reduce_areas{2,area}(:,1),dim_reduce_areas{2,area}(:,2),dim_reduce_areas{2,area}(:,3),5,repmat([1:90],[1 30]),'filled')
    view(view_list_pseudoB(area,:))
    set(gca,'xtick',[],'ytick',[],'ztick',[])
    grid off
end
colormap(new_jet_colormap)

figure('units','normalized','position',[0.1 0 0.8 0.4]) % visualize the internal structures in their correlation matrix form
for area = 1:6 % loop over visual areas
    current_area_mouseA = mean(cell2mat(pseudo_mouseA(area)),3,'omitnan'); % calculate the average neuronal activity across movie repeats for a single visual area of pseudo-mouse A
    pseudo_mouseA_strcut = corr(current_area_mouseA); % calculate the internal struture for current pseudo-area (correlation between time bins)
    
    current_area_mouseB = mean(cell2mat(pseudo_mouseB(area)),3,'omitnan'); % calculate the average neuronal activity across movie repeats for a single visual area of pseudo-mouse A
    pseudo_mouseB_strcut = corr(current_area_mouseB); % calculate the internal struture for current pseudo-area (correlation between time bins)
    
    subplot(2,6,area)
    imagesc(pseudo_mouseA_strcut)
    subplot(2,6,area+6)
    imagesc(pseudo_mouseB_strcut)
end
colormap(newmap3)

%% Figure 6G - between pseudo-mice permutation decoder

nat_movie = 1; % natural movie 1

% subsample for each mouse in the neuropixels dataset the neuronal activity during the first 20 movie repeats
subset_population_vectors = {}; % define an empty variable that will store the subsmaple neuronal activity of indevidual mice
for area = 1:6 % loop over areas
    for mouse = 1:size(neuropixels_population_vectors,1) % loop over mice
        if ~isempty(neuropixels_population_vectors{mouse,area,nat_movie}) % test if current mouse was recorded from current area
            subset_population_vectors{mouse,area} = neuropixels_population_vectors{mouse,area,nat_movie}(:,:,1:20); % subset the neuronal activity during the first 20 movie repeats
        end
    end
end

% decode the identity of each visual area based on the
% similarity between the internal structures of the two psedo-mice
similarity_between_pseudomice = []; % define an empty variable that will store the decoder results for the non-shuffled pseudo-mice
similarity_between_shuffled_pseudomice = []; % define an empty variable that will store the decoder results for the huffled pseudo-mice

num_shuffles = 1000; % number of pseudo-mice realizations
for shuffle = 1:num_shuffles % loop over realizations
    clc;
    disp(['Performing between pseudo-mice decoding. Realization: ',num2str(shuffle),'/',num2str(num_shuffles)])
    
    % split the dataset into two independent group of mice
    pseudo_mouseA_ind = sort(randperm(size(subset_population_vectors,1),size(subset_population_vectors,1)./2)); % indices for group A (pseudo-mouse A)
    pseudo_mouseB_ind = find(~ismember([1:size(subset_population_vectors,1)],pseudo_mouseA_ind)); % indices for group B (pseudo-mouse B)
    
    % for each visual area in each pseudo mouse, pool all units across mice
    % to create 12 pseudo areas (6 areas x 2 pseudo-mice)
    pseudo_area_cell_num = []; % define an empty variable that will store the number of units in each pseudo-area
    pseudo_mouseA = {}; % define an empty variable that will store the pooled units across mice for each pseudo-area of pseudo-mouse A
    pseudo_mouseB = {}; % define an empty variable that will store the pooled units across mice for each pseudo-area of pseudo-mouse B
    for area = 1:6 % loop over areas
        pseudo_mouseA{area} = cell2mat(subset_population_vectors(pseudo_mouseA_ind,area,nat_movie)); % pool units across mice for a single visual area (pseudo-mouse A)  
        pseudo_mouseB{area} = cell2mat(subset_population_vectors(pseudo_mouseB_ind,area,nat_movie)); % pool units across mice for a single visual area (pseudo-mouse B) 
        pseudo_area_cell_num(1,area) = size(pseudo_mouseA{area},1); % store number of cells in pseudo-area of pseudo-mouse A
        pseudo_area_cell_num(2,area) = size(pseudo_mouseB{area},1); % store number of cells in pseudo-area of pseudo-mouse B
    end
    
    % min_cell_num - the minimal number of cells across pseudo-areas and pseudo-mice.
    % this number will be used to randomly subsample the same number of cells for all pseudo-ares 
    min_cell_num = min(pseudo_area_cell_num(:));
    
    % subsampling randomly the same number of units for all pseudo-areas
    pseudo_mouseA_subset = {}; % define an empty variable that will store the subsampled units for each pseudo-area of pseudo-mouse A
    pseudo_mouseB_subset = {}; % define an empty variable that will store the subsampled units for each pseudo-area of pseudo-mouse B
    for area = 1:6 % loop over areas
        % random sampling #min_cell_num of units from each pseudo-area of pseudo-mouse A
        subset_cell_ids_mouseA = sort(randperm(pseudo_area_cell_num(1,area),min_cell_num));
        pseudo_mouseA_subset{area} = pseudo_mouseA{area}(subset_cell_ids_mouseA,:,:);
        
        % random sampling #min_cell_num of units from each pseudo-area of pseudo-mouse B
        subset_cell_ids_mouseB = sort(randperm(pseudo_area_cell_num(2,area),min_cell_num));
        pseudo_mouseB_subset{area} = pseudo_mouseB{area}(subset_cell_ids_mouseB,:,:);
    end
    
    % creating shuffled pseudo mice
    all_cells_pseudo_mouseA = cell2mat(pseudo_mouseA_subset'); % pooling all the cells across all pseudo-areas of pseudo-mouse A
    all_cells_pseudo_mouseB = cell2mat(pseudo_mouseB_subset'); % pooling all the cells across all pseudo-areas of pseudo-mouse B

     % random redistribution of cells across pseudo-areas within a given pseudo-mouse
    rand_cells_id_mouseA = randperm(size(all_cells_pseudo_mouseA,1)); % random permutation of cells indices for pseudo-mouse A
    rand_cells_id_mouseB = randperm(size(all_cells_pseudo_mouseB,1)); % random permutation of cells indices for pseudo-mouse B

    shuffled_pseudo_mouseA_subset = {}; % define an empty variable that will store the redistributed cells across areas for pseudo-mouse A
    shuffled_pseudo_mouseB_subset = {}; % define an empty variable that will store the redistributed cells across areas for pseudo-mouse B
    for area = 1:6 % loop over areas
        current_pseudo_area = [1:min_cell_num] + min_cell_num*(area-1); % define range of cell indices for the current area
        
        shuffled_pseudo_mouseA_subset{area} = all_cells_pseudo_mouseA(rand_cells_id_mouseA(current_pseudo_area),:,:); % randomly subsample #min_cell_num cells to current pseud-area of pseudo-mouse A
        shuffled_pseudo_mouseB_subset{area} = all_cells_pseudo_mouseB(rand_cells_id_mouseB(current_pseudo_area),:,:); % randomly subsample #min_cell_num cells to current pseud-area of pseudo-mouse B
    end
    
    % calculating the internal structures for each movie repeat for all
    % pseudo-areas of both example pseudo-mice
    internal_structures_pseudoA = []; % define an empty variable the will store the internal structures of pseudo-mouse A
    internal_structures_pseudoB = []; % define an empty variable the will store the internal structures of pseudo-mouse B
    internal_structures_pseudoA_shuffle = []; % define an empty variable the will store the internal structures of shuffled pseudo-mouse A
    internal_structures_pseudoB_shuffle = []; % define an empty variable the will store the internal structures of shuffled pseudo-mouse B
    internal_structures_labels = []; % define an empty variable that will store the label of each internal structure
    triu_ind = boolean(triu(ones(30),1)); % define a boolian matrix with true values in the upper half of the matrix
    for area = 1:6 % loop over areas
        current_area_mouseA = mean(pseudo_mouseA_subset{area},3,'omitnan'); % calculating the average neuronal activity across movie repeats for a single area in pseudo-mouse A
        current_structure_mouseA = corr(current_area_mouseA); % calculate the internal struture for current area of pseudo-mouse A
        internal_structures_pseudoA(:,area) = current_structure_mouseA(triu_ind); % vectorize and store the internal structure
        
        current_area_mouseB = mean(pseudo_mouseB_subset{area},3,'omitnan'); % calculating the average neuronal activity across movie repeats for a single area in pseudo-mouse B
        current_structure_mouseB = corr(current_area_mouseB); % calculate the internal struture for current area of pseudo-mouse B
        internal_structures_pseudoB(:,area) = current_structure_mouseB(triu_ind); % vectorize and store the internal structure
        
        current_area_mouseA_shuffle = mean(shuffled_pseudo_mouseA_subset{area},3,'omitnan'); % calculating the average neuronal activity across movie repeats for a single area in shuffle pseudo-mouse A
        current_structure_mouseA_shuffle = corr(current_area_mouseA_shuffle); % calculate the internal struture for current area of shuffle pseudo-mouse A
        internal_structures_pseudoA_shuffle(:,area) = current_structure_mouseA_shuffle(triu_ind); % vectorize and store the internal structure
        
        current_area_mouseB_shuffle = mean(shuffled_pseudo_mouseB_subset{area},3,'omitnan'); % calculating the average neuronal activity across movie repeats for a single area in shuffle pseudo-mouse B
        current_structure_mouseB_shuffle = corr(current_area_mouseB_shuffle); % calculate the internal struture for current area of shuffle pseudo-mouse B
        internal_structures_pseudoB_shuffle(:,area) = current_structure_mouseB_shuffle(triu_ind); % vectorize and store the internal structure
    end
    
   
    % calculating the similarity between the internal structures of a
    % reference mouse (pseudo-mouse A) and all 720 permutations of the test
    % mouse (pseudo-mouse B)
    similarity_sum = []; % define an empty variable that will store the total similarity (correlation sum) between non-shuffled pseudo-mice
    similarity_sum_shuffled = []; % define an empty variable that will store the total similarity (correlation sum) between shuffled pseudo-mice
    permutations = flipud(perms([1:6])); % define a matrix with all possible permutations
    for perm = 1:size(permutations,1) % loop over permutations
        current_perm = permutations(perm,:); % set the current permutations
        similarity_between_pseudomice = corr(internal_structures_pseudoA,internal_structures_pseudoB(:,current_perm)); % calculate  the correlation between the internal structures of reference and permutated non-shuffled pseudo-mice
        similarity_between_pseudomice_shuffled = corr(internal_structures_pseudoA_shuffle,internal_structures_pseudoB_shuffle(:,current_perm)); % calculate  the correlation between the internal structures of reference and permutated shuffled pseudo-mice
        
        similarity_sum(perm) = sum(diag(similarity_between_pseudomice)); % calculate the sum of correlations between corresponding pseudo-areas for non-shuffled pseudo-mice
        similarity_sum_shuffled(perm) = sum(diag(similarity_between_pseudomice_shuffled)); % calculate the sum of correlations between corresponding pseudo-areas for shuffled pseudo-mice
    end
    
    [B,I] = max(similarity_sum); % find the permutation with the highest similarity between non-shuffled pseudo-mice
    between_similarity_acc(shuffle,:) = permutations(I,:) == [1:6]; % assess decoder prediction for non-shuffled pseudo-mice
    
    [B,I] = max(similarity_sum_shuffled); % find the permutation with the highest similarity between shuffled pseudo-mice
    between_similarity_acc_shuffled(shuffle,:) = permutations(I,:) == [1:6]; % assess decoder prediction for shuffled pseudo-mice
    
end

% calculate decoder overall performance across all realizations
between_pseudo_mice_decoder_acc = (sum(between_similarity_acc)./size(between_similarity_acc,1))*100; % non-shuffled pseudo-mice
between_shuffled_pseudo_mice_decoder_acc = (sum(between_similarity_acc_shuffled)./size(between_similarity_acc_shuffled,1))*100; % shuffled pseudo-mice

plt = []; % define empty variable to store plot information for decoder performance
figure('units','normalized','position',[0.3 0.3 0.2 0.275]) % visualize decoder performance
hold on
plt(1) = bar(between_pseudo_mice_decoder_acc,'facecolor',[0.8 0.8 0.8],'edgecolor','none');
plt(2) = bar(between_shuffled_pseudo_mice_decoder_acc,'facecolor',[0.5 0.5 0.5],'edgecolor','none');
hold on
plot([0 6.75],[100 100]./6,'--','color',[0.2 0.2 0.2],'linewidth',2)
text(0.85,0.175,'Chance','Units','normalized','FontSize',11)
xlim([0 8])
set(gca,'xtick',1:6,'xticklabel',brain_areas(1:6))
ylabel({'Successful classifications ';'between pseudo-mice (%)'})
lgd = legend(plt,{'Real','Shuffle'});
legend('boxoff')
lgd.Position = [0.85 0.85 0.05 0.05];

%% Figure 6H - Between pseudo-mice permutation decoder - cell count - takes approximatly 2 hours to run

nat_movie = 1; % natural movie 1

% subsample for each mouse in the neuropixels dataset the neuronal activity during the first 20 movie repeats
subset_population_vectors = {}; % define an empty variable that will store the subsmaple neuronal activity of indevidual mice
for area = 1:6 % loop over areas
    for mouse = 1:size(neuropixels_population_vectors,1) % loop over mice
        if ~isempty(neuropixels_population_vectors{mouse,area,nat_movie}) % test if current mouse was recorded from current area
            subset_population_vectors{mouse,area} = neuropixels_population_vectors{mouse,area,nat_movie}(:,:,1:20); % subset the neuronal activity during the first 20 movie repeats
        end
    end
end

% decode the identity of each visual area based on the
% similarity between the internal structures of the two psedo-mice
cell_cutoff_list = [5,10,25,50,75,100,150,200,250,300,350,400,500,600,700]; % define list of how many cells are included in the analysis
between_similarity_acc_all_cutoffs = []; % define an empty variable that will store the decoder results for the non-shuffled pseudo-mice
between_similarity_acc_shuffled_all_cutoffs = []; % define an empty variable that will store the decoder results for the suffled pseudo-mice

num_shuffles = 3000; % number of pseudo-mice realizations for each cell count threshold
min_cell_num_all_shuffles = []; % define empty variable to store number of cells sampled in leach realization of pseudo-mice
for cells_included = 1:length(cell_cutoff_list)+1 % loop over cell count thresholds
    between_similarity_acc = nan(num_shuffles,6); % define empty variable that store decoders predictions for non-shuffled pseudo-mice
    between_similarity_acc_shuffled = nan(num_shuffles,6); % define empty variable that store decoders predictions for shuffled pseudo-mice
    
    for shuffle = 1:num_shuffles % loop over realizations
        
        % determine the number of cells to be incuded in the analysis
        if cells_included <= length(cell_cutoff_list) % if below 700 cells than use 'cell_cutoff_list'
            breaker = cell_cutoff_list(cells_included);
        else %used the entire dataset
            breaker = 1;
        end
        
        % create realizations of pseudo-mice until all pseudo-areas contain 
        % more cells relative to 'cell_cutoff_list'
        pseudo_area_cell_num = zeros(2,6); % define a variable for of zeros that will store the number of units in each pseudo-area
        while   breaker > min(pseudo_area_cell_num(:)) % check if all pseudo-areas contain more cells relative to 'cell_cutoff_list'
            
            % split the dataset into two independent group of mice
            pseudo_mouseA_ind = sort(randperm(size(subset_population_vectors,1),size(subset_population_vectors,1)./2)); % indices for group A (pseudo-mouse A)
            pseudo_mouseB_ind = find(~ismember([1:size(subset_population_vectors,1)],pseudo_mouseA_ind)); % indices for group B (pseudo-mouse B)
            
            % for each visual area in each pseudo mouse, pool all units across mice
            % to create 12 pseudo areas (6 areas x 2 pseudo-mice)
            pseudo_mouseA = {}; % define an empty variable that will store the pooled units across mice for each pseudo-area of pseudo-mouse A
            pseudo_mouseB = {}; % define an empty variable that will store the pooled units across mice for each pseudo-area of pseudo-mouse B
            for area = 1:6 % loop over areas
                pseudo_mouseA{area} = cell2mat(subset_population_vectors(pseudo_mouseA_ind,area,nat_movie)); % pool units across mice for a single visual area (pseudo-mouse A)  
                pseudo_mouseB{area} = cell2mat(subset_population_vectors(pseudo_mouseB_ind,area,nat_movie)); % pool units across mice for a single visual area (pseudo-mouse B) 
                pseudo_area_cell_num(1,area) = size(pseudo_mouseA{area},1); % store number of cells in pseudo-area of pseudo-mouse A
                pseudo_area_cell_num(2,area) = size(pseudo_mouseB{area},1); % store number of cells in pseudo-area of pseudo-mouse B
            end
            % if the minimum value in "pseudo_area_cell_num" is higher than
            % the value in "breaker" than the loop breaks and the algorithm continues
        end
        
        % min_cell_num - the minimal number of cells across pseudo-areas and pseudo-mice.
        % this number will be used to randomly subsample the same number of cells for all pseudo-ares 
        if cells_included <= length(cell_cutoff_list) % when below 700
            min_cell_num = cell_cutoff_list(cells_included);
        else % when using the entire dataset
            min_cell_num = min(pseudo_area_cell_num(:));
            min_cell_num_all_shuffles(shuffle) = min_cell_num;
        end
        
         clc;
        disp(['Performing between Neuropixels pseudo-mice decoding. Cells included: ',num2str(min_cell_num),' | Realization: ',num2str(shuffle),'/',num2str(num_shuffles)])
       
        % subsampling randomly the same number of units for all pseudo-areas
        pseudo_mouseA_subset = {}; % define an empty variable that will store the subsampled units for each pseudo-area of pseudo-mouse A
        pseudo_mouseB_subset = {}; % define an empty variable that will store the subsampled units for each pseudo-area of pseudo-mouse B
        for area = 1:6 % loop over areas
            % random sampling #min_cell_num of units from each pseudo-area of pseudo-mouse A
            subset_cell_ids_mouseA = sort(randperm(pseudo_area_cell_num(1,area),min_cell_num));
            pseudo_mouseA_subset{area} = pseudo_mouseA{area}(subset_cell_ids_mouseA,:,:);
            
            % random sampling #min_cell_num of units from each pseudo-area of pseudo-mouse B
            subset_cell_ids_mouseB = sort(randperm(pseudo_area_cell_num(2,area),min_cell_num));
            pseudo_mouseB_subset{area} = pseudo_mouseB{area}(subset_cell_ids_mouseB,:,:);
        end
        
        % creating shuffled pseudo mice
        all_cells_pseudo_mouseA = cell2mat(pseudo_mouseA_subset'); % pooling all the cells across all pseudo-areas of pseudo-mouse A
        all_cells_pseudo_mouseB = cell2mat(pseudo_mouseB_subset'); % pooling all the cells across all pseudo-areas of pseudo-mouse B
        
        % random redistribution of cells across pseudo-areas within a given pseudo-mouse
        rand_cells_id_mouseA = randperm(size(all_cells_pseudo_mouseA,1)); % random permutation of cells indices for pseudo-mouse A
        rand_cells_id_mouseB = randperm(size(all_cells_pseudo_mouseB,1)); % random permutation of cells indices for pseudo-mouse B
        
        shuffled_pseudo_mouseA_subset = {}; % define an empty variable that will store the redistributed cells across areas for pseudo-mouse A
        shuffled_pseudo_mouseB_subset = {}; % define an empty variable that will store the redistributed cells across areas for pseudo-mouse B
        for area = 1:6 % loop over areas
            current_pseudo_area = [1:min_cell_num] + min_cell_num*(area-1); % define range of cell indices for the current area
            
            shuffled_pseudo_mouseA_subset{area} = all_cells_pseudo_mouseA(rand_cells_id_mouseA(current_pseudo_area),:,:); % randomly subsample #min_cell_num cells to current pseud-area of pseudo-mouse A
            shuffled_pseudo_mouseB_subset{area} = all_cells_pseudo_mouseB(rand_cells_id_mouseB(current_pseudo_area),:,:); % randomly subsample #min_cell_num cells to current pseud-area of pseudo-mouse B
        end
        
        % calculating the internal structures for each movie repeat for all
        % pseudo-areas of both example pseudo-mice
        internal_structures_pseudoA = []; % define an empty variable the will store the internal structures of pseudo-mouse A
        internal_structures_pseudoB = []; % define an empty variable the will store the internal structures of pseudo-mouse B
        internal_structures_pseudoA_shuffle = []; % define an empty variable the will store the internal structures of shuffled pseudo-mouse A
        internal_structures_pseudoB_shuffle = []; % define an empty variable the will store the internal structures of shuffled pseudo-mouse B
        internal_structures_labels = []; % define an empty variable that will store the label of each internal structure
        triu_ind = boolean(triu(ones(30),1)); % define a boolian matrix with true values in the upper half of the matrix
        for area = 1:6 % loop over areas
            current_area_mouseA = mean(pseudo_mouseA_subset{area},3,'omitnan'); % calculating the average neuronal activity across movie repeats for a single area in pseudo-mouse A
            current_structure_mouseA = corr(current_area_mouseA); % calculate the internal struture for current area of pseudo-mouse A
            internal_structures_pseudoA(:,area) = current_structure_mouseA(triu_ind); % vectorize and store the internal structure
            
            current_area_mouseB = mean(pseudo_mouseB_subset{area},3,'omitnan'); % calculating the average neuronal activity across movie repeats for a single area in pseudo-mouse B
            current_structure_mouseB = corr(current_area_mouseB); % calculate the internal struture for current area of pseudo-mouse B
            internal_structures_pseudoB(:,area) = current_structure_mouseB(triu_ind); % vectorize and store the internal structure
            
            current_area_mouseA_shuffle = mean(shuffled_pseudo_mouseA_subset{area},3,'omitnan'); % calculating the average neuronal activity across movie repeats for a single area in shuffle pseudo-mouse A
            current_structure_mouseA_shuffle = corr(current_area_mouseA_shuffle); % calculate the internal struture for current area of shuffle pseudo-mouse A
            internal_structures_pseudoA_shuffle(:,area) = current_structure_mouseA_shuffle(triu_ind); % vectorize and store the internal structure
            
            current_area_mouseB_shuffle = mean(shuffled_pseudo_mouseB_subset{area},3,'omitnan'); % calculating the average neuronal activity across movie repeats for a single area in shuffle pseudo-mouse B
            current_structure_mouseB_shuffle = corr(current_area_mouseB_shuffle); % calculate the internal struture for current area of shuffle pseudo-mouse B
            internal_structures_pseudoB_shuffle(:,area) = current_structure_mouseB_shuffle(triu_ind); % vectorize and store the internal structure
        end
        
        % calculating the similarity between the internal structures of a
    % reference mouse (pseudo-mouse A) and all 720 permutations of the test
    % mouse (pseudo-mouse B)
    similarity_sum = []; % define an empty variable that will store the total similarity (correlation sum) between non-shuffled pseudo-mice
    similarity_sum_shuffled = []; % define an empty variable that will store the total similarity (correlation sum) between shuffled pseudo-mice
    permutations = flipud(perms([1:6])); % define a matrix with all possible permutations
    for perm = 1:size(permutations,1) % loop over permutations
        current_perm = permutations(perm,:); % set the current permutations
        similarity_between_pseudomice = corr(internal_structures_pseudoA,internal_structures_pseudoB(:,current_perm)); % calculate  the correlation between the internal structures of reference and permutated non-shuffled pseudo-mice
        similarity_between_pseudomice_shuffled = corr(internal_structures_pseudoA_shuffle,internal_structures_pseudoB_shuffle(:,current_perm)); % calculate  the correlation between the internal structures of reference and permutated shuffled pseudo-mice
        
        similarity_sum(perm) = sum(diag(similarity_between_pseudomice)); % calculate the sum of correlations between corresponding pseudo-areas for non-shuffled pseudo-mice
        similarity_sum_shuffled(perm) = sum(diag(similarity_between_pseudomice_shuffled)); % calculate the sum of correlations between corresponding pseudo-areas for shuffled pseudo-mice
    end
    
    
        [B,I] = max(similarity_sum); % find the permutation with the highest similarity between non-shuffled pseudo-mice
    between_similarity_acc(shuffle,:) = permutations(I,:) == [1:6]; % assess decoder prediction for non-shuffled pseudo-mice
    
    [B,I] = max(similarity_sum_shuffled); % find the permutation with the highest similarity between shuffled pseudo-mice
    between_similarity_acc_shuffled(shuffle,:) = permutations(I,:) == [1:6]; % assess decoder prediction for shuffled pseudo-mice
    
    end
    between_similarity_acc_all_cutoffs(cells_included,:) = sum(between_similarity_acc)./num_shuffles; % store calculate decoder overall performance across all realizations of non-shuflled pseudo-mice
    between_similarity_acc_shuffled_all_cutoffs(cells_included,:) = sum(between_similarity_acc_shuffled)./num_shuffles; % store calculate decoder overall performance across all realizations of shuffled pseudo-mice
    
end

plt = []; % define empty variable to store plot information for decoder performance
xvalues = [cell_cutoff_list,round(mean(min_cell_num_all_shuffles,'omitnan'))]; % define x axis values
figure('units','normalized','position',[0.35 0.4 0.25 0.375]) % visualize decoder performance as function of the number of cells included
for area = 1:6 % loop over areas
    hold on
    plt(area) = plot(xvalues,between_similarity_acc_all_cutoffs(:,area)*100,'color',colors(area,:),'linewidth',3);
    plot(xvalues,between_similarity_acc_shuffled_all_cutoffs(:,area)*100,'color',[0.2 0.2 0.2]+0.1*(area-1),'linewidth',3)
end
text(0.55, 0.225,['Shuffled pseudo-mice'],'Units','normalized','color',[0.4 0.4 0.4],'fontsize',12)
set(gca,'xtick',[0:100:800],'xticklabels',[0:100:800])
ylim([0 100])
xlim([-10 820])
ylabel('Successful classifications (%)')
xlabel('Number of cells included')
legend(plt,brain_areas(1:6),'Location','best')
legend('boxoff')

%% Figure 7A - Internal structure in reduced space - single animal v1 example across sessions

load('Figure7A.mat','state') % load seed for tsne visualization
rng(state) % set seed

area = 1; % area V1
mouse = 89; % exmaple mouse #89
current_area = calcium_excitatory_population_vectors{area}{mouse,3}; % subset the neuronal activity of example mouse

% reshape the neuronal activity from 3D (475 cells x 90 time bins x 10 repeats) into 2D (475 cells x 900 time bins across repeats)
reshaped_current_area = reshape(current_area,[size(current_area,1),size(current_area,2)*size(current_area,3)]);
timebin_colors = repmat([1:90],[1,30]); % define a colorscheme vector for time bins

dim_reduce_sess1 = tsne(reshaped_current_area(:,1:900)','Algorithm','exact','Distance','Cosine',...
    'NumDimensions',2,'NumPCAComponents',10,'Perplexity',200); % perform tsne dim reduction on neuronal activity of session 1

dim_reduce_sess2 = tsne(reshaped_current_area(:,901:1800)','Algorithm','exact','Distance','Cosine',...
    'NumDimensions',2,'NumPCAComponents',10,'Perplexity',200); % perform tsne dim reduction on neuronal activity of session 2

dim_reduce_sess3 = tsne(reshaped_current_area(:,1801:2700)','Algorithm','exact','Distance','cosine',...
    'NumDimensions',2,'NumPCAComponents',10,'Perplexity',200); % perform tsne dim reduction on neuronal activity of session 3

figure('units','normalized','position',[0.25 0.4 0.45 0.225]) % visualize internal structures across sessions
subplot(1,3,1) % visualize the internal structure of session 1
scatter(dim_reduce_sess1(:,2),dim_reduce_sess1(:,1),5,timebin_colors(:,1:900),'filled')
set(gca, 'xdir', 'reverse')
xlim([-10 10])
ylim([-10 10])
text(0.1, 0.1,['Natural movie 1'],'Units','normalized','color',[0.6 0.6 0.6],'fontsize',12)
xlabel('Component 1')
ylabel('Component 2')
title({'Session 1 - Day 93'})

subplot(1,3,2) % visualize the internal structure of session 2
scatter(dim_reduce_sess2 (:,1),dim_reduce_sess2(:,2),5,timebin_colors(:,1:900),'filled')
xlim([-10 10])
ylim([-10 10])
xlabel('Component 1')
ylabel('Component 2')
title({'Session 2 - Day 94'})

subplot(1,3,3) % visualize the internal structure of session 3
scatter(dim_reduce_sess3(:,2),dim_reduce_sess3(:,1),5,timebin_colors(:,1:900),'filled')
xlim([-10 10])
ylim([-8 8])
xlabel('Component 1')
ylabel('Component 2')
title({'Session 3 - Day 98'})
set(gca, 'xdir', 'reverse', 'ydir', 'reverse')
colormap(new_jet_colormap)

%% Figure 7B - Internal structure versus PV for pseudo-AL

area = 3; % area AL
current_area = cell2mat(calcium_excitatory_population_vectors{area}(:,1)); % subset all mice recorded from area AL

active_sess1 = mean(mean(current_area(:,:,1:10),3,'omitnan'),2,'omitnan')>0; % active cells in session 1 (repeats 1-10)
active_sess2 = mean(mean(current_area(:,:,11:20),3,'omitnan'),2,'omitnan')>0; % active cells in session 1 (repeats 11-20)
active_sess3 = mean(mean(current_area(:,:,21:30),3,'omitnan'),2,'omitnan')>0; % active cells in session 1 (repeats 21-30)
valid_cells = find([active_sess1 & active_sess2 & active_sess3]); % active in all three sessions

valid_current_area = current_area(valid_cells,:,:); % subset only cells tat were active in all imaging sessions

rep_list = [1:10;21:30]; % list of sessions to include in visualization (session 1: repeats 1-10; session 3: repeats 21-30)
figure('units','normalized','position',[0.4 0.4 0.225 0.4]) % visualize internal structure vs single cells' responses across sessions
for sess = 1:2
    current_sess = mean(valid_current_area(:,:,rep_list(sess,:)),3,'omitnan'); % subset neuronal responses for current session
    
    subplot(2,2,sess,'units','normalized','position',[0.2+0.3*(sess-1) 0.65 0.275  0.26])
    internal_structure = corr(current_sess); % calculate the internal structure (correlation between time bins)
    internal_structure (boolean(eye(size(internal_structure ,1)))) = NaN; % convert main diagonal (all values of 1 by definition) into NaNs (for visuzalization)
    internal_structure(isnan(internal_structure)) = max(internal_structure(:),[],'omitnan'); % convert main diagonal to maximal correlation value (for visualiztion)
    imagesc(internal_structure) % visualize internal structure matrix
    if sess == 1
        title('Session 1')
        set(gca,'xtick',[])
        ylabel('Time bin')
    else
        title('Session 3')
        set(gca,'xtick',[],'ytick',[])
        set(gca,'ytick',[])
        cb = colorbar;
        set(cb,'position',[0.8 0.65 0.04 0.26])
        cb.Label.String = 'Correlation';
        cb.FontSize = 12;
    end
    
    % apply gaussian smoothing to single cells' neuronal responses
    smooth_current_rep = []; 
    for cell = 1:size(current_sess,1) % loop over cells
        smooth_current_rep(cell,:) = imgaussfilt(current_sess(cell,:),3); % sigma = 3
    end
    
    smooth_current_rep_norm = smooth_current_rep ./ max(smooth_current_rep,[],2); % normalize activity based on the time bin with maximal values
    if sess == 1
        [row,col] = find(smooth_current_rep_norm == 1);
        [B,I] = sort(col);
    end
    subplot(2,2,sess+2,'units','normalized','position',[0.2+0.3*(sess-1) 0.15 0.275  0.45])
    imagesc(smooth_current_rep_norm(row(I),:)) % visualize single cells'  neuronal responses
    
    if sess == 1
        ylabel('Sorted cell ID')
    else
        set(gca,'ytick',[])
        cb = colorbar;
        set(cb,'position',[0.8 0.15 0.04 0.45])
        cb.Label.String = 'Normalized activity';
        cb.FontSize = 12;
        
    end
    xlabel('Time bin')
    
end
colormap(newmap3)


%% Figure 7C - Stability of internal structure compared to population vector correlation

subset_population_vectors = {}; % define an empty variable that will store all mice recorded from each visual area
for area = 1:6 % loop over areas
    subset_population_vectors{area} = calcium_excitatory_population_vectors{area}(:,1); % subset all mice recorded from the current area
end

triu_id = boolean(triu(ones(30),1)); % define a boolean matrix in the size of 30x30 with true values in its upper half
num_shuffles = 1000; % define the number of pseudo-mice realizations

elapsed_session_all_measurments = {}; % define an empty variable that will store the pv corr and internal structure similarity values for all areas
for area = 1:6
    elapsed_session_pv = nan(num_shuffles,3); % define a NaN matrix that will store the pv corr values across pseudo-mice realizations
    elapsed_session_struc = nan(num_shuffles,3); % define a NaN matrix that will store the internal structure similarity values across pseudo-mice realizations
    
    current_area = cell2mat(subset_population_vectors{area}); % subset and pool the neuronal activity of all cells recorded from a single visual area
    
    active_sess1 = mean(mean(current_area(:,:,1:10),3,'omitnan'),2,'omitnan')>0; % active cells in session 1 (repeats 1-10)
    active_sess2 = mean(mean(current_area(:,:,11:20),3,'omitnan'),2,'omitnan')>0; % active cells in session 2 (repeats 11-20)
    active_sess3 = mean(mean(current_area(:,:,21:30),3,'omitnan'),2,'omitnan')>0; % active cells in session 3 (repeats 21-30)
    valid_cells = find([active_sess1 & active_sess2 & active_sess3]); % active in all three sessions
    
    valid_current_area = current_area(valid_cells,:,:); % subset only cells tat were active in all imaging sessions

    for shuffle = 1:num_shuffles % loop over realizations
        clc;
        disp(['Generating calcium imaging pseudo-mice and calculating'])
        disp(['internal structure and population vectors stability:'])
        disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Realization: ',num2str(shuffle),'/',num2str(num_shuffles)])
        
        
        valid_cells = sort(randperm(size(valid_current_area,1),round(size(valid_current_area,1)*0.7))); % subsample 70% of the pooled cell population
        subset_valid_current_area = valid_current_area(valid_cells,:,:); % included only cells that were chosen in "valid_cells"
        
        %calculate the average activity rate for session halves across movie repeats
            mean_activity_per_half = []; % define an empty variable that will store the averaged neuronal responses of each session half
            for half = 1:6 % loop over session halves
                current_half = [1:5] + 5*(half-1); % define range of movie repeats for current half
                mean_activity_per_half(:,:,half) = mean(subset_valid_current_area(:,:,current_half),3,'omitnan');  % average the neuronal activity across chosen movie repeats
            end
     
        % calculate the internal structure of each session half and the mean
        % pv correlation between session halves
        all_structures = []; % define an empty variable that will store the internal structure of each session half
        mean_pv = []; % define an empty variale that will store the mean pv correlation values between sessions halves
        for half1 = 1:6 % loop over halves
            current_half1 = mean_activity_per_half(:,:,half1); % subset the neruonal activity of a single session half
            for half2 = 1:6 % loop over halves
                current_half2 = mean_activity_per_half(:,:,half2); % subset the neruonal activity of a single session half
                
                current_structure = corr(current_half1,current_half2); % calculate the pv correlation between time bins of session halves
                 mean_pv(half1,half2) = mean(diag(current_structure),'omitnan'); % calculate the average pv corr value across corresponding time bins
                
                if half1==half2 % calculate internal structure for each half
                    all_structures(:,half1) = current_structure(triu_id); % store the internal structure of each half
                end
            end
        end

        structure_stability = corr(all_structures); % calculate the correlation across all the internal structures of all sessions halves
       
        % calculate the mean pv correlation as function of elapsed time
        within_session_pv_values = [mean_pv(1,2),mean_pv(3,4),mean_pv(5,6)]; % pv corr values between halves of the same session (within session)
            proximal_sessions_pv_values = [mean_pv(1:2,3:4),mean_pv(3:4,5:6)]; % pv corr values between halves of proximal sessions (sessions 1&2 and sessions 2&3)
            distal_session_pv_values = mean_pv(1:2,5:6); % pv corr values between halves of distal sessions (sessions 1&3)
 
        elapsed_session_pv(shuffle,1) = mean(within_session_pv_values(:),'omitnan'); % average pv corr values for within session
        elapsed_session_pv(shuffle,2) = mean(proximal_sessions_pv_values(:),'omitnan'); % average pv corr values for proximal sessions
        elapsed_session_pv(shuffle,3) = mean(distal_session_pv_values(:),'omitnan'); % average pv corr values for distal sessionsend
        
        % calculate the similarity beween internal structures as function of elapsed time
         within_session_struc_values = [structure_stability(1,2),structure_stability(3,4),structure_stability(5,6)]; % internal structure similarity values between halves of the same session (within session)
            proximal_sessions_struc_values = [structure_stability(1:2,3:4),structure_stability(3:4,5:6)]; % internal structure similarity values between halves of proximal sessions (sessions 1&2 and sessions 2&3)
            distal_session_struc_values = structure_stability(1:2,5:6); % internal structure similarity values between halves of distal sessions (sessions 1&3)
 
        elapsed_session_struc(shuffle,1) = mean(within_session_struc_values(:),'omitnan'); % average internal structure similarity values for within session
        elapsed_session_struc(shuffle,2) = mean(proximal_sessions_struc_values(:),'omitnan'); % average internal structure similarity corr values for proximal sessions
        elapsed_session_struc(shuffle,3) = mean(distal_session_struc_values(:),'omitnan'); % average internal structure similarity corr values for distal sessionsend
       
    end
    elapsed_session_all_measurments(area,:) = {elapsed_session_pv,elapsed_session_struc}; % store the pv corr and internal structure similarity values of each visual area
    
end

plt = []; % define an empty variable to store the plot information of each area
figure('units','normalized','position',[0.35 0.4 0.25 0.35]) % visualization of pv vs internal structure stability
for area = 1:6 % loop over areas
    current_area_pv = elapsed_session_all_measurments{area,1}; % subset pv corr values of a single area
    current_area_structure = elapsed_session_all_measurments{area,2}; % subset internal structure similarity values of a single area
    
    elapsed_session_struc_norm = current_area_structure./current_area_structure(:,1); % normalize internal structure similarity values based on 'within session' values
    elapsed_session_pv_norm = current_area_pv./current_area_pv(:,1); % normalize pv corr values based on 'within session' values
    
    mean_struct = mean(elapsed_session_struc_norm,'omitnan'); % calculate mean internal structure similarity across pseudo-mice realizations
    std_struct = std(elapsed_session_struc_norm,'omitnan'); % calculate standard diviation across pseudo-mice realizations
    
    mean_pv_corr = mean(elapsed_session_pv_norm,'omitnan'); % calculate mean pv corr across pseudo-mice realizations
    std_pv_corr  = std(elapsed_session_pv_norm,'omitnan'); % calculate standard diviation across pseudo-mice realizations
        
    hold on
    plt(area) = errorbar(mean_pv_corr ,std_pv_corr ,'o','color',[0.2 0.2 0.2]+0.1*(area-1),...
        'markerfacecolor',[0.2 0.2 0.2]+0.1*(area-1),'capsize',0,'linestyle','-','linewidth',3,'markersize',3);
    
    plt(area+6) = errorbar(mean_struct,std_struct,'o','color',colors(area,:),...
        'markerfacecolor',colors(area,:),'capsize',0,'linestyle','-','linewidth',3,'markersize',3);
    
end
xlim([0.5 3.5])
ylim([0.7 1.01])
set(gca,'xtick',1:3,'xticklabel',{'Within session','Proximal sessions','Distal sessions'})
ylabel('Correlation (normalized)')
legend(plt,brain_areas(1:6),'Location','best')
legend('boxoff')


%% Figure 7D - Internal structure stability as a function of the number of cells

subset_population_vectors = {}; % define an empty variable that will store all mice recorded from each visual area
for area = 1:6 % loop over areas
    subset_population_vectors{area} = calcium_excitatory_population_vectors{area}(:,1); % subset all mice recorded from the current area
end

triu_id = boolean(triu(ones(30),1));  % define a boolean matrix in the size of 30x30 with true values in its upper half
cell_count_list_full = [25, 50, 100, 250,500,1000,2000,4000,6000, 7900]; % cell count threshold list
elapsed_session_all_measurments = {}; % define an empty variable that will store the pv corr and internal structure similarity values for all areas
for area = 1:6 % loop over areas
    current_area = cell2mat(subset_population_vectors{area}); % subset and pool the neuronal activity of all cells recorded from a single visual area
    
     active_sess1 = mean(mean(current_area(:,:,1:10),3,'omitnan'),2,'omitnan')>0; % active cells in session 1 (repeats 1-10)
    active_sess2 = mean(mean(current_area(:,:,11:20),3,'omitnan'),2,'omitnan')>0; % active cells in session 2 (repeats 11-20)
    active_sess3 = mean(mean(current_area(:,:,21:30),3,'omitnan'),2,'omitnan')>0; % active cells in session 3 (repeats 21-30)
    valid_cells = find([active_sess1 & active_sess2 & active_sess3]); % active in all three sessions
    
    valid_current_area = current_area(valid_cells,:,:); % subset only cells tat were active in all imaging sessions

    cell_count_list = cell_count_list_full(cell_count_list_full<=size(valid_current_area,1)); % subsample cell count list based on the total number of cells recorded in the current area
    
    shuffle_num = 100; % define the number of pseudo-mice realizations
    elapsed_session_pv = nan(length(cell_count_list),3,shuffle_num); % define a NaN matrix that will store the pv corr values across pseudo-mice realizations
    elapsed_session_struc = nan(length(cell_count_list),3,shuffle_num); % define a NaN matrix that will store the internal structure similarity values across pseudo-mice realizations
    
    for cell_count = 1:length(cell_count_list) % loop over cell count treshold
        for shuffle = 1:shuffle_num % loop over realizations
            clc;
            disp(['Generating calcium imaging pseudo-mice and calculating'])
            disp(['internal structure and population vectors stability:'])
            disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Cells included: ',num2str(cell_count_list_full(cell_count)),'/',num2str(cell_count_list(end)),...
                ' | Realization: ',num2str(shuffle),'/',num2str(num_shuffles)])
            
            
            valid_cells = sort(randperm(size(valid_current_area,1),cell_count_list(cell_count))); % subsample #cell_count_list(cell_count) of cells
            subset_valid_current_area = valid_current_area(valid_cells,:,:); % included only cells that were chosen in "valid_cells"

            %calculate the average activity rate for session halves across movie repeats
            mean_activity_per_half = []; % define an empty variable that will store the averaged neuronal responses of each session half
            for half = 1:6 % loop over session halves
                current_half = [1:5] + 5*(half-1); % define range of movie repeats for current half
                mean_activity_per_half(:,:,half) = mean(subset_valid_current_area(:,:,current_half),3,'omitnan');  % average the neuronal activity across chosen movie repeats
            end
     
           % calculate the internal structure of each session half and the mean
        % pv correlation between session halves
        all_structures = []; % define an empty variable that will store the internal structure of each session half
        mean_pv = []; % define an empty variale that will store the mean pv correlation values between sessions halves
        for half1 = 1:6 % loop over halves
            current_half1 = mean_activity_per_half(:,:,half1); % subset the neruonal activity of a single session half
            for half2 = 1:6 % loop over halves
                current_half2 = mean_activity_per_half(:,:,half2); % subset the neruonal activity of a single session half
                
                current_structure = corr(current_half1,current_half2); % calculate the pv correlation between time bins of session halves
                 mean_pv(half1,half2) = mean(diag(current_structure),'omitnan'); % calculate the average pv corr value across corresponding time bins
                
                if half1==half2 % calculate internal structure for each half
                    all_structures(:,half1) = current_structure(triu_id); % store the internal structure of each half
                end
            end
        end

        structure_stability = corr(all_structures); % calculate the correlation across all the internal structures of all sessions halves
       
            % calculate the mean pv correlation as function of elapsed time
        within_session_pv_values = [mean_pv(1,2),mean_pv(3,4),mean_pv(5,6)]; % pv corr values between halves of the same session (within session)
            proximal_sessions_pv_values = [mean_pv(1:2,3:4),mean_pv(3:4,5:6)]; % pv corr values between halves of proximal sessions (sessions 1&2 and sessions 2&3)
            distal_session_pv_values = mean_pv(1:2,5:6); % pv corr values between halves of distal sessions (sessions 1&3)
 
        elapsed_session_pv(cell_count,1,shuffle) = mean(within_session_pv_values(:),'omitnan'); % average pv corr values for within session
        elapsed_session_pv(cell_count,2,shuffle) = mean(proximal_sessions_pv_values(:),'omitnan'); % average pv corr values for proximal sessions
        elapsed_session_pv(cell_count,3,shuffle) = mean(distal_session_pv_values(:),'omitnan'); % average pv corr values for distal sessionsend
        
        % calculate the similarity beween internal structures as function of elapsed time
         within_session_struc_values = [structure_stability(1,2),structure_stability(3,4),structure_stability(5,6)]; % internal structure similarity values between halves of the same session (within session)
            proximal_sessions_struc_values = [structure_stability(1:2,3:4),structure_stability(3:4,5:6)]; % internal structure similarity values between halves of proximal sessions (sessions 1&2 and sessions 2&3)
            distal_session_struc_values = structure_stability(1:2,5:6); % internal structure similarity values between halves of distal sessions (sessions 1&3)
 
        elapsed_session_struc(cell_count,1,shuffle) = mean(within_session_struc_values(:),'omitnan'); % average internal structure similarity values for within session
        elapsed_session_struc(cell_count,2,shuffle) = mean(proximal_sessions_struc_values(:),'omitnan'); % average internal structure similarity corr values for proximal sessions
        elapsed_session_struc(cell_count,3,shuffle) = mean(distal_session_struc_values(:),'omitnan'); % average internal structure similarity corr values for distal sessionsend
        end
    end
    elapsed_session_all_measurments(area,:) = {elapsed_session_pv,elapsed_session_struc}; % store the pv corr and internal structure similarity values of each visual area
    
end


figure('units','normalized','position',[0.35 0.4 0.45 0.45]) % visualization of pv vs internal structure stability
for area =  1:6 % loop over areas
    curent_area_pv = elapsed_session_all_measurments{area,1}; % subset pv corr values of a single area
    curent_area_pv = mean(curent_area_pv,3,'omitnan'); % average pv corr values across realizations
    elapsed_session_pv_norm = [curent_area_pv./curent_area_pv(:,1,:)]; % normalize pv corr values based on 'within session' values
    mean_elapsed_session_pv_norm = nanmean(elapsed_session_pv_norm,3);
    
    curent_area_struct = elapsed_session_all_measurments{area,2}; % calculate mean internal structure similarity across pseudo-mice realizations
    curent_area_struct = mean(curent_area_struct,3,'omitnan'); % average internal structure similarity values across realizations
    elapsed_session_struct_norm = [curent_area_struct./curent_area_struct(:,1,:)];  % normalize internal structure similarity values based on 'within session' values
    mean_elapsed_session_struct_norm = mean(elapsed_session_struct_norm,3);

    
    
    subplot(2,3,area)
    for cell_count = 1:size(mean_elapsed_session_struct_norm,1)
        hold on
        plot(mean_elapsed_session_pv_norm(cell_count,:),'color',[0.2 0.2 0.2] + 0.075*(cell_count-1),'linewidth',3)
        plot(mean_elapsed_session_struct_norm(cell_count,:),'color',newmap3(end-6*(cell_count-1),:),'linewidth',3)
        
    end
    text(0.1,0.175,brain_areas{area},'Units','normalized','FontSize',15)
    
    ylim([0.7 1.01])
    xlim([0.75 3.25])
    if area == 3
        cb = colorbar('north','xtick',[0 1],'xticklabel',[25 7900]);
        cb.Position = [0.7 0.9 0.2 0.05];
        colormap(flipud(newmap3))
    end
    if area >3
        set(gca,'xtick',1:3,'xticklabel',{'Within session','Proximal sessions','Distal sessions'})
        xtickangle(15)
    else
        set(gca,'xtick',[])
    end
    if area ==4 || area == 1
        ylabel('Correlation (normalized)')
    else
        set(gca,'ytick',[])
    end
    
end
suptitle('Number of cells included:')

%% Figure S1D - Single unit example of within and between block rate and tuning stability

nat_movie = 1; % natural movie 1
mouse = 23; % example mouse #23
area = 1; % area V1
unit = 9; % example unit #9
example_mouse = neuropixels_population_vectors{mouse,area,nat_movie}; % subset the neuronal activity of example mouse
current_unit = squeeze(example_mouse(unit,:,:))'.*30; % subset neuronal activity of example unit

figure('units','normalized','position',[0.575 0.4 0.1 0.2]) % visualize tuning curve of example unit
hold on
mean_current_unit = [];
mean_current_unit(1,:) = mean(current_unit(1:10,:),'omitnan'); % average tuning curve of block A
mean_current_unit(2,:) = mean(current_unit(11:20,:),'omitnan'); % average tuning curve of block B
plot(mean_current_unit(1,:),'-','markersize',20,'linewidth',2,'color',[0.2 0.2 0.2])
plot(mean_current_unit(2,:)','-','markersize',20,'linewidth',2,'color',colors(area,:))
[r,p] = corr(mean_current_unit(1,:)',mean_current_unit(2,:)'); % tuning curve similarity across blocks
text(0.325, 0.675,['r: ',num2str(r)],'Units','normalized','color',[0 0 0])
legend({'Block A','Block B'},'Location','best')
legend('boxoff')
xlim([1 30])
ylabel('Mean activity rate (HZ)')
xlabel('Time in movie (sec)')


current_unit = squeeze(example_mouse(unit,:,11:20))*30'; % subset neuronal activity of example unit during block B 
figure('units','normalized','position',[0.3 0.3 0.25 0.425]) % visualize tuning curve of example unit across movie repeats
for repeat = 1:10 % loop over movie repeats
    hold on
    plot3(ones(1,30)*repeat,[1:30],current_unit(:,repeat),'linewidth',4,'color',magma_colormap(1+25*(repeat-1),:))
end
set(gca,'ytick',5:5:30,'xtick',1:10)
ylabel('Time in movie (sec)')
xlabel('Movie repeat')
zlabel('Activity rate (Hz)')
view(138.8,45.2)
title('Neuropixels - example unit:')

%% Figure S1E - Dissociation between tuning curve stability and activity rate stability for three V1 units

cell_cutoff = 15; % minimal number of cells recorded
nat_movie = 1; % natual movie 1
repeats_num = 30; % number of movie repeats in each stimuli block
area = 1; % area V1
current_area = neuropixels_population_vectors(:,area,nat_movie); % subset neuronal activity of all mice recorded from V1
valid_mice = neuropixels_cell_count(:,area,nat_movie) >= cell_cutoff & movie_repeats(:,nat_movie) == repeats_num; % subset only mice from the 'Functional connectivity' group
all_cells = cell2mat(current_area(valid_mice))*30; % pool all V1 cells across mice

% calculate the average activity rate and tuning curve across blocks for
% each unit.
activity_rate = [];
tuning_curve_across_blocks = [];
for unit = 1:size(all_cells,1) % loop over units
    current_cellA = squeeze(all_cells(unit,:,1:30)); % neuronal activity during block A (repeats 1-30)
    current_cellB = squeeze(all_cells(unit,:,31:60));% neuronal activity during block B (repeats 31-60)
    
    activity_rate(unit,:) = [mean(mean(current_cellA,'omitnan')),mean(mean(current_cellB,'omitnan'))]; % average activity rate across blocks
    tuning_curve_across_blocks(unit) = corr(mean(current_cellA,2,'omitnan'),mean(current_cellB,2,'omitnan')); % tuning curve corr across blocks
end

text_pos_list = [0.1, 0.95; 0.5, 0.95; 0.5, 0.95];
cell_list = [216,1650,72];
figure('units','normalized','position',[0.25 0.4 0.45 0.25]) % visulize ativity patterns across blocks for 3 example units
for unit = 1:length(cell_list)
    
    current_cellA = nanmean(squeeze(all_cells(cell_list(unit),:,1:30)),2);
    current_cellB = nanmean(squeeze(all_cells(cell_list(unit),:,31:60)),2);
    subplot(1,3,unit)
    
    hold on
    plot(current_cellA,'color',[0.2 0.2 0.2],'linewidth',2)
    plot(current_cellB,'color',colors(1,:),'linewidth',2)
    
    text(text_pos_list(unit,1), text_pos_list(unit,2),['A: ',num2str(nanmean(current_cellA)),' Hz'],'Units','normalized','color',[0.2 0.2 0.2])
    text(text_pos_list(unit,1), text_pos_list(unit,2)-0.05,['B: ',num2str(nanmean(current_cellB)),' Hz'],'Units','normalized','color',colors(1,:))
    text(text_pos_list(unit,1), text_pos_list(unit,2)-0.1,['r: ',num2str(tuning_curve_across_blocks(cell_list(unit)))],'Units','normalized','color',[0 0 0])
    
    if unit ==1
        ylabel('Activity rate (spike/sec)')
        title({'Stable tuning';'stable firing rate'})
    elseif unit == 2
        xlabel('Time in movie (sec)')
        title({'Stable tuning';'unstable firing rate'})
    else
        title({'Unstable tuning';'stable firing rate'})
    end
    xlim([1 30])
end


%% Figure S1F - Relationship between single cell tuning curve stability and single cell activity rate - all V1 units

cell_cutoff = 15; % minimal number of cells recorded
nat_movie = 1; % natual movie 1
repeats_num = 30; % number of movie repeats in each stimuli block
area = 1; % area V1
current_area = neuropixels_population_vectors(:,area,nat_movie); % subset neuronal activity of all mice recorded from V1
valid_mice = neuropixels_cell_count(:,area,nat_movie) >= cell_cutoff & movie_repeats(:,nat_movie) == repeats_num; % subset only mice from the 'Functional connectivity' group
all_cells = cell2mat(current_area(valid_mice))*30; % pool all V1 cells across mice

% calculate the average activity rate and tuning curve across blocks for
% each unit.
activity_rate = [];
tuning_curve_across_blocks = [];
for unit = 1:size(all_cells,1) % loop over units
    current_cellA = squeeze(all_cells(unit,:,1:30)); % neuronal activity during block A (repeats 1-30)
    current_cellB = squeeze(all_cells(unit,:,31:60));% neuronal activity during block B (repeats 31-60)
    
    activity_rate(unit,:) = [mean(mean(current_cellA,'omitnan')),mean(mean(current_cellB,'omitnan'))]; % average activity rate across blocks
    tuning_curve_across_blocks(unit) = corr(mean(current_cellA,2,'omitnan'),mean(current_cellB,2,'omitnan')); % tuning curve corr across blocks
end
r  = corr(tuning_curve_across_blocks.^2',nanmean(activity_rate,2)); % explained variance between average activity rate and tuning curve stability across blocks


figure('units','normalized','position',[0.35 0.4 0.25 0.35]) % visulize the relationship between activity rate and tuning curve stability 
hold on
dscatter(tuning_curve_across_blocks',nanmean(activity_rate,2))
colormap(magma_colormap)
xlabel({'Tuning curve correlation'; 'between Natural movie 1 blocks'})
ylabel('Mean firing rate (spike/sec)')
cb = colorbar;
set(cb,'xtick',[0.01 1],'xticklabel',{'min','max'})
cb.Label.String = 'Density';
cb.FontSize = 12;

text(0.1,0.95,['N = ',num2str(length(tuning_curve_across_blocks))],'Units','normalized','color',[0 0 0])
text(0.1,0.9,['R^2 = ',num2str(r)],'Units','normalized','color',[0 0 0])

title({'Neuropixels - all V1 units from the';"'Function Connectivity' group:"})

%% Figure S1G - Relationship between PV, ensemble rate and tuning curve correlations - repeat resolution
nat_movie = 1; % natural movie 1
cell_cutoff = 15; % minimal number of recorded cells
num_repeats = 30; % number of movie repeats in each stimulus block
for area = 1:6 % loop over 
    valid_mice = [neuropixels_cell_count(:,area,nat_movie) >= cell_cutoff] & [movie_repeats(:,nat_movie) == num_repeats]; % include only mice from the function connectivity group with more than 15 cells
    current_area = neuropixels_population_vectors(valid_mice,area,nat_movie); % subset all mice that passed the requirments
    
    % calculate the pv, rate and tuning corr across pairs of movie repeats
    % andmeasure the linear relationship between the different measurments
    rate_tuning_relationship = [];
    rate_pv_relationship = [];
    tuning_pv_relationship = [];
    pv_rate_tuning_relationship = [];
    for mouse = 1:size(current_area,1) % loop over mice
        
        clc;
        disp(['Calculating the relationship between population vectors,'])
        disp(['ensemble rate and tuning curve correlations'])
        disp(['across movie repeats within a block:'])
        disp(['Dataset: Neuropixels | Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse}; % subset neuronal activity of a single mouse (#units x 30 time bins x 60 movie repeats)
        
        rate_corr_blockA = corr(squeeze(mean(current_mouse(:,:,1:30),2,'omitnan'))); % calculate rate correlation across pairs of repeats in block A
        rate_corr_blockB = corr(squeeze(mean(current_mouse(:,:,31:60),2,'omitnan'))); % calculate rate correlation across pairs of repeats in block B
        tuning_corr_blockA = []; 
        tuning_corr_blockB = [];
        pv_corr_blockA = [];
        pv_corr_blockB = [];
        for repeat1 = 1:30 % loop over movie repeats
            for repeat2 = 1:30 % loop over movie repeats
                tuning_corr_blockA(repeat1,repeat2) = median(diag(corr(current_mouse(:,:,repeat1)',current_mouse(:,:,repeat2)')),'omitnan'); % calculate median tuning curve correlation across pair of movie repeats of block A
                tuning_corr_blockB(repeat1,repeat2) = median(diag(corr(current_mouse(:,:,repeat1+num_repeats)',current_mouse(:,:,repeat2+num_repeats)'))); % calculate median tuning curve correlation across pair of movie repeats of block B
                
                pv_corr_blockA(repeat1,repeat2) = mean(diag(corr(current_mouse(:,:,repeat1),current_mouse(:,:,repeat2))),'omitnan'); % calculate the pv corr between movie repeats of block A and average across corresponding time bins 
                pv_corr_blockB(repeat1,repeat2) = mean(diag(corr(current_mouse(:,:,repeat1+num_repeats),current_mouse(:,:,repeat2+num_repeats))),'omitnan'); % calculate the pv corr between movie repeats of block B and average across corresponding time bins 
            end
        end
        
        % calculate the linear relationship between rate,pv and tuning
        % while controling for the time interval between movie repeats
        for diagonal = 1:20 % loop over time intervals
            mdl = fitlm(diag(rate_corr_blockA,diagonal),diag(tuning_corr_blockA,diagonal));
            rate_tuning_relationship(mouse,diagonal,1) = mdl.Rsquared.Ordinary; % explained variance between rate and tuning corr block A
            
            mdl = fitlm(diag(rate_corr_blockB,diagonal),diag(tuning_corr_blockB,diagonal));
            rate_tuning_relationship(mouse,diagonal,2) = mdl.Rsquared.Ordinary; % explained variance between rate and tuning corr block B
            
            mdl = fitlm(diag(rate_corr_blockA,diagonal),diag(pv_corr_blockA,diagonal));
            rate_pv_relationship(mouse,diagonal,1) = mdl.Rsquared.Ordinary; % explained variance between rate and pv corr block A
            
            mdl = fitlm(diag(rate_corr_blockB,diagonal),diag(pv_corr_blockB,diagonal));
            rate_pv_relationship(mouse,diagonal,2) = mdl.Rsquared.Ordinary; % explained variance between rate and pv corr block B
            
            mdl = fitlm(diag(tuning_corr_blockA,diagonal),diag(pv_corr_blockA,diagonal));
            tuning_pv_relationship(mouse,diagonal,1) = mdl.Rsquared.Ordinary; % explained variance between pv and tuning corr block A
            
            mdl = fitlm(diag(tuning_corr_blockB,diagonal),diag(pv_corr_blockB,diagonal));
            tuning_pv_relationship(mouse,diagonal,2) = mdl.Rsquared.Ordinary; % explained variance between pv and tuning corr block B
            
            predictors_blockA = [diag(rate_corr_blockA,diagonal),diag(tuning_corr_blockA,diagonal)];
            mdl = fitlm(predictors_blockA,diag(pv_corr_blockA,diagonal));
            pv_rate_tuning_relationship(mouse,diagonal,1) = mdl.Rsquared.Ordinary; % explained variance between pv, rate and tuning corr block A
            
            predictors_blockB = [diag(rate_corr_blockB,diagonal),diag(tuning_corr_blockB,diagonal)];
            mdl = fitlm(predictors_blockB,diag(pv_corr_blockB,diagonal));
            pv_rate_tuning_relationship(mouse,diagonal,2) = mdl.Rsquared.Ordinary; % explained variance between pv, rate and tuning corr block B
            
        end
    end
    
    % average explined variance across time intervals and blocks
    rate_pv_explained = mean(mean(rate_pv_relationship,2,'omitnan'),3,'omitnan');
    tuning_pv_explained = mean(mean(tuning_pv_relationship,2,'omitnan'),3,'omitnan');
    rate_tuning_explained = mean(mean(rate_tuning_relationship,2,'omitnan'),3,'omitnan');
    pv_rate_tuning_explained = mean(mean(pv_rate_tuning_relationship,2,'omitnan'),3,'omitnan');
    
    all_models{area} = [pv_rate_tuning_explained(:),rate_pv_explained(:),tuning_pv_explained(:),rate_tuning_explained(:)];
end

pos_list = [0.1 0.575 0.25 0.325;0.4 0.575 0.25 0.325;0.7 0.575 0.25 0.325;...
    0.1 0.2 0.25 0.325;0.4 0.2 0.25 0.325;0.7 0.2 0.25 0.325];
figure('units','normalized','position',[0.25 0.2 0.45 0.55]) % visualize the linear relationships between the different measurments
for area = 1:6
    current_area = all_models{area};
    subplot(2,3,area,'units','normalized','position',pos_list(area,:))
    
    figure_boxplot(current_area);
    text(0.775,0.935,brain_areas{area},'Units','normalized','FontSize',15)
    text(0.75,0.85,['N=',num2str(size(current_area,1))],'Units','normalized','FontSize',13)
    
    set(gca,'ActivePositionProperty','position')
    if area == 4 || area == 1
        ylabel('Explained variance (R^2)')
    end
    if area >3
        set(gca,'xticklabels',{'PV ~ Rate & Tuning','PV ~ Rate','PV ~ Tuning','Rate ~ Tuning'})
        xtickangle(30)
    else
        set(gca,'xtick',[])
    end
    ylim([-0.1 1.1])
end
suptitle('Neuropixels - across movie repeats within a block:')

%% Figure S1H - Relationship between PV, ensemble rate and tuning curve correlations - days resolution

cell_cutoff = 20; % minimal number of cells recorded in each FOV
for area = 1:6 % loop over areas
    valid_mice = min(calcium_excitatory_cell_count{area},[],2) >= cell_cutoff; % include only mice with atleast 20 cells in each of the three sessions
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,1); % subset only mice that passed the cell count requirments
    pv_rate_tuning_relationships = [];
    for mouse = 1:length(current_area) % loop over mice
        clc;
        disp(['Calculating the relationship between population vectors,'])
        disp(['ensemble rate and tuning curve correlationsacross sessions:'])
        disp(['Dataset: Calcium imaging | Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse}; % subset neuronal activity of a single mouse (#cells x 30 time bins x 30 movie repeats)
        
        % calculate the pv, rate and tuning curve corr across sessions
        % halves using only the cells that were active in both time point
        rate_corr = [];
        pv_corr = [];
        tuning_corr = [];
        for sess1 = 1:6 % loop over halves
            current_sess1 = mean(current_mouse(:,:,[1:5]+5*(sess1-1)),3,'omitnan'); % average activity over movie repeats
            for sess2 = 1:6 % loop over halves
                current_sess2 = mean(current_mouse(:,:,[1:5]+5*(sess2-1)),3,'omitnan');% average activity over movie repeats
                
                valid_cells = mean(current_sess1,2) > 0 & mean(current_sess2,2) > 0; % only cells that were active in both session halves
                
                rate_corr(sess1,sess2) = corr(mean(current_sess1(valid_cells,:),2,'omitnan'),mean(current_sess2(valid_cells,:),2,'omitnan'),'rows','pairwise');
                pv_corr(sess1,sess2) = mean(diag(corr(current_sess1(valid_cells,:),current_sess2(valid_cells,:),'rows','pairwise')),'omitnan');
                tuning_corr(sess1,sess2) = median(diag(corr(current_sess1(valid_cells,:)',current_sess2(valid_cells,:)','rows','pairwise')),'omitnan');
            end
        end
        
        % subset only values off the main diagonal
        rate_temp = [rate_corr(1:2,3:6),rate_corr(3:4,5:6)]; 
        pv_temp = [pv_corr(1:2,3:6),pv_corr(3:4,5:6)];
        tuning_temp = [tuning_corr(1:2,3:6),tuning_corr(3:4,5:6)];
        
        % calculate the explained variance between the different measurments
        mdl = fitlm([rate_temp(:),tuning_temp(:)],pv_temp(:)); % fit model to explain pv using rate and tuning
         pv_rate_tuning_relationships(mouse,1) = mdl.Rsquared.Ordinary;
        pv_rate_tuning_relationships(mouse,2) = corr(rate_temp(:),pv_temp(:),'rows','pairwise').^2; % explained variance of pv using rate corr
        pv_rate_tuning_relationships(mouse,3) = corr(pv_temp(:),tuning_temp(:),'rows','pairwise').^2; % explained variance of pv using tuning corr
        pv_rate_tuning_relationships(mouse,4) = corr(rate_temp(:),tuning_temp(:),'rows','pairwise').^2; % explained variance of tuning using rate corr
        
        
    end
    pv_rate_tuning_relationships_areas{area} = pv_rate_tuning_relationships;
    
end

% visualize the distribution of explained variance between the different measurments
pos_list = [0.1 0.575 0.25 0.325;0.4 0.575 0.25 0.325;0.7 0.575 0.25 0.325;...
    0.1 0.2 0.25 0.325;0.4 0.2 0.25 0.325;0.7 0.2 0.25 0.325];
figure('units','normalized','position',[0.25 0.2 0.45 0.55])
for area = 1:6
    current_area = pv_rate_tuning_relationships_areas{area};
    subplot(2,3,area,'units','normalized','position',pos_list(area,:))
    
    figure_boxplot(current_area);
    text(0.775,0.935,brain_areas{area},'Units','normalized','FontSize',15)
    text(0.75,0.85,['N=',num2str(size(current_area,1))],'Units','normalized','FontSize',13)
    
    set(gca,'ActivePositionProperty','position')
    if area == 4 || area == 1
        ylabel('Explained variance (R^2)')
    end
    if area >3
        set(gca,'xticklabels',{'PV ~ Rate & Tuning','PV ~ Rate','PV ~ Tuning','Rate ~ Tuning'})
        xtickangle(30)
    else
        set(gca,'xtick',[])
    end
    ylim([-0.1 1.1])
end
suptitle('Calcium imaging - across days:')

%% Figure S1I - relationship between activity rate, activity diff and tuning - repeat resolution
nat_movie = 1; % natural movie 1
cell_cutoff = 15; % minimal number of cells recorded critria
num_repeats = 30; % number of movie repeats in each stimuli block
for area = 1:6 % loop over areas
    valid_mice = neuropixels_cell_count(:,area,nat_movie) >= cell_cutoff & movie_repeats(:,nat_movie) == num_repeats; % include only mice from the functional connectivity group with at least 15 recorded cells
    current_area = neuropixels_population_vectors(valid_mice,area,nat_movie); % subset only mice that passed the requirements of 'valid_mice'
    
    all_mice_relationships = []; % loop over mice
    for mouse = 1:size(current_area,1)
        
        clc;
        disp(['Calculating the relationship between activity rate,'])
        disp(['absolute activity rate difference or absolute activity rate difference score'])
        disp(['and tuning curve correlation across movie repeats within a block:'])
        disp(['Dataset: Neuropixels | Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])

        current_mouse = current_area{mouse}; % subset neuronal activity of a single mouse (#units x 30 time bins x 60 movie repeats)
        
        for repeat1 = 1:30 % loop over movie repeats
            for repeat2 = 1:30 % loop over movie repeats
                tuning_corr_blockA = diag(corr(current_mouse(:,:,repeat1)',current_mouse(:,:,repeat2)')); % tuning curve corr between pair of movie repeats of block A across corresponding units
                tuning_corr_blockB = diag(corr(current_mouse(:,:,repeat1+num_repeats)',current_mouse(:,:,repeat2+num_repeats)')); % tuning curve corr between pair of movie repeats of block B across corresponding units
                
                activity_rate_blockA = mean(mean(current_mouse(:,:,[repeat1,repeat2]),3,'omitnan'),2,'omitnan'); % calculate the average activity rate across pair of movie repeats of block A
                activity_rate_blockB = mean(mean(current_mouse(:,:,[repeat1+num_repeats,repeat2+num_repeats]),3,'omitnan'),2,'omitnan'); % calculate the average activity rate across pair of movie repeats of block B
                
                activity_diff_blockA = abs(mean(current_mouse(:,:,repeat1),2,'omitnan') - mean(current_mouse(:,:,repeat2),2,'omitnan')); % calculate the absolute activity rate difference between pair of movie repeats in block A
                activity_diff_blockB =  abs(mean(current_mouse(:,:,repeat1+num_repeats),2,'omitnan') - mean(current_mouse(:,:,repeat2+num_repeats),2,'omitnan')); % calculate the absolute activity rate difference between pair of movie repeats in block A

                activity_diff_score_blockA = abs([mean(current_mouse(:,:,repeat1),2,'omitnan') - mean(current_mouse(:,:,repeat2),2,'omitnan')]./...
                    [mean(current_mouse(:,:,repeat1),2,'omitnan') + mean(current_mouse(:,:,repeat2),2,'omitnan')]); % calculate the absolute activity rate difference score between pair of movie repeats in block A
                activity_diff_score_blockB = abs([mean(current_mouse(:,:,repeat1+num_repeats),2,'omitnan') - mean(current_mouse(:,:,repeat2+num_repeats),2,'omitnan')]./...
                    [mean(current_mouse(:,:,repeat1+num_repeats),2,'omitnan') + mean(current_mouse(:,:,repeat2+num_repeats),2,'omitnan')]); % calculate the absolute activity rate difference score between pair of movie repeats in block B
                
                % linear relationship between tuning curve stability and activity rate
                tuning_VS_activity_rate_blockA(repeat1,repeat2) = corr(activity_rate_blockA,tuning_corr_blockA,'rows','pairwise'); 
                tuning_VS_activity_rate_blockB(repeat1,repeat2) = corr(activity_rate_blockB,tuning_corr_blockB,'rows','pairwise'); 
                
                % linear relationship between tuning curve stability and
                % absolute activity rate diffence
                tuning_VS_activity_diff_blockA(repeat1,repeat2) = corr(activity_diff_blockA,tuning_corr_blockA,'rows','pairwise');
                tuning_VS_activity_diff_blockB(repeat1,repeat2) = corr(activity_diff_blockB,tuning_corr_blockB,'rows','pairwise');
                
                 % linear relationship between tuning curve stability and
                % absolute activity rate diffence score
                tuning_VS_activity_diff_score_blockA(repeat1,repeat2) = corr(activity_diff_score_blockA,tuning_corr_blockA,'rows','pairwise');
                tuning_VS_activity_diff_score_blockB(repeat1,repeat2) = corr(activity_diff_score_blockB,tuning_corr_blockB,'rows','pairwise');
                
            end
        end

        triu_ind = boolean(triu(ones(30),1)); % boolean matrix to be used to extract the upper half values of correlation matrices
        
        current_mouse_relationships = [[tuning_VS_activity_rate_blockA(triu_ind);tuning_VS_activity_rate_blockB(triu_ind)],...
            [tuning_VS_activity_diff_blockA(triu_ind);tuning_VS_activity_diff_blockB(triu_ind)],...
            [tuning_VS_activity_diff_score_blockA(triu_ind);tuning_VS_activity_diff_score_blockB(triu_ind)]].^2; % convert correlations into explained varince
        
        % average the explained variance values across pairs of movie repeats
        % resulting in a single value for each mouse
        all_mice_relationships = [all_mice_relationships;mean(current_mouse_relationships,'omitnan')]; 

    end    
    all_mice_relationships_areas{area} = all_mice_relationships;
end

% visulalize the distribution of explained variance values of the different
% models across visual brain areas
pos_list = [0.1 0.575 0.25 0.325;0.4 0.575 0.25 0.325;0.7 0.575 0.25 0.325;...
    0.1 0.2 0.25 0.325;0.4 0.2 0.25 0.325;0.7 0.2 0.25 0.325];
figure('units','normalized','position',[0.25 0.2 0.45 0.55])
for area = 1:6
    current_area = all_mice_relationships_areas{area};
    subplot(2,3,area,'units','normalized','position',pos_list(area,:))
    
    figure_boxplot(current_area);
    text(0.775,0.935,brain_areas{area},'Units','normalized','FontSize',15)
    text(0.75,0.85,['N=',num2str(size(current_area,1))],'Units','normalized','FontSize',13)
    
    set(gca,'ActivePositionProperty','position')
    if area == 4 || area == 1
        ylabel({'Explained variance by the';'tuning curve correlation (R^2)'})
    end
    if area >3
        set(gca,'xticklabels',{'Mean activity rate','Absolute activity diff.','Absolute activity diff. score'})
        
        xtickangle(30)
    else
        set(gca,'xtick',[])
    end
    ylim([-0.1 1.1])
end
suptitle('Neuropixels - across movie repeats within a block:')


%% Figure S1J - Relationship between activity rate, activity diff and tuning - days resolution
cell_cutoff = 20; % minimal number of cells recorded in a given session
for area = 1:6 % loop over areas
    valid_mice = min(calcium_excitatory_cell_count{area},[],2) >= cell_cutoff; % include only mice with atleast 20 recorded cells in each of the three sessions
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,1); % subset mice that passed the requirments of 'valid_mice'
    
    
    all_relationships = [];
    for mouse = 1:length(current_area) % loop over mice
        clc;
        disp(['Calculating the relationship between activity rate,'])
        disp(['absolute activity rate difference or absolute activity rate '])
        disp(['difference score and tuning curve correlation across days:'])
        disp(['Dataset: Calcium imaging | Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse}; % subset neuronal activity of a single mouse (#cells x 30 time bins x 30 movie repeats)
        
        % calculate the linear relationship between tuning curve stability
        % and activity rate measurments
        tuning_vs_activity_rate = [];
        tuning_rate_vs_activity_diff = [];
        tuning_rate_vs_activity_diff = [];
        for sess1 = 1:6 % loop over session halves
            current_sess1 = mean(current_mouse(:,:,[1:5]+5*(sess1-1)),3,'omitnan'); % average neuronal activity over movie repeats
            for sess2 = 1:6 % loop over session halves
                current_sess2 = mean(current_mouse(:,:,[1:5]+5*(sess2-1)),3,'omitnan'); % average neuronal activity over movie repeats
                valid_cells = mean(current_sess1,2,'omitnan') > 0 & mean(current_sess2,2,'omitnan') > 0; % include only cells that were active in both compared time points
                
                activity_rate_both_sess = [mean(current_sess1(valid_cells,:),2,'omitnan'),mean(current_sess2(valid_cells,:),2,'omitnan')]; % averge activity across time bins for each half
                mean_activity_rate = mean(activity_rate_both_sess,2,'omitnan'); % average activity across sesssions halves
                abs_activity_diff = abs(activity_rate_both_sess(:,1) - activity_rate_both_sess(:,2)); % absolute activity rate diff
                abs_activity_dff_score =  abs([activity_rate_both_sess(:,1) - activity_rate_both_sess(:,2)] ./ [activity_rate_both_sess(:,1) + activity_rate_both_sess(:,2)]);  % absolute activity rate diff score
                tuning_corr = diag(corr(current_sess1(valid_cells,:)',current_sess2(valid_cells,:)','rows','pairwise')); % tuning curve correlation of corresponding cells across session halves
                
                tuning_vs_activity_rate(sess1,sess2) = corr(mean_activity_rate,tuning_corr,'rows','pairwise'); % relationship between tuning corr and activity rate
                tuning_rate_vs_activity_diff(sess1,sess2) = corr(tuning_corr,abs_activity_diff,'rows','pairwise'); % relationship between tuning corr and activity rate difference
                tuning_rate_vs_activity_diff(sess1,sess2) = corr(tuning_corr,abs_activity_diff,'rows','pairwise'); % relationship between tuning corr and activity rate difference score
                
            end
        end
        triu_id = boolean(triu(ones(6),1));
        relationships = [tuning_vs_activity_rate(triu_id),...
            tuning_rate_vs_activity_diff(triu_id),tuning_rate_vs_activity_diff(triu_id)].^2; % converte correlations into explained variance
        all_relationships(mouse,:) = mean(relationships,'omitnan'); % average explained variances values across session halves pairs
    end
    all_relationships_areas{area} = all_relationships;
end

% visualize the distribution of explained variance of the different models
pos_list = [0.1 0.575 0.25 0.325;0.4 0.575 0.25 0.325;0.7 0.575 0.25 0.325;...
    0.1 0.2 0.25 0.325;0.4 0.2 0.25 0.325;0.7 0.2 0.25 0.325];
figure('units','normalized','position',[0.25 0.2 0.45 0.55])
for area = 1:6
    current_area = all_relationships_areas{area};
    subplot(2,3,area,'units','normalized','position',pos_list(area,:))
    
    figure_boxplot(current_area);
    text(0.775,0.935,brain_areas{area},'Units','normalized','FontSize',15)
    text(0.75,0.85,['N=',num2str(size(current_area,1))],'Units','normalized','FontSize',13)
    
    set(gca,'ActivePositionProperty','position')
    if area == 4 || area == 1
        ylabel({'Explained variance by the';'tuning curve correlation (R^2)'})
    end
    if area >3
        set(gca,'xticklabels',{'Mean activity rate','Absolute activity diff.','Absolute activity diff. score'})
        
        xtickangle(30)
    else
        set(gca,'xtick',[])
    end
    ylim([-0.1 1.1])
end
suptitle('Calcium imaging - across days:')



%% Figure S1K - Between movie repeats knn decoder (Neuropixels and calcium imaging)

% KNN decoding across movie repeats using the neuropixels dataset
neuropixels_elapse_repeat_decoder_acc_areas = {};
neuropixels_elapse_repeat_decoder_acc_areas_shuffle = {};
calcium_elapse_repeat_decoder_acc_areas = {};
calcium_elapse_repeat_decoder_acc_areas_shuffle = {};
cell_cutoff = 15; % minimal number of recorded units
nat_movie = 1; % natural movie 1
num_repeats = 30; % number of movie repeats in each block
for area = 1:6 % loop over areas
    valid_mice = neuropixels_cell_count(:,area,nat_movie) >= cell_cutoff & movie_repeats(:,nat_movie) == num_repeats; % include only mice from function connectivity group with at least 15 units
    current_area = neuropixels_population_vectors(valid_mice,area,nat_movie); % subset mice that passed the requirments
    
    elapse_repeat_decoder_acc = [];
    elapse_repeat_decoder_acc_shuffle  = [];
    for mouse = 1:size(current_area,1) % loop over mice
        
        current_mouse = current_area{mouse}; % subset neuronal activity of a single mouse (#units x 30 time bins x 60 movie repeats)
        current_mouse_blockA = current_mouse(:,:,1:30); % activity during block A (repeats 1-30)
        current_mouse_blockB = current_mouse(:,:,31:60); % activity during block B (repeats 31-60)
        
        % perform knn decoding between pairs of movie repeats using data
        % from one repeat as train data and that of a proceding repeat as test
        mean_knn_decoding_accuracy = nan(30);
        mean_knn_decoding_accuracy_shuffle  = nan(30);
        for repeat1 = 1:num_repeats % loop over repeats
            for repeat2 = 1:num_repeats % loop over repeats
                if repeat1 < repeat2 % test only proceding repeats
                    clc;
                    disp(['Performing time-lapse decoding across repeats:'])
                    disp(['Dataset: Neuropixls | Stimulus: Natural movie 1 | Area: ',brain_areas{area},...
                        ' | Mouse: ',num2str(mouse),'/',num2str(length(current_area)),' | Movie repeat: ',num2str(repeat1),'/29'])
                    
                    % decoding using non-shuffled data
                    % decoding on block A
                    train_data = current_mouse_blockA(:,:,repeat1); 
                    test_data = current_mouse_blockA(:,:,repeat2);
                    mdl = fitcknn(train_data',[1:30]); % fit model to predict the time bin label based on neuronal activity in each time bin
                    prediction = [];
                    prediction(1,:) = predict(mdl,test_data'); % predict labels for each time bin in test data
                    % decoding on block B
                    train_data = current_mouse_blockB(:,:,repeat1);
                    test_data = current_mouse_blockB(:,:,repeat2);
                    mdl = fitcknn(train_data',[1:30]);% fit model to predict the time bin label based on neuronal activity in each time bin
                    prediction(2,:) = predict(mdl,test_data'); % predict labels for each time bin in test data
                    
                    accuracy = (sum(prediction == [1:30],2)./30)*100; % calculate the persentage of correctly label time bins
                    mean_knn_decoding_accuracy(repeat1,repeat2) = mean(accuracy,'omitnan'); % average performance across blocks
                    
                    % decoding using shuffled data
                    % decoding on block A
                    train_data = current_mouse_blockA(:,:,repeat1);
                    test_data = current_mouse_blockA(:,:,repeat2);
                    mdl = fitcknn(train_data',randperm(30)); % fit model to predict the time bin label based on neuronal activity in each time bin using shuffled labels
                    prediction = [];
                    prediction(1,:) = predict(mdl,test_data'); % predict labels for each time bin in test data
                    
                     % decoding on block A
                    train_data = current_mouse_blockB(:,:,repeat1);
                    test_data = current_mouse_blockB(:,:,repeat2);
                    mdl = fitcknn(train_data',randperm(30));% fit model to predict the time bin label based on neuronal activity in each time bin using shuffled labels
                    prediction(2,:) = predict(mdl,test_data');% predict labels for each time bin in test data
                    
                    accuracy = (sum(prediction == [1:30],2)./30)*100; % calculate the persentage of correctly label time bins
                    mean_knn_decoding_accuracy_shuffle(repeat1,repeat2) = mean(accuracy,'omitnan'); % average performance across blocks
                end
                
            end
        end
        
         % calculate the average performance as function of elapsed time
        for diagonal = 1:29
            elapse_repeat_decoder_acc(mouse,diagonal) = mean(diag(mean_knn_decoding_accuracy,diagonal),'omitnan'); % average performance across values of a given time interval
            elapse_repeat_decoder_acc_shuffle(mouse,diagonal) = mean(diag(mean_knn_decoding_accuracy_shuffle,diagonal),'omitnan'); % average performance across values of a given time interval
        end
        
    end
    neuropixels_elapse_repeat_decoder_acc_areas{area} = elapse_repeat_decoder_acc;
    neuropixels_elapse_repeat_decoder_acc_areas_shuffle{area} = elapse_repeat_decoder_acc_shuffle;
end

% KNN decoding across movie repeats using the calcium imaging dataset
nat_movie = 1;% natural movie 1
cell_cutoff = 20; % minimum of 20 recorded cells
for area = 1:6 % loop over areas
    valid_mice = min(calcium_excitatory_cell_count{area},[],2) >= cell_cutoff; % include only mice with atleast 20 cells in each of the recorded sessions
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie);
    
    
    elapse_repeat_decoder_acc = [];
    elapse_repeat_decoder_acc_shuffle = [];
    for mouse = 1:length(current_area) % loop over mice
        
        current_mouse = current_area{mouse}; % subset neuronal activity of a single mouse
        
        % subset only the cells that were active during session 1 (repeats 1-10)
        current_mouse_sess1 =  current_mouse(:,:,1:10);
        valid_cells = mean(squeeze(mean(current_mouse_sess1,2,'omitnan')),2,'omitnan')>0;
        valid_current_mouse_sess1 = current_mouse_sess1(valid_cells,:,:);
        
        % subset only the cells that were active during session 2 (repeats 11-20)
        current_mouse_sess2 =  current_mouse(:,:,11:20);
        valid_cells = mean(squeeze(mean(current_mouse_sess2,2,'omitnan')),2,'omitnan')>0;
        valid_current_mouse_sess2 = current_mouse_sess2(valid_cells,:,:);
        
        % subset only the cells that were active during session 3 (repeats 21-30)
        current_mouse_sess3 =  current_mouse(:,:,21:30);
        valid_cells = mean(squeeze(mean(current_mouse_sess3,2,'omitnan')),2,'omitnan')>0;
        valid_current_mouse_sess3 = current_mouse_sess3(valid_cells,:,:);
        
        % KNN decoding of time bins between pairs of movie repeats using
        % data from one movie repeat as train data to predict the time bin
        % of a proceding movie repeat
        knn_decoder_across_sessions = nan(10);
        knn_decoder_across_sessions_shuffle = nan(10);
        for repeat1 = 1:10 % loop over movie repeats
            for repeat2 = 1:10 % loop over movie repeats
                if repeat1 < repeat2 % decode proceding movie repeat
                    clc;
                    disp(['Performing time-lapse decoding across repeats:'])
                    disp(['Dataset: Neuropixls | Stimulus: Natural movie 1 | Area: ',brain_areas{area},...
                        ' | Mouse: ',num2str(mouse),'/',num2str(length(current_area)),' |  Movie repeat: ',num2str(repeat1),'/9'])
                    
                    % decoding using non-shuffled data
                    % decoding on session 1
                    train_data = valid_current_mouse_sess1(:,:,repeat1);
                    test_data = valid_current_mouse_sess1(:,:,repeat2);
                    mdl = fitcknn(train_data',[1:30]); % fit model to predict the time bin label based on neuronal activity in each time bin
                    pred = predict(mdl,test_data'); % predict labels for each time bin in test data
                    knn_decoder_across_sessions(repeat1,repeat2,1) = (sum(pred == [1:30]')./30)*100; % calculate the persentage of correctly label time bins
                    
                    % decoding on session 2
                    train_data = valid_current_mouse_sess2(:,:,repeat1);
                    test_data = valid_current_mouse_sess2(:,:,repeat2);
                    mdl = fitcknn(train_data',[1:30]); % fit model to predict the time bin label based on neuronal activity in each time bin
                    pred = predict(mdl,test_data'); % predict labels for each time bin in test data
                    knn_decoder_across_sessions(repeat1,repeat2,2) = (sum(pred == [1:30]')./30)*100; % calculate the persentage of correctly label time bins
                    
                    % decoding on session 3
                    train_data = valid_current_mouse_sess3(:,:,repeat1);
                    test_data = valid_current_mouse_sess3(:,:,repeat2);
                    mdl = fitcknn(train_data',[1:30]); % fit model to predict the time bin label based on neuronal activity in each time bin
                    pred = predict(mdl,test_data'); % predict labels for each time bin in test data
                    knn_decoder_across_sessions(repeat1,repeat2,3) = (sum(pred == [1:30]')./30)*100; % calculate the persentage of correctly label time bins
                    
                       % decoding using non-shuffled data
                    % decoding on session 1
                    train_data = valid_current_mouse_sess1(:,:,repeat1);
                    test_data = valid_current_mouse_sess1(:,:,repeat2);
                    mdl = fitcknn(train_data',randperm(30));
                    pred = predict(mdl,test_data');
                    knn_decoder_across_sessions_shuffle(repeat1,repeat2,1) = (sum(pred == [1:30]')./30)*100;
                    
                    % decoding on session 2
                    train_data = valid_current_mouse_sess2(:,:,repeat1);
                    test_data = valid_current_mouse_sess2(:,:,repeat2);
                    mdl = fitcknn(train_data',randperm(30));
                    pred = predict(mdl,test_data');
                    knn_decoder_across_sessions_shuffle(repeat1,repeat2,2) = (sum(pred == [1:30]')./30)*100;
                    
                    % decoding on session 1
                    train_data = valid_current_mouse_sess3(:,:,repeat1);
                    test_data = valid_current_mouse_sess3(:,:,repeat2);
                    mdl = fitcknn(train_data',randperm(30));
                    pred = predict(mdl,test_data');
                    knn_decoder_across_sessions_shuffle(repeat1,repeat2,3) = (sum(pred == [1:30]')./30)*100;
                end
            end
        end
        
        %average decoder performace across sessions
        mean_acc_across_mice = mean(knn_decoder_across_sessions,3,'omitnan');
        mean_acc_across_mice_shuffle = mean(knn_decoder_across_sessions_shuffle,3,'omitnan');
        
        % calculate performace as function of elapsed time
        for diagonal = 1:9 % loop over time intervals
            elapse_repeat_decoder_acc(mouse,diagonal) = nanmean(diag(mean_acc_across_mice,diagonal));
            elapse_repeat_decoder_acc_shuffle(mouse,diagonal) = nanmean(diag(mean_acc_across_mice_shuffle,diagonal));
        end
        
    end
    calcium_elapse_repeat_decoder_acc_areas{area} = elapse_repeat_decoder_acc;
    calcium_elapse_repeat_decoder_acc_areas_shuffle{area} = elapse_repeat_decoder_acc_shuffle;
end

% visualize knn decoder performance as function of elapsed time for each dataset
p =[];
z = [];
plt = [];
figure('units','normalized','position',[0.3 0.2 0.2 0.5])
for area = 1:6 % loop over area
    % neuropixels dataset - non-shuffled data
    neuropixels_current_area = neuropixels_elapse_repeat_decoder_acc_areas{area};
    neuropixels_mean_acc = mean(neuropixels_current_area,'omitnan'); % average across mice
    neuropixels_std_acc = std(neuropixels_current_area,'omitnan'); % standard deviation across mice
    neuropixels_ste_acc = neuropixels_std_acc./sqrt(size(neuropixels_current_area,1)); % standard error across mice
    
    % neuropixels dataset - shuffled data
    neuropixels_current_area_shuffle = neuropixels_elapse_repeat_decoder_acc_areas_shuffle{area};
    neuropixels_mean_acc_shuffle = mean(neuropixels_current_area_shuffle,'omitnan'); % average across mice
    neuropixels_std_acc_shuffle = std(neuropixels_current_area_shuffle,'omitnan'); % standard deviation across mice
    neuropixels_ste_acc_shuffle = neuropixels_std_acc_shuffle./sqrt(size(neuropixels_current_area_shuffle,1)); % standard error across mice
    
    % two-tailed Wilcoxon signed-rank test between minimal and maximal interval
    [p(1,area),~,stats] = signrank(neuropixels_current_area(:,1),neuropixels_current_area(:,29));
    z(1,area) = stats.zval;
    
    subplot(2,1,1)
    x = [1:length(neuropixels_mean_acc)]';
    y = neuropixels_mean_acc';
    dy = neuropixels_ste_acc';
    hold on
    fill([x;flipud(x)],[y-dy;flipud(y+dy)],colors(area,:),'linestyle','none','facealpha',0.4);
    plot(neuropixels_mean_acc,'color',colors2(area,:),'linewidth',2)
    
    x = [1:length(neuropixels_mean_acc_shuffle)]';
    y = neuropixels_mean_acc_shuffle';
    dy = neuropixels_ste_acc_shuffle';
    hold on
    fill([x;flipud(x)],[y-dy;flipud(y+dy)],[0.8 0.8 0.8],'linestyle','none','facealpha',0.4);
    plot(neuropixels_mean_acc_shuffle,'color',[0.7 0.7 0.7],'linewidth',2);
    text(0.975,0.075,'Shuffle','Units','normalized','FontSize',10)
    text(0.45,0.9,'1 repeat = 30 seconds','Units','normalized','FontSize',10)
    
    if area == 6
        plot([9 9],ylim,'--','color',[0.2 0.2 0.2],'linewidth',1.5)
    end
    
    xlim([0 30])
    ylabel({'Neuropixels';'Correct classifications (%)'})
    xlabel('Elapsed time (# of movie repeats)')
    set(gca,'xtick',[1,5,10,15,20,25,29])
    
    
    % calcium imaging dataset - non-shuffled data
    calcium_current_area = calcium_elapse_repeat_decoder_acc_areas{area};
    calcium_mean_acc = mean(calcium_current_area,'omitnan'); % average across mice
    calcium_std_acc = std(calcium_current_area,'omitnan'); % standard deviation across mice
    calcium_ste_acc = calcium_std_acc./sqrt(size(calcium_current_area,1)); % standard error across mice
    
    % calcium imaging dataset - shuffled data
    calcium_current_area_shuffle = calcium_elapse_repeat_decoder_acc_areas_shuffle{area};
    calcium_mean_acc_shuffle = nanmean(calcium_current_area_shuffle,'omitnan'); % average across mice
    calcium_std_acc_shuffle = nanstd(calcium_current_area_shuffle,'omitnan'); % standard deviation across mice
    calcium_ste_acc_shuffle = calcium_std_acc_shuffle./sqrt(size(calcium_current_area_shuffle,1)); % standard error across mice
    
    % two-tailed Wilcoxon signed-rank test between minimal and maximal interval
    [p(2,area),~,stats] = signrank(calcium_current_area(:,1),calcium_current_area(:,9));
    z(2,area) = stats.zval;
    
    subplot(2,1,2)
    x = [1:length(calcium_mean_acc)]';
    y = calcium_mean_acc';
    dy = calcium_ste_acc';
    hold on
    fill([x;flipud(x)],[y-dy;flipud(y+dy)],colors(area,:),'linestyle','none','facealpha',0.4);
    plt(area) = plot(calcium_mean_acc,'color',colors(area,:),'linewidth',2);
    plot(calcium_mean_acc,'color',colors2(area,:),'linewidth',2)
    
    x = [1:length(calcium_mean_acc_shuffle)]';
    y = calcium_mean_acc_shuffle';
    dy = calcium_ste_acc_shuffle';
    hold on
    fill([x;flipud(x)],[y-dy;flipud(y+dy)],[0.8 0.8 0.8],'linestyle','none','facealpha',0.4);
    plot(calcium_mean_acc_shuffle,'color',[0.7 0.7 0.7],'linewidth',2);
    text(0.315,0.115,'Shuffle','Units','normalized','FontSize',10)
    if area == 6
        plot([9 9],ylim,'--','color',[0.2 0.2 0.2],'linewidth',1.5)
    end
    
    xlim([0 30])
    ylabel({'Calcium imaging';'Correct classifications (%)'})
    xlabel('Elapsed time (# of movie repeats)')
    set(gca,'xtick',[1,5,10,15,20,25,29])
end
legend(plt,brain_areas(1:6),'Location','Best')
legend('boxoff')
corrected_p = bonf_holm(p);

subplot(2,1,1,'units','normalized','position',[0.2 0.6 0.6 0.35])
title('Time-lapse decoder:')
subplot(2,1,2,'units','normalized','position',[0.2 0.1 0.6 0.35])

% TODO - statistics display!!!%

%% Figure S1L - PV, ensemble rate and tuning curve correlation difference between movie repeats (Neuropixels and calcium imaging)

% analysis of neuropixels dataset
cell_cutoff = 15; % minimum number of recorded units from each area
nat_movie = 1; % natural movie 1
num_repeats = 30; % number of movie repeats in each block
neuropixels_elapse_repeat_areas = {};
for area = 1:6 % loop over areas
    valid_mice = neuropixels_cell_count(:,area,nat_movie) >= cell_cutoff & movie_repeats(:,nat_movie) == num_repeats; % include only mice from the functional connectivity group with at least 15 units in the selected area
    current_area = neuropixels_population_vectors(valid_mice,area,nat_movie); % subset only mice that passed the requirments of 'valid_mice'
    
    elapse_repeat_pv = [];
    elapse_repeat_rate  = [];
    elapse_repeat_tuning  = [];
    for mouse = 1:size(current_area,1) % loop over mice
        clc;
        disp(['Calculating PV, ensemble rate and tuning curve correlations across movie repeats:'])
        disp(['Dataset: Neuropixls | Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse}; % subset neuronal activity of a single mouse (#units x 30 time bins x 60 movie repeats)
        current_mouse_blockA = current_mouse(:,:,1:30); % neuronal activity during block A (repeats 1-30)
        current_mouse_blockB = current_mouse(:,:,31:60); % neuronal activity during block B (repeats 31-60)
        
        % calculate pv, rate and tuning corr across pairs of movie repeats
        rate_corr = [];
        tuning_corr = [];
        pv_corr = [];
        for repeat1 = 1:num_repeats % loop over movie repeats
            for repeat2 = 1:num_repeats % loop over movie repeats
                
                pv_corr(repeat1,repeat2,1) = mean(diag(corr(current_mouse_blockA(:,:,repeat1),current_mouse_blockA(:,:,repeat2))),'omitnan'); % mean pv corr across corresponding time bins between pair of movie repeats in block A
                pv_corr(repeat1,repeat2,2) = mean(diag(corr(current_mouse_blockB(:,:,repeat1),current_mouse_blockB(:,:,repeat2))),'omitnan'); % mean pv corr across corresponding time bins between pair of movie repeats in block B
                
                rate_corr(repeat1,repeat2,1) = corr(mean(current_mouse_blockA(:,:,repeat1),2,'omitnan'),mean(current_mouse_blockA(:,:,repeat2),2,'omitnan')); % average the activity across time bins and calculate the ensemble rate corr between pair of movie repeats in block A
                rate_corr(repeat1,repeat2,2) = corr(mean(current_mouse_blockB(:,:,repeat1),2,'omitnan'),mean(current_mouse_blockB(:,:,repeat2),2,'omitnan')); % average the activity across time bins and calculate the ensemble rate corr between pair of movie repeats in block B
                
                tuning_corr(repeat1,repeat2,1) = median(diag(corr(current_mouse_blockA(:,:,repeat1)',current_mouse_blockA(:,:,repeat2)')),'omitnan'); % median tuning curve corr across corresponding units between pair of movie repeats in block A
                tuning_corr(repeat1,repeat2,2) = median(diag(corr(current_mouse_blockB(:,:,repeat1)',current_mouse_blockB(:,:,repeat2)')),'omitnan'); % median tuning curve corr across corresponding units between pair of movie repeats in block B
                
            end
        end
        
        % average correlations across blocks
        mean_pv_corr = mean(pv_corr,3,'omitnan');
        mean_rate_corr = mean(rate_corr,3,'omitnan');
        mean_tuning_corr = mean(tuning_corr,3,'omitnan');
        
        % calculate pv,rate and tuning corr as function of elapsed time
        for diagonal = 1:29
            elapse_repeat_pv(mouse,diagonal) = mean(diag(mean_pv_corr,diagonal),'omitnan'); % average across pv corr values for a given time interval
            elapse_repeat_rate(mouse,diagonal) = mean(diag(mean_rate_corr,diagonal),'omitnan'); % average across rate corr values for a given time interval
            elapse_repeat_tuning(mouse,diagonal) = mean(diag(mean_tuning_corr,diagonal),'omitnan'); % average across tuning corr values for a given time interval
        end
        
    end
    neuropixels_elapse_repeat_areas(area,:) = {elapse_repeat_pv,elapse_repeat_rate,elapse_repeat_tuning};
end

% analysis of calcium imaging dataset
nat_movie = 1; % natural movie 1
cell_cutoff = 20; % minimum number of cells per imaging session
calcium_elapse_repeat_areas = {};
for area = 1:6 % loop over areas
    valid_mice = min(calcium_excitatory_cell_count{area},[],2) >= cell_cutoff; % subset only mice with atleast 20 cells in each session
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie); % subset mice that passed the requirements of 'valid_mice'
    
   
    elapse_repeat_pv = [];
    elapse_repeat_rate = [];
    elapse_repeat_tuning = [];
    for mouse = 1:length(current_area) % loop over mice
        clc;
        disp(['Calculating PV, ensemble rate and tuning curve correlations across movie repeats:'])
        disp(['Dataset: Calcium imaging | Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse}; % subset neuronal activity of a single mouse (#cells x 30 time bins x 30 movie repeats
        
        current_mouse_sess1 =  current_mouse(:,:,1:10); % neuronal activity of session 1 (repeats 1-10)
        valid_cells = mean(squeeze(mean(current_mouse_sess1,2,'omitnan')),2,'omitnan')>0; % find cells that were active in session 1
        valid_current_mouse_sess1 = current_mouse_sess1(valid_cells,:,:); % subset cells that were active in session 1
        
        current_mouse_sess2 =  current_mouse(:,:,11:20); % neuronal activity of session 2 (repeats 11-20)
        valid_cells = mean(squeeze(mean(current_mouse_sess2,2,'omitnan')),2,'omitnan')>0; % find cells that were active in session 2
        valid_current_mouse_sess2 = current_mouse_sess2(valid_cells,:,:); % subset cells that were active in session 2
        
        current_mouse_sess3 =  current_mouse(:,:,21:30); % neuronal activity of session 3 (repeats 21-30)
        valid_cells = mean(squeeze(mean(current_mouse_sess3,2,'omitnan')),2,'omitnan')>0; % find cells that were active in session 3
        valid_current_mouse_sess3 = current_mouse_sess3(valid_cells,:,:); % subset cells that were active in session 3
        
        % calculate the pv, rate and tuning correlation between movie repeats
        % within the same sessions  
        pv_corr = [];
        rate_corr = [];
        tuning_corr = [];
        for repeat1 = 1:10 % loop over movie repeats
            for repeat2 = 1:10 % loop over movie repeats
                pv_corr(repeat1,repeat2,1) = mean(diag(corr(current_mouse_sess1(:,:,repeat1),current_mouse_sess1(:,:,repeat2))),'omitnan'); % calculate the pv corr between time bins of different movie repeats in session 1 and average across corresponding time bins
                pv_corr(repeat1,repeat2,2) = mean(diag(corr(current_mouse_sess2(:,:,repeat1),current_mouse_sess2(:,:,repeat2))),'omitnan'); % calculate the pv corr between time bins of different movie repeats in session 2 and average across corresponding time bins
                pv_corr(repeat1,repeat2,3) = mean(diag(corr(current_mouse_sess3(:,:,repeat1),current_mouse_sess3(:,:,repeat2))),'omitnan'); % calculate the pv corr between time bins of different movie repeats in session 3 and average across corresponding time bins
                
                rate_corr(repeat1,repeat2,1) = corr(mean(current_mouse_sess1(:,:,repeat1),2,'omitnan'),mean(current_mouse_sess1(:,:,repeat2),2,'omitnan')); % average the neuronal activity across time bins and calculate the rate corr between movie repeats in session 1
                rate_corr(repeat1,repeat2,2) = corr(mean(current_mouse_sess2(:,:,repeat1),2,'omitnan'),mean(current_mouse_sess2(:,:,repeat2),2,'omitnan')); % average the neuronal activity across time bins and calculate the rate corr between movie repeats in session 2
                rate_corr(repeat1,repeat2,3) = corr(mean(current_mouse_sess3(:,:,repeat1),2,'omitnan'),mean(current_mouse_sess3(:,:,repeat2),2,'omitnan')); % average the neuronal activity across time bins and calculate the rate corr between movie repeats in session 3
                
                tuning_corr(repeat1,repeat2,1) = mean(diag(corr(current_mouse_sess1(:,:,repeat1)',current_mouse_sess1(:,:,repeat2)')),'omitnan'); % calculate the tuning curve between pair of movie repeats in session 1 and calculate the median across corresponding cells
                tuning_corr(repeat1,repeat2,2) = mean(diag(corr(current_mouse_sess2(:,:,repeat1)',current_mouse_sess2(:,:,repeat2)')),'omitnan'); % calculate the tuning curve between pair of movie repeats in session 2 and calculate the median across corresponding cells
                tuning_corr(repeat1,repeat2,3) = mean(diag(corr(current_mouse_sess3(:,:,repeat1)',current_mouse_sess3(:,:,repeat2)')),'omitnan'); % calculate the tuning curve between pair of movie repeats in session 3 and calculate the median across corresponding cells
                
            end
        end
        % average correlation values across sessions
        pv_corr = mean(pv_corr,3,'omitnan'); 
        rate_corr = mean(rate_corr,3,'omitnan');
        tuning_corr = mean(tuning_corr,3,'omitnan');
        
        % calculate the pv, rate and tuning corr as function of elapsed time
        for diagonal = 1:9 % loop over time intervals
            elapse_repeat_pv(mouse,diagonal) = mean(diag(pv_corr,diagonal),'omitnan'); % average pv corr values of a single time interval
            elapse_repeat_rate(mouse,diagonal) = mean(diag(rate_corr,diagonal),'omitnan'); % average rate corr values of a single time interval
            elapse_repeat_tuning(mouse,diagonal) = mean(diag(tuning_corr,diagonal),'omitnan'); % average tuning corr values of a single time interval
        end
        
    end
    calcium_elapse_repeat_areas(area,:) = {elapse_repeat_pv,elapse_repeat_rate,elapse_repeat_tuning};
end

% visualize the pv, rate and tuning corr as function of elapsed time for
% all visual areas and datasets
ylims = {[-0.095 0],[-0.1 0],[-0.12 0];...
    [-0.06 0],[-0.125 0],[-0.025 0.005]};
figure('units','normalized','position',[0.2 0.2 0.475 0.45])
for measurment = 1:3 % loop over stability measurments
    for area = 1:6 % loop over areas
        % visualize neuropixels data
        neuropixels_current_area = neuropixels_elapse_repeat_areas{area,measurment}; % subset values of a single area
        neuropixels_current_area_diff = neuropixels_current_area-neuropixels_current_area(:,1); % scale the values based on the values of the minimal interval
        neuropixels_mean_stability = mean(neuropixels_current_area_diff,'omitnan'); % calculate the mean across mice
        neuropixels_std_stability = std(neuropixels_current_area_diff,'omitnan'); % calculate the standard deviation across mice
        neuropixels_ste_stability = neuropixels_std_stability./sqrt(size(neuropixels_current_area,1)); % calculate the standard error across mice
        
        subplot(2,3,measurment)
        x = [1:length(neuropixels_mean_stability)]';
        y = neuropixels_mean_stability';
        dy = neuropixels_ste_stability';
        hold on
        fill([x;flipud(x)],[y-dy;flipud(y+dy)],colors(area,:),'linestyle','none','facealpha',0.4);
        if measurment == 1
            plt(area) = plot(neuropixels_mean_stability,'color',colors(area,:),'linewidth',2);
        end
        plot(neuropixels_mean_stability,'color',colors2(area,:),'linewidth',2)
        
        ylim(ylims{1,measurment})
        xlim([0 30])
        if measurment == 1
            ylabel({'Neuropixels';'Correlation difference'})
            title('Population vectors')
        elseif measurment == 2
            title('Ensemble rate')
        elseif measurment == 3
            title('Tuning curve')
        end
        set(gca,'xtick',[1,5,10,15,20,25,29])
        
        
        
        % visualize calcium imaging data
        calcium_current_area = calcium_elapse_repeat_areas{area,measurment}; % subset values of a single area
        calcium_current_area_diff = calcium_current_area-calcium_current_area(:,1); % scale the values based on the values of the minimal interval
        calcium_mean_stability = mean(calcium_current_area_diff,'omitnan'); % calculate the mean across mice
        calcium_std_stability = std(calcium_current_area_diff,'omitnan'); % calculate the standard deviation across mice
        calcium_ste_stability = calcium_std_stability./sqrt(size(calcium_current_area_diff,1)); % calculate the standard error across mice
        
        subplot(2,3,measurment+3)
        x = [1:length(calcium_mean_stability)]';
        y = calcium_mean_stability';
        dy = calcium_ste_stability';
        hold on
        fill([x;flipud(x)],[y-dy;flipud(y+dy)],colors(area,:),'linestyle','none','facealpha',0.4);
        
        plot(calcium_mean_stability,'color',colors2(area,:),'linewidth',2)
        xlim([0.5 9.5])
        ylim(ylims{2,measurment})
        
        xlabel('Elapsed time (# of movie repeats)')
        set(gca,'xtick',[1:9])
    end
    if measurment == 1
        text(0.05,0.1,'1 repeat = 30 seconds','Units','normalized','FontSize',10)
        ylabel({'Calcium imaging';'Correlation difference'})
    end
end
subplot(2,3,1)
legend(plt,brain_areas(1:6),'Location','Best')
legend('boxoff')


%% Figure S2A - Running speed across movie repeats within a block
repeat_num = 30; % number of movie repeats
nat_movie = 1; % natural movie 1
valid_mice = neuropixels_running_speed(movie_repeats(:,nat_movie) == repeat_num,nat_movie); % subset mice from the 'Functional connectivity' group

running_speed_across_mice = [];
for mouse = 1:size(valid_mice,1) % loop over mice
    running_speed_across_mice(:,:,mouse) = valid_mice{mouse,nat_movie}; % subset running speed of a single mouse
end

mean_running_speed_across_blocks = squeeze(mean(running_speed_across_mice,2,'omitnan'))'; % calculate mean running speed across blocks 
mean_running_speed = mean(mean_running_speed_across_blocks,'omitnan'); % calculate mean running speed for each movie repeat across mice
std_running_speed = std(mean_running_speed_across_blocks,'omitnan'); % calculate standard deviation across mice
ste_running_speed = std_running_speed./sqrt(size(mean_running_speed_across_blocks,1)); % calculate standard error across mice

% visualize the mean running speed in each movie repeat within blocks
figure('units','normalized','position',[0.3 0.3 0.25 0.35])
hold on
y = [ones(9,1).*0;ones(9,1).*25];
x = [0:8,8:-1:0];
fill(x',y',[0.8 0.8 0.8],'linestyle','none','facealpha',0.1);
errorbar(mean_running_speed,ste_running_speed,'color',[0.8 0.8 0.8],'capsize',0','linewidth',3)
plot(mean_running_speed,'color',[0.4 0.4 0.4],'linewidth',3)
plot([8 8],ylim,'--','color',[0.2 0.2 0.2],'linewidth',2)
set(gca,'xtick',[1,5:5:30])
ylabel('Mean running speed (cm/sec)')
xlabel('Movie repeat')
ylim([0 25])
xlim([0 31])

% perform statistical tests for the difference between running speeds
% across different movie repeats
sig_mat = [];
for trial1 = 1:30 % loop over movie repeats
    for trial2 = 1:30 % loop over movie repeats
        sig_mat(trial1,trial2) = ttest(mean_running_speed_across_blocks(:,trial1),mean_running_speed_across_blocks(:,trial2)); % perform two-sided paired t-test between movie repeats
    end
end

figure('units','normalized','position',[0.575 0.425 0.1 0.175])
imagesc(sig_mat)
hold on
plot([8 8],ylim,'--','color',[0.2 0.2 0.2],'linewidth',2)
plot(xlim,[8 8],'--','color',[0.2 0.2 0.2],'linewidth',2)
ylabel('Movie repeat')
xlabel('Movie repeat')
colormap(sig_colormap)


%% Figure S2B - Pupil area across movie repeats within a block
repeat_num = 30; % number of movie repeats
nat_movie = 1; % natural movie 1
valid_mice = neuropixels_pupil_size(movie_repeats(:,nat_movie) == repeat_num,nat_movie); % subset mice from the 'Functional connectivity' group

pupil_area_across_mice = [];
sub = 1;
for mouse = 1:size(valid_mice,1) % loop over mice
    if ~isempty(valid_mice{mouse,nat_movie}) % check if eye tracking data is available for current mouse
        pupil_area_across_mice(:,:,sub) = valid_mice{mouse,nat_movie}; % subset pupil size of single mouse
        sub = sub + 1;
    end
end

mean_pupil_area_across_blocks = squeeze(mean(pupil_area_across_mice,2,'omitnan'))'; % calculate mean pupil size across blocks 
mean_pupil_area = mean(mean_pupil_area_across_blocks,'omitnan'); % calculate mean pupil size for each movie repeat across mice
std_pupil_area = std(mean_pupil_area_across_blocks,'omitnan'); % calculate standard deviation across mice
ste_pupil_area = std_pupil_area./sqrt(size(mean_pupil_area_across_blocks,1)); % calculate standard error across mice

% visualize the mean pupil area in each movie repeat within blocks
figure('units','normalized','position',[0.3 0.3 0.25 0.35])
hold on
y = [ones(9,1).*250;ones(9,1).*650];
x = [0:8,8:-1:0];
fill(x',y',[0.8 0.8 0.8],'linestyle','none','facealpha',0.1);
errorbar(mean_pupil_area,ste_pupil_area,'color',[0.8 0.8 0.8],'capsize',0','linewidth',3)
plot(mean_pupil_area,'color',[0.4 0.4 0.4],'linewidth',3)
plot([8 8],ylim,'--','color',[0.2 0.2 0.2],'linewidth',2)
set(gca,'xtick',[1,5:5:30])
ylabel('Mean pupil area (a.u.)')
xlabel('Movie repeat')
ylim([250 650])
xlim([0 31])

% perform statistical tests for the difference between pupil area
% across different movie repeats
sig_mat = [];
for trial1 = 1:30 % loop over movie repeats
    for trial2 = 1:30 % loop over movie repeats
        sig_mat(trial1,trial2) = ttest(mean_pupil_area_across_blocks(:,trial1),mean_pupil_area_across_blocks(:,trial2)); % perform two-sided paired t-test between movie repeats
    end
end

figure('units','normalized','position',[0.575 0.425 0.1 0.175])
imagesc(sig_mat)
hold on
plot([8 8],ylim,'--','color',[0.2 0.2 0.2],'linewidth',2)
plot(xlim,[8 8],'--','color',[0.2 0.2 0.2],'linewidth',2)
ylabel('Movie repeat')
xlabel('Movie repeat')
colormap(sig_colormap)

%% Figure S2C - mean activity rate across movie repeats within a block
cell_cutoff = 15; % minimum number of units recorded from each area
nat_movie = 1; % natural movie 1
num_repeats = 30; % number of movie repeats in each block
mean_activity_all_areas = {};
for area = 1:6 % loop over areas
    valid_mice = neuropixels_cell_count(:,area,nat_movie) >= cell_cutoff & movie_repeats(:,nat_movie) == num_repeats; % include only mice from the 'Functional connectivity' grou
    current_area = neuropixels_population_vectors(valid_mice,area,nat_movie); % subset mive that passed the requirments of 'valid_mice'
    
    % calculate the mean activity rate across units in each movie repeat
    mean_activity_all_mice = [];
    for mouse = 1:size(current_area,1) % loop over mice
        current_mouse = current_area{mouse}*30; % subset neuronal activity of a single mouse
        
        mean_activity_blockA = squeeze(mean(mean(current_mouse(:,:,1:30),2,'omitnan'),'omitnan')); % calculate the average activity across time bins and units for each movie repeat in block A (repeats 1-30)
        mean_activity_blockB = squeeze(mean(mean(current_mouse(:,:,31:60),2,'omitnan'),'omitnan')); % calculate the average activity across time bins and units for each movie repeat in block A (repeats 1-30)
        mean_activity_all_mice(mouse,:) = mean([mean_activity_blockA,mean_activity_blockB]','omitnan'); % average activity rates across blocks
    end
    
    mean_activity_all_areas{area} = mean_activity_all_mice;
end

% visualize the mean activity rates for each movie repeat across units
figure('units','normalized','position',[0.3 0.3 0.35 0.375])
for area = 1:6 % loop over areas
    current_area =  mean_activity_all_areas{area}; % subset neuronal activity of a single area
    mean_activity = mean(current_area,'omitnan'); % mean across mice
    std_activity = std(current_area,'omitnan'); % standard deviation across mice
    ste_activity = std_activity./sqrt(size(current_area,1)); % standard error across mice
    
    subplot(2,3,area)
    hold on
    y = [ones(9,1).*5.35;ones(9,1).*9];
    x = [0:8,8:-1:0];
    fill(x',y',[0.8 0.8 0.8],'linestyle','none','facealpha',0.1);
    errorbar(mean_activity,ste_activity,'color',[0.8 0.8 0.8],'capsize',0','linewidth',3)
    plot(mean_activity,'color',[0.4 0.4 0.4],'linewidth',3)
    plot([8 8],ylim,'--','color',[0.2 0.2 0.2],'linewidth',2)
    text(0.65,0.875,brain_areas{area},'Units','normalized','FontSize',15)
    text(0.6,0.75,['N=',num2str(size(current_area,1))],'Units','normalized','FontSize',13)
    
    if area >3
        set(gca,'xtick',[1,5:5:30])
    else
        set(gca,'xtick',[])
    end
    if area == 4
        ylabel('Mean activity rate (spike/sec)')
    elseif area == 5
        xlabel('Movie repeat')
    end
    ylim([5.35 9])
    xlim([0 31])
end

%% Figure S2D - Ensemble rate correlation on subset of movie repeats (repeats 9-30)
nat_movie = 1; % natural movie 1
cell_cutoff = 15; % cell count threshold of 15 units
num_repeats = 30; % number of movie repeats in each block
repeat_lim = {1:30,9:30}; % range of repeats for each subset
subset_list = {'Full','Repeats 9-30'};
elapse_repeat_rate_corr_area = {}; ;% define an empty variable for ensemble rate corr across areas

% calculate for each visual area of each mouse the ensemble rate correlation between
% pairs of movie repeats and as a function of elapsed time
for area = 1:6 % loop over areas
    % valid_mice - only mice from the 'Functional connectivity' group with at least 15 units
    valid_mice = neuropixels_cell_count(:,area,nat_movie) >= cell_cutoff & movie_repeats(:,nat_movie) == num_repeats;
    current_area = neuropixels_population_vectors(valid_mice,area,nat_movie);  % subset all mice that met the requirements of 'valid mice'
    for repeat_span = 1:2 % perform analysis on different subsets of neuronal data
        current_span = repeat_lim{repeat_span}; % range of repeats for current subset
        
        elapse_repeat_rate_corr = []; % define an empty variable for ensemble rate corr matrices across mice
        for mouse = 1:size(current_area,1) % loop over mice
            
            clc;
            disp(['Calculating ensemble rate correlation between movie repeats:'])
            disp(['Dataset: Neuropixels | Stimulus: Natural movie 1 | Area: ',brain_areas{area},...
                ' | Subset: ',subset_list{repeat_span},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
            
            current_mouse = current_area{mouse}; % subset the neuronal activity of a single mouse (size of #units by 30 time bins by 60 movie repeats)
            
            current_mouse_blockA = current_mouse(:,:,1:30); % subset neuronal activity in each movie repeat during block A (repeats 1-30)
            current_mouse_blockA = squeeze(mean(current_mouse_blockA(:,:,current_span),2,'omitnan')); % calculate aaverage activity across time bins of selected movie repeats
             
            current_mouse_blockB = current_mouse(:,:,31:60); % subset neuronal activity in each movie repeat during block A (repeats 31-60)
            current_mouse_blockB = squeeze(mean(current_mouse_blockB(:,:,current_span),2,'omitnan')); % calculate aaverage activity across time bins of selected movie repeats
            
             % calculate ensemble rate correlation between pairs of movie
            % repeats in each block
            mean_rate_corr = [];
            mean_rate_corr(:,:,1) = corr(current_mouse_blockA); % ensemble rate correlation for block A
            mean_rate_corr(:,:,2) = corr(current_mouse_blockB); % ensemble rate correlation for block B
            mean_rate_corr = mean(mean_rate_corr,3,'omitnan'); % average ensemble rate correlation across blocks
            
            % calculate ensemble rate correlation as function of elapsed time
            for diagonal = 1:length(current_span)-1 % loop over diagonals
                elapse_repeat_rate_corr(mouse,diagonal) = mean(diag(mean_rate_corr,diagonal),'omitnan'); % average pv corr values across diagonals
            end
            
        end
        elapse_repeat_rate_corr_area{repeat_span,area} = elapse_repeat_rate_corr; % store ensemble rate corr values across mice for each area
    end
end

% friedman test parameters
pvalues = []; % define empty variable for the pvalues
df = []; % define empty variable for the degrees of freedom
chi = []; % define empty variable for the chi square values
plt = [];
ylims = [0.88 0.98;0.87 1;0.86 0.98;0.88 0.98;0.87 0.98;0.87 0.98];
figure('units','normalized','position',[0.3 0.3 0.35 0.385]) % visualization of ensemble rate corr as function of time across areas and mice
for area = 1:6 % loop over areas
    for repeat_span = 1:2 % loop over repeats ranges
        current_area = elapse_repeat_rate_corr_area{repeat_span,area}; % subset values of specific visual area
        
        [pvalues(repeat_span,area),tbl,stats] = friedman(current_area,[1],'off'); % perform friedman test for main effect of elapsed time
        df(repeat_span,area) = tbl{2,3}; % save degrees of freedom
        chi(repeat_span,area) = tbl{2,5}; % save chi square values
        
        subplot(2,3,area)
        hold on
        mean_elapse_repeat = mean(current_area,'omitnan'); % calculate mean ensemble rate corr across mice
        std_elapse_repeat = std(current_area,'omitnan'); % calculate standard deviation across mice
        ste_elapse_repeat = std_elapse_repeat./sqrt(size(current_area,1)); % calculate standard error across mice
        
        x = [1:length(mean_elapse_repeat)]';
        y = mean_elapse_repeat';
        dy = ste_elapse_repeat';
        if repeat_span == 1
            fill([x;flipud(x)],[y-dy;flipud(y+dy)],[0.6 0.6 0.6],'linestyle','none','facealpha',0.4);
            plt(repeat_span) = plot(y,'color',[0.4 0.4 0.4],'linewidth',3);
        elseif repeat_span == 2
            fill([x;flipud(x)],[y-dy;flipud(y+dy)],colors(area,:),'linestyle','none','facealpha',0.4);
            plt(repeat_span) = plot(y,'color',colors2(area,:),'linewidth',3);
        end
        text(0.7,0.875,brain_areas{area},'Units','normalized','FontSize',15)
        
        if area == 4 || area == 1
            ylabel('Ensemble rate correlation')
        elseif area == 5
            xlabel({'Elapsed time (# of movie repeats)';'1 repeat = 30 seconds'})
        end
        
        if area >3
            set(gca,'xtick',[1 10 20 29])
        else
            set(gca,'xtick',[])
        end
        ylim(ylims(area,:))
        xlim([0 30])
        
    end
    legend(plt,{'All reps','Reps 9-30'},'Location','Best');
    legend('boxoff')
end

% correct for multiple comparisons using bonferroni-holm method
corrected_pval = bonf_holm(pvalues(2,:));

% define statistics summary table
VarNames = {'area','df','chi','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas(1:6)'),df(2,:)',chi(2,:)',pvalues(2,:)',corrected_pval(:),'VariableNames',VarNames);

% display statistics summary table
clc;
disp(['Ensemble rate correlation as a function of elapsed time for subsampled data (repeats 9-30)'])
disp(['Friedman�s tests with Holm�Bonferroni correction:'])
disp(statistics)

%% Figure S2E - Tuning curve correlation on subset of movie repeats (repeats 9-30)
nat_movie = 1; % natural movie 1
cell_cutoff = 15; % cell count threshold of 15 units
num_repeats = 30; % number of movie repeats in each block
repeat_lim = {1:30,9:30};
subset_list = {'Full','Repeats 9-30'};
elapse_repeat_tuning_corr_area = {}; % define an empty variable for tuning curve corr across areas

% calculate for each visual area of each mouse the tuning curve correlation between
% pairs of movie repeats as function of elapsed time
for area = 1:6 % loop over brain areas
    % valid_mice - only mice from the 'Functional connectivity' group with at least 15 units
    valid_mice = neuropixels_cell_count(:,area,nat_movie) >= cell_cutoff & movie_repeats(:,nat_movie) == num_repeats;
    current_area = neuropixels_population_vectors(valid_mice,area,nat_movie); % subset all mice that met the requirements of 'valid mice'
    for repeat_span = 1:2
        current_span = repeat_lim{repeat_span};
        
        elapse_repeat_tuning_corr = []; % define an empty variable for tuning curve corr matrices across mice 
        for mouse = 1:size(current_area,1)  % loop over mice
            
            clc;
            disp(['Calculating tuning curve correlation between movie repeats:'])
            disp(['Dataset: Neuropixels | Stimulus: Natural movie 1 | Area: ',brain_areas{area},...
                ' | Subset: ',subset_list{repeat_span},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
            
            current_mouse = current_area{mouse};  % subset a single mouse
            
            current_mouse_blockA = current_mouse(:,:,1:30); % subset neuronal activity during the first block (repeats 1-30)
            current_mouse_blockB = current_mouse(:,:,31:60); % subset neuronal activity during the second block (repeats 31-60)
            
            tuning_corr = [];  % define an empty variable for average tuning curve across movie repeats
            for repeat1 = 1:length(current_span) % loop over movie repeats
                for repeat2 = 1:length(current_span) % loop over movie repeats
                    tuning_corr(repeat1,repeat2,1) = median(diag(corr(current_mouse_blockA(:,:,repeat1)',current_mouse_blockA(:,:,repeat2)')),'omitnan'); % calculate tuning curve corr batween movie repeats of block A and calculate the median across corresponding units
                    tuning_corr(repeat1,repeat2,2) = median(diag(corr(current_mouse_blockB(:,:,repeat1)',current_mouse_blockB(:,:,repeat2)')),'omitnan'); % calculate tuning curve corr batween movie repeats of block B and calculate the median across corresponding units
                end
            end
            mean_tuning_corr = mean(tuning_corr,3,'omitnan'); % average across blocks
            
            % calculate tuning curve correlation as function of elapsed time
            for diagonal = 1:length(current_span)-1 % loop over diagonals
                elapse_repeat_tuning_corr(mouse,diagonal) = mean(diag(mean_tuning_corr,diagonal),'omitnan'); % average tuning curve corr values across diagonals
            end
            
        end
        elapse_repeat_tuning_corr_area{repeat_span,area} = elapse_repeat_tuning_corr; % store mean tuning curve values across mice for each area
    end
end


% friedman test parameters
pvalues = []; % define empty variable for the pvalues
df = []; % define empty variable for the degrees of freedom
chi = []; % define empty variable for the chi square values
ylims = [0.3 0.45;0.325 0.475;0.225 0.425;0.25 0.4;0.15 0.3;0.275 0.4];
figure('units','normalized','position',[0.3 0.3 0.35 0.385])
for area = 1:6
    plt = [];
    for repeat_span = 1:2
        current_area = elapse_repeat_tuning_corr_area{repeat_span,area};
        
        [pvalues(repeat_span,area),tbl,stats] = friedman(current_area,[1],'off');
        df(repeat_span,area) = tbl{2,3};
        chi(repeat_span,area) = tbl{2,5};
        
        subplot(2,3,area)
        hold on
        mean_elapse_repeat = nanmean(current_area);
        std_elapse_repeat = nanstd(current_area);
        ste_elapse_repeat = std_elapse_repeat./sqrt(size(current_area,1));
        
        x = [1:length(mean_elapse_repeat)]';
        y = mean_elapse_repeat';
        dy = ste_elapse_repeat';
        if repeat_span == 1
            fill([x;flipud(x)],[y-dy;flipud(y+dy)],[0.6 0.6 0.6],'linestyle','none','facealpha',0.4);
            plt(repeat_span) = plot(y,'color',[0.4 0.4 0.4],'linewidth',3);
        elseif repeat_span == 2
            fill([x;flipud(x)],[y-dy;flipud(y+dy)],colors(area,:),'linestyle','none','facealpha',0.4);
            plt(repeat_span) = plot(y,'color',colors2(area,:),'linewidth',3);
        end
        text(0.7,0.875,brain_areas{area},'Units','normalized','FontSize',15)
        
        if area == 4 || area == 1
            ylabel('Tuning curve correlation')
        elseif area == 5
            xlabel({'Elapsed time (# of movie repeats)';'1 repeat = 30 seconds'})
        end
        
        if area >3
            set(gca,'xtick',[1 10 20 29])
        else
            set(gca,'xtick',[])
        end
        ylim(ylims(area,:))
        xlim([0 30])
        
    end
    legend(plt,{'All reps','Reps 9-30'},'Location','Best');
    legend('boxoff')
end

corrected_pval = bonf_holm(pvalues(2,:));

VarNames = {'area','df','chi','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas(1:6)'),df(2,:)',chi(2,:)',pvalues(2,:)',corrected_pval(:),'VariableNames',VarNames);

clc;
disp(['Tuning curve correlation as a function of elapsed time for subsampled data (repeats 9-30)'])
disp(['Friedman�s tests with Holm�Bonferroni correction:'])
disp(statistics)


%% Figure S2F - Ensemble rate correlation on subset of non-adapted units
cell_cutoff = 15;
nat_movie = 1;
num_repeats = 30;
elapsed_rate_corr_area = {};
for area = 1:6
    valid_mice = neuropixels_cell_count(:,area,nat_movie) >= cell_cutoff & movie_repeats(:,nat_movie) == num_repeats;
    current_area = neuropixels_population_vectors(valid_mice,area,nat_movie);
    
    rate_corr_subset = [];
    rate_corr = [];
    included_units = [];
    for mouse = 1:size(current_area,1)
        clc;
        disp(['Calculating ensemble rate correlation between movie repeats:'])
        disp(['Dataset: Neuropixels | Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse};
        
        pval_across_units = [];
        for unit = 1:size(current_mouse,1)
            current_unit = squeeze(current_mouse(unit,:,:));
            
            current_unit_blockA = nanmean(current_unit(:,1:30));
            current_unit_blockB = nanmean(current_unit(:,31:60));
            
            [~,pval_across_units(unit,1)] = ttest2(current_unit_blockA(:,1:15),current_unit_blockA(:,16:30),'tail','right');
            [~,pval_across_units(unit,2)] = ttest2(current_unit_blockB(:,1:15),current_unit_blockB(:,16:30),'tail','right');
            
        end
        
        valid_units = pval_across_units > 0.05;
        
        if sum(sum(valid_units) == 0) == 0
            included_units(mouse,:)= sum(valid_units)./size(pval_across_units,1);
            
            current_mouse_blockA = squeeze(nanmean(current_mouse(:,:,1:30),2));
            current_mouse_blockB = squeeze(nanmean(current_mouse(:,:,31:60),2));
            
            valid_current_mouse_blockA = current_mouse_blockA(valid_units(:,1),:);
            valid_current_mouse_blockB = current_mouse_blockB(valid_units(:,2),:);
            
            rate_corr = [];
            rate_corr(:,:,1) = corr(current_mouse_blockA);
            rate_corr(:,:,2) = corr(current_mouse_blockB);
            rate_corr = nanmean(rate_corr,3);
            
            rate_corr_subsample = [];
            rate_corr_subsample(:,:,1) = corr(valid_current_mouse_blockA);
            rate_corr_subsample(:,:,2) = corr(valid_current_mouse_blockB);
            rate_corr_subsample = nanmean(rate_corr_subsample,3);
            
            for diagonal = 1:29
                elapsed_rate_corr(mouse,diagonal)  = nanmean(diag(rate_corr,diagonal));
                elapsed_rate_corr_subset(mouse,diagonal) = nanmean(diag(rate_corr_subsample,diagonal));
            end
        end
    end
    elapsed_rate_corr_area(area,:) = {elapsed_rate_corr,elapsed_rate_corr_subset};
end


pvalues = [];
df = [];
chi = [];
ylims = [0.88 0.98;0.88 0.99;0.86 0.98;0.87 0.98;0.87 0.98;0.87 0.98];
figure('units','normalized','position',[0.3 0.3 0.35 0.385])
for area = 1:6
    for subset = 1:2
        current_area = elapsed_rate_corr_area{area,subset};
        
        [pvalues(subset,area),tbl,stats] = friedman(current_area,[1],'off');
        df(subset,area) = tbl{2,3};
        chi(subset,area) = tbl{2,5};
        
        subplot(2,3,area)
        hold on
        mean_elapse_repeat = nanmean(current_area);
        std_elapse_repeat = nanstd(current_area);
        ste_elapse_repeat = std_elapse_repeat./sqrt(size(current_area,1));
        
        x = [1:length(mean_elapse_repeat)]';
        y = mean_elapse_repeat';
        dy = ste_elapse_repeat';
        if subset == 1
            fill([x;flipud(x)],[y-dy;flipud(y+dy)],[0.6 0.6 0.6],'linestyle','none','facealpha',0.4);
            plt(subset) = plot(y,'color',[0.4 0.4 0.4],'linewidth',3);
        elseif subset == 2
            fill([x;flipud(x)],[y-dy;flipud(y+dy)],colors(area,:),'linestyle','none','facealpha',0.4);
            plt(subset) = plot(y,'color',colors2(area,:),'linewidth',3);
        end
        text(0.7,0.875,brain_areas{area},'Units','normalized','FontSize',15)
        
        if area == 4 || area == 1
            ylabel('Ensemble rate correlation')
        elseif area == 5
            xlabel({'Elapsed time (# of movie repeats)';'1 repeat = 30 seconds'})
        end
        
        if area >3
            set(gca,'xtick',[1 10 20 29])
        else
            set(gca,'xtick',[])
        end
        ylim(ylims(area,:))
        xlim([0 30])
        
    end
    legend(plt,{'All cells','Non-suppressed'},'Location','Best');
    legend('boxoff')
end

corrected_pval = bonf_holm(pvalues(2,:));

VarNames = {'area','df','chi','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas(1:6)'),df(2,:)',chi(2,:)',pvalues(2,:)',corrected_pval(:),'VariableNames',VarNames);

clc;
disp(['Ensemble rate correlation as a function of elapsed time for non-suppressed units'])
disp(['Friedman�s tests with Holm�Bonferroni correction:'])
disp(statistics)

%% Figure S2G - Tuning curve correlation on subset of non-adapted units
cell_cutoff = 15;
nat_movie = 1;
num_repeats = 30;
elapsed_tuning_corr_area = {};
for area = 1:6
    valid_mice = neuropixels_cell_count(:,area,nat_movie) >= cell_cutoff & movie_repeats(:,nat_movie) == num_repeats;
    current_area = neuropixels_population_vectors(valid_mice,area,nat_movie);
    
    elapsed_tuning_corr = [];
    elapsed_tuning_corr_subset = [];
    for mouse = 1:size(current_area,1)
        clc;
        disp(['Calculating tuning curve correlation between movie repeats:'])
        disp(['Dataset: Neuropixels | Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse};
        
        pval_across_units = [];
        for unit = 1:size(current_mouse,1)
            current_unit = squeeze(current_mouse(unit,:,:));
            
            current_unit_blockA = nanmean(current_unit(:,1:30));
            current_unit_blockB = nanmean(current_unit(:,31:60));
            
            [~,pval_across_units(unit,1)] = ttest2(current_unit_blockA(:,1:15),current_unit_blockA(:,16:30),'tail','right');
            [~,pval_across_units(unit,2)] = ttest2(current_unit_blockB(:,1:15),current_unit_blockB(:,16:30),'tail','right');
            
        end
        
        valid_units = pval_across_units > 0.05;
        
        if sum(sum(valid_units) == 0) == 0
            
            current_mouse_blockA = current_mouse(:,:,1:30);
            current_mouse_blockB = current_mouse(:,:,31:60);
            
            valid_current_mouse_blockA = current_mouse_blockA(valid_units(:,1),:,:);
            valid_current_mouse_blockB = current_mouse_blockB(valid_units(:,2),:,:);
            
            tuning_corr = [];
            tuning_corr_subset = [];
            for repeat1 = 1:30
                for repeat2 = 1:30
                    tuning_corr(repeat1,repeat2,1) = nanmedian(diag(corr(current_mouse_blockA(:,:,repeat1)',current_mouse_blockA(:,:,repeat2)')));
                    tuning_corr(repeat1,repeat2,2) = nanmedian(diag(corr(current_mouse_blockB(:,:,repeat1)',current_mouse_blockB(:,:,repeat2)')));
                    
                    tuning_corr_subset(repeat1,repeat2,1) = nanmedian(diag(corr(valid_current_mouse_blockA(:,:,repeat1)',valid_current_mouse_blockA(:,:,repeat2)')));
                    tuning_corr_subset(repeat1,repeat2,2) = nanmedian(diag(corr(valid_current_mouse_blockB(:,:,repeat1)',valid_current_mouse_blockB(:,:,repeat2)')));
                end
            end
            tuning_corr = nanmean(tuning_corr,3);
            tuning_corr_subset = nanmean(tuning_corr_subset,3);
            
            for diagonal = 1:29
                elapsed_tuning_corr(mouse,diagonal)  = nanmean(diag(tuning_corr,diagonal));
                elapsed_tuning_corr_subset(mouse,diagonal) = nanmean(diag(tuning_corr_subset,diagonal));
            end
        end
    end
    elapsed_tuning_corr_area(area,:) = {elapsed_tuning_corr,elapsed_tuning_corr_subset};
end


pvalues = [];
df = [];
chi = [];
ylims = [0.28 0.45;0.325 0.475;0.225 0.425;0.25 0.4;0.13 0.3;0.275 0.4];
figure('units','normalized','position',[0.3 0.3 0.35 0.385])
for area = 1:6
    for subset = 1:2
        current_area = elapsed_tuning_corr_area{area,subset};
        
        [pvalues(subset,area),tbl,stats] = friedman(current_area,[1],'off');
        df(subset,area) = tbl{2,3};
        chi(subset,area) = tbl{2,5};
        
        subplot(2,3,area)
        hold on
        mean_elapse_repeat = nanmean(current_area);
        std_elapse_repeat = nanstd(current_area);
        ste_elapse_repeat = std_elapse_repeat./sqrt(size(current_area,1));
        
        x = [1:length(mean_elapse_repeat)]';
        y = mean_elapse_repeat';
        dy = ste_elapse_repeat';
        if subset == 1
            fill([x;flipud(x)],[y-dy;flipud(y+dy)],[0.6 0.6 0.6],'linestyle','none','facealpha',0.4);
            plt(subset) = plot(y,'color',[0.4 0.4 0.4],'linewidth',3);
        elseif subset == 2
            fill([x;flipud(x)],[y-dy;flipud(y+dy)],colors(area,:),'linestyle','none','facealpha',0.4);
            plt(subset) = plot(y,'color',colors2(area,:),'linewidth',3);
        end
        text(0.7,0.875,brain_areas{area},'Units','normalized','FontSize',15)
        
        if area == 4 || area == 1
            ylabel('Tuning curve correlation')
        elseif area == 5
            xlabel({'Elapsed time (# of movie repeats)';'1 repeat = 30 seconds'})
        end
        
        if area >3
            set(gca,'xtick',[1 10 20 29])
        else
            set(gca,'xtick',[])
        end
        ylim(ylims(area,:))
        xlim([0 30])
        
    end
    legend(plt,{'All cells','Non-suppressed'},'Location','Best');
    legend('boxoff')
end

corrected_pval = bonf_holm(pvalues(2,:));

VarNames = {'area','df','chi','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas(1:6)'),df(2,:)',chi(2,:)',pvalues(2,:)',corrected_pval(:),'VariableNames',VarNames);

clc;
disp(['Tuning curve correlation as a function of elapsed time for non-suppressed units'])
disp(['Friedman�s tests with Holm�Bonferroni correction:'])
disp(statistics)

%% Figure S2H - Activity rate difference index for individual cells
cell_cutoff = 15;
nat_movie = 1;
num_repeats = 30;

delta_activity_all_areas = {};
for area = 1:6
    valid_mice = neuropixels_cell_count(:,area,nat_movie) >= cell_cutoff & movie_repeats(:,nat_movie) == num_repeats;
    current_area = neuropixels_population_vectors(valid_mice,area,nat_movie);
    
    delta_activity_all_cells  = [];
    for mouse = 1:size(current_area,1)
        current_mouse = current_area{mouse};
        
        mean_activity_all_cells = [];
        mean_activity_all_cells(:,:,1) = squeeze(nanmean(current_mouse(:,:,1:30),2));
        mean_activity_all_cells(:,:,2) = squeeze(nanmean(current_mouse(:,:,31:60),2));
        
        activity_diff_index = [];
        activity_diff_index(:,1) = [nanmean(mean_activity_all_cells(:,26:30,1),2) - nanmean(mean_activity_all_cells(:,1:5,1),2)]./...
            [nanmean(mean_activity_all_cells(:,26:30,1),2) + nanmean(mean_activity_all_cells(:,1:5,1),2)];
        activity_diff_index(:,2) = [nanmean(mean_activity_all_cells(:,26:30,2),2) - nanmean(mean_activity_all_cells(:,1:5,2),2)]./...
            [nanmean(mean_activity_all_cells(:,26:30,2),2) + nanmean(mean_activity_all_cells(:,1:5,2),2)];
        delta_activity_all_cells = [delta_activity_all_cells;activity_diff_index(:)];
    end
    delta_activity_all_areas{area} = delta_activity_all_cells;
end

delta_activity_mat= nan(3330,6);
for area = 1:6
    delta_activity_mat(1:length(delta_activity_all_areas{area}),area) = delta_activity_all_areas{area};
end

figure('units','normalized','position',[0.3 0.3 0.225 0.35])
hold on
for area = 1:6
    current_area = delta_activity_all_areas{area};
    unit_num_area = length(current_area);
    temp = linspace(2,-2,unit_num_area);
    temp = temp(randperm(unit_num_area));
    scatter(ones([unit_num_area 1])*(area*10)+temp',current_area,20,'filled','MarkerFaceColor',[0.95 0.95 0.95],'MarkerEdgeColor',[0.85 0.85 0.85])
    violin(area*10,current_area,'facealpha',1,'facecolor',[0.8 0.8 0.8],'linecolor',[0.6 0.6 0.6],'linewidth',1.75,'withmdn','true')
    
end
ylim([-1.2 1.2])

ylabel('Activity rate difference index')
xtickangle(90)
set(gca,'xtick',[1:6]*10,'xticklabel',brain_areas(1:6))

%% Figure S2I - responsiveness of four V1 units across movie repeats and blocks
nat_movie = 1;
mouse = 23;
area = 1;
example_mouse = neuropixels_population_vectors{mouse,area,nat_movie}*30;
units_list = [9,18,25,63];

figure('units','normalized','position',[0.3 0.3 0.325 0.375])
for unit = 1:length(units_list)
    current_unit = squeeze(example_mouse(units_list(unit),:,:))';
    smooth_current_unit = [];
    for repeat = 1:size(current_unit,1)
        smooth_current_unit(repeat,:) = imgaussfilt(current_unit(repeat,:),2);
    end
    norm_current_unit = smooth_current_unit ./ max(smooth_current_unit,[],2);
    [row,col] = find(norm_current_unit == 1);
    [B,I] = sort(row);
    subplot(2,4,unit)
    imagesc(norm_current_unit)
    hold on
    plot(xlim, ([size(current_unit,1) , size(current_unit,1)]./2)+0.5,'--',...
        'linewidth',2,'color','w')
    
    
    if unit ==1
        ylabel('Movie repeat')
        text(0.475, 0.925,['Block A'],'Units','normalized','color','w')
        text(0.475, 0.4,['Block B'],'Units','normalized','color','w')
    end
    colormap(newmap3)
    title(['Unit #',num2str(units_list(unit))])
    
    subplot(2,4,unit+4)
    hold on
    mean_current_unit = [];
    mean_current_unit(1,:) = nanmean(current_unit(1:10,:));
    mean_current_unit(2,:) = nanmean(current_unit(11:20,:));
    plot(mean_current_unit(1,:),'-','markersize',20,'linewidth',2,'color',[0.4 0.4 0.4])
    plot(mean_current_unit(2,:)','-','markersize',20,'linewidth',2,'color',[0.8 0.8 0.8])
    [r,p] = corr(mean_current_unit(1,:)',mean_current_unit(2,:)');
    
    title(['r = ',num2str(r)])
    xlim([1 30])
    if unit ==1
        ylabel({'Mean activity rate';'(spike/frame)'})
        legend({'Block A','Block B'},'Location','Best')
        legend('boxoff')
    elseif unit ==2
        xlabel('Time in movie (sec)')
    end
end

%% Figure S2J - Tuning curve correlation between blocks and across V1 units
nat_movie = 1;
area = 1;
mouse = 23;
current_mouse = neuropixels_population_vectors{mouse,area,nat_movie};
mean_activity_blockA = nanmean(current_mouse(:,:,1:10),3);
mean_activity_blockB = nanmean(current_mouse(:,:,11:20),3);
between_blocks_reliability = corr(mean_activity_blockA',mean_activity_blockB');

figure('units','normalized','position',[0.3 0.3 0.22 0.3])
imagesc(between_blocks_reliability)
colormap(newmap3(1:end-10,:))
colorbar
ylabel('Unit ID (Block A)')
xlabel('Unit ID (Block B)')
cb = colorbar;
cb.Label.String = 'Tuning curve correlation';
cb.FontSize = 12;
title('Single animal example:')

%% Figure S2K - Tuning curve correlation of the same units across blocks
nat_movie = 1;
area = 1;
mouse = 23;
current_mouse = neuropixels_population_vectors{mouse,area,nat_movie};
mean_activity_blockA = nanmean(current_mouse(:,:,1:10),3);
mean_activity_blockB = nanmean(current_mouse(:,:,11:20),3);
between_blocks_reliability = corr(mean_activity_blockA',mean_activity_blockB');
tuning_corr = diag(between_blocks_reliability);

figure('units','normalized','position',[0.3 0.3 0.2 0.3])
histogram(tuning_corr,25,'facecolor',[0.8 0.8 0.8],'edgecolor',[0.8 0.8 0.8])
hold on
plot([0.605 0.605],ylim,'--','linewidth',1.5,'color',[0.2 0.2 0.2])
% One arrow from left to right with text on left side
x = [0.575 0.85];    % adjust length and location of arrow
y = [0.87 0.87];      % adjust hieght and width of arrow
annotation('textarrow',x,y,'FontSize',13,'Linewidth',2)
text(0.6,0.975,['Units included'],'Units','normalized')
h=text(0.54,0.3,['Sliding threshold'],'Units','normalized');
set(h,'Rotation',90);
ylim([0 13])
ylabel('Unit count')
xlabel({'Tuning curve correlation';'between blocks'})
xlim([0 1.05])
set(gca,'box','off')

%% Figure S2K - Fraction of cells included in the analysis
repeats = 30;
cell_cutoff = 15;
nat_movie = 1;
reliability_cutoff_list = [-1,0.6,0.7,0.8,0.9];
fraction_of_valid_areas = {};
for area = 1:6
    valid_mice = neuropixels_cell_count(:,area,nat_movie) >= cell_cutoff & movie_repeats(:,nat_movie) == repeats;
    current_area = neuropixels_population_vectors(valid_mice,area,nat_movie);
    
    
    for reliability = 1:length(reliability_cutoff_list)
        reliability_cutoff = reliability_cutoff_list(reliability);
        
        sub = 1;
        fraction_of_valid = [];
        for mouse = 1:size(current_area,1)
            current_mouse = current_area{mouse};
            activity_blockA = nanmean(current_mouse(:,:,1:30),3);
            activity_blockB = nanmean(current_mouse(:,:,31:60),3);
            tuning_corr_between_blocks = diag(corr(activity_blockA',activity_blockB'));
            
            valid_units = tuning_corr_between_blocks >= reliability_cutoff;
            if sum(valid_units) >= cell_cutoff
                fraction_of_valid(sub) = sum(valid_units)./length(valid_units);
                sub = sub + 1;
            end
            
        end
        fraction_of_valid_areas{area,reliability} = fraction_of_valid;
    end
    
end

figure('units','normalized','position',[0.3 0.3 0.315 0.35])
for area = 1:6
    mean_fraction = [];
    ste_fraction = [];
    for reliability = 1:length(reliability_cutoff_list)
        current_area = fraction_of_valid_areas{area,reliability};
        mean_fraction(reliability) = nanmean(current_area);
        std_fraction = nanstd(current_area);
        ste_fraction(reliability) =  std_fraction./sqrt(length(current_area));
    end
    
    subplot(2,3,area)
    hold on
    bar(mean_fraction,'facecolor',[0.85 0.85 0.85],'edgecolor','none')
    errorbar(mean_fraction,ste_fraction,'color',[0.75 0.75 0.75],'linewidth',3,'linestyle','none','capsize',0)
    text(0.7,0.9,brain_areas{area},'Units','normalized','FontSize',12)
    
    set(gca,'xtick',1:5,'xticklabels',{'-1','0.6','0.7','0.8','0.9'})
    if area ==4
        ylabel('Fraction of units included')
    elseif area == 5
        xlabel('Tuning curve correlation threshold')
    end
    ylim([0 1.05])
end

%% Figure S2M - Ensemble rate correlation between repeats within a block including only stable units across blocks
repeats = 30;
cell_cutoff = 15;
nat_movie = 1;
reliability_cutoff_list = [-1,0.6,0.7,0.8,0.9];
elapsed_rate_corr_areas = {};
for area = 1:6
    valid_mice = neuropixels_cell_count(:,area,nat_movie) >= cell_cutoff & movie_repeats(:,nat_movie) == repeats;
    current_area = neuropixels_population_vectors(valid_mice,area,nat_movie);
    
    for reliability = 1:length(reliability_cutoff_list)
        reliability_cutoff = reliability_cutoff_list(reliability);
        
        sub = 1;
        elapsed_rate_corr = [];
        for mouse = 1:size(current_area,1)
            clc;
            disp(['Calculating ensemble rate correlation between movie repeats'])
            disp(['Dataset: Neuropixels | Stimulus: Natural movie 1 | Area: ',brain_areas{area},...
                ' | Threshold: ',num2str(reliability),'/',num2str(length(reliability_cutoff_list)),' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
            
            current_mouse = current_area{mouse};
            activity_blockA = current_mouse(:,:,1:30);
            activity_blockB = current_mouse(:,:,31:60);
            tuning_corr_between_blocks = diag(corr(nanmean(activity_blockA,3)',nanmean(activity_blockB,3)'));
            
            valid_units = tuning_corr_between_blocks >= reliability_cutoff;
            if sum(valid_units) >= cell_cutoff
                valid_activity_blockA = activity_blockA(valid_units,:,:);
                valid_activity_blockB = activity_blockB(valid_units,:,:);
                
                mean_activity_blockA = squeeze(nanmean(valid_activity_blockA,2));
                mean_activity_blockB = squeeze(nanmean(valid_activity_blockB,2));
                
                rate_corr = [];
                rate_corr(:,:,1) = corr(mean_activity_blockA);
                rate_corr(:,:,2) = corr(mean_activity_blockB);
                mean_rate_corr = nanmean(rate_corr,3);
                
                for diagonal = 1:29
                    elapsed_rate_corr(sub,diagonal) = nanmean(diag(mean_rate_corr,diagonal));
                end
                
                sub = sub + 1;
            end
        end
        
        elapsed_rate_corr_areas{area,reliability} = elapsed_rate_corr;
    end
    
end

figure('units','normalized','position',[0.3 0.3 0.315 0.35])
for area = 1:6
    plt = [];
    for reliability = 1:length(reliability_cutoff_list)
        current_area =  elapsed_rate_corr_areas{area,reliability};
        mean_stability = nanmean(current_area);
        mean_stability_norm = mean_stability./mean_stability(1);
        
        subplot(2,3,area)
        hold on
        plt(reliability) = plot(mean_stability_norm,'color',colors(area,:)- 0.1*(5-reliability)*colors(area,:),'linewidth',2);
        text(0.7,0.9,brain_areas{area},'Units','normalized','FontSize',12)
        
        ylim([0.9 1])
        xlim([0 30])
        
        if area ==4
            ylabel('Ensemble rate correlation')
        elseif area == 5
            xlabel({'Elapsed time (# of movie repeats)';'1 repeat = 30 seconds'})
        end
        
        if area > 3
            set(gca,'xtick',[1,10,20,29])
        else
            set(gca,'xtick',[])
        end
        if ~(area == 4 || area == 1)
            set(gca,'ytick',[])
        end
        
        lgd = legend(plt,{'All','0.6','0.7','0.8','0.9'},'Location','Best');
        lgd.FontSize = 8;
        legend('boxoff')
    end
end

%% Figure S2N - Tuning correlation between repeats within a block including only stable units across blocks
repeats = 30;
cell_cutoff = 15;
nat_movie = 1;
reliability_cutoff_list = [-1,0.6,0.7,0.8,0.9];
elapsed_tuning_corr_areas = {};
for area = 1:6
    valid_mice = neuropixels_cell_count(:,area,nat_movie) >= cell_cutoff & movie_repeats(:,nat_movie) == repeats;
    current_area = neuropixels_population_vectors(valid_mice,area,nat_movie);
    
    for reliability = 1:length(reliability_cutoff_list)
        reliability_cutoff = reliability_cutoff_list(reliability);
        
        sub = 1;
        elapsed_tuning_corr = [];
        for mouse = 1:size(current_area,1)
            clc;
            disp(['Calculating tuning curve correlation between movie repeats'])
            disp(['Dataset: Neuropixels | Stimulus: Natural movie 1 | Area: ',brain_areas{area},...
                ' | Threshold: ',num2str(reliability),'/',num2str(length(reliability_cutoff_list)),' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
            
            current_mouse = current_area{mouse};
            activity_blockA = current_mouse(:,:,1:30);
            activity_blockB = current_mouse(:,:,31:60);
            tuning_corr_between_blocks = diag(corr(nanmean(activity_blockA,3)',nanmean(activity_blockB,3)'));
            
            valid_units = tuning_corr_between_blocks >= reliability_cutoff;
            if sum(valid_units) >= cell_cutoff
                valid_activity_blockA = activity_blockA(valid_units,:,:);
                valid_activity_blockB = activity_blockB(valid_units,:,:);
                
                tuning_corr = [];
                for repeat1 = 1:30
                    for repeat2 = 1:30
                        tuning_corr(repeat1,repeat2,1) = nanmedian(diag(corr(valid_activity_blockA(:,:,repeat1)',valid_activity_blockA(:,:,repeat2)')));
                        tuning_corr(repeat1,repeat2,2) = nanmedian(diag(corr(valid_activity_blockB(:,:,repeat1)',valid_activity_blockB(:,:,repeat2)')));
                    end
                end
                mean_tuning_corr = nanmean(tuning_corr,3);
                
                for diagonal = 1:29
                    elapsed_tuning_corr(sub,diagonal) = nanmean(diag(mean_tuning_corr,diagonal));
                end
                
                sub = sub + 1;
            end
        end
        
        elapsed_tuning_corr_areas{area,reliability} = elapsed_tuning_corr;
    end
    
end

figure('units','normalized','position',[0.3 0.3 0.315 0.35])
for area = 1:6
    plt = [];
    for reliability = 1:length(reliability_cutoff_list)
        current_area =  elapsed_tuning_corr_areas{area,reliability};
        mean_stability = nanmean(current_area);
        mean_stability_norm = mean_stability./mean_stability(1);
        
        subplot(2,3,area)
        hold on
        plt(reliability) = plot(mean_stability_norm,'color',colors(area,:)- 0.1*(5-reliability)*colors(area,:),'linewidth',2);
        text(0.7,0.9,brain_areas{area},'Units','normalized','FontSize',12)
        
        ylim([0.715 1])
        xlim([0 30])
        
        if area ==4
            ylabel('Tuning curve correlation')
        elseif area == 5
            xlabel({'Elapsed time (# of movie repeats)';'1 repeat = 30 seconds'})
        end
        
        if area > 3
            set(gca,'xtick',[1,10,20,29])
        else
            set(gca,'xtick',[])
        end
        if ~(area == 4 || area == 1)
            set(gca,'ytick',[])
        end
        
        lgd = legend(plt,{'All','0.6','0.7','0.8','0.9'},'Location','Best');
        lgd.FontSize = 8;
        legend('boxoff')
    end
end

%% Figure S3A - Between blocks stability - PV correlation for area AM across mice

nat_movie = 2; % natural movie 3
area = 6; % area AM
cell_cutoff = 20; % minimun number of cells recorded 
valid_mice = calcium_excitatory_cell_count{area}(:,1) >= cell_cutoff; % including only mice with atleast 20 cells in session A
current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie); % subset mice that passed requiremnets of 'valid_mice'

pv_corr_across_mice = [];
for mouse = 1:length(current_area) % loop over mice
    current_mouse = current_area{mouse}; % subset neuronal activity of a single mouse (#cells x 30 time bins x 10 movie repeats)
    
    current_mouse_blockA =  current_mouse(:,:,1:5); % neuronal activity in block A (repeats 1-5)
    current_mouse_blockA_half1 = mean(current_mouse_blockA(:,:,1:2),3,'omitnan'); % average across movie repeats of first half of block A (repeats 1-2)
    current_mouse_blockA_half2 = mean(current_mouse_blockA(:,:,3:5),3,'omitnan'); % average across movie repeats of second half of block A (repeats 3-5)
    
    current_mouse_blockB =  current_mouse(:,:,6:10);  % neuronal activity in block B (repeats 6-10)
    current_mouse_blockB_half1 = mean(current_mouse_blockB(:,:,1:2),3,'omitnan'); % average across movie repeats of first half of block B (repeats 6-7)
    current_mouse_blockB_half2 = mean(current_mouse_blockB(:,:,3:5),3,'omitnan'); % average across movie repeats of second half of block B (repeats 8-10)
    
    pv_corr_across_mice(:,:,mouse) = corr([current_mouse_blockA_half1,current_mouse_blockA_half2,...
        current_mouse_blockB_half1,current_mouse_blockB_half2]); % calculate the pv correlation between time bins of the four blocks halves
    
end

mean_pv_corr_across_mice = mean(pv_corr_across_mice,3,'omitnan'); % average pv corr matrices across mice
mean_pv_corr_across_mice(boolean(eye(size(mean_pv_corr_across_mice,1)))) = NaN; % for visualization - convert the main diagonal (values of 1 by definition) into NaNs
mean_pv_corr_across_mice(isnan(mean_pv_corr_across_mice)) = max(mean_pv_corr_across_mice(:)); % for visualization - convert the main diagonal into maximal corr value

figure('units','normalized','position',[0.3 0.3 0.275 0.4]) % visualize pv corr between blocks for example area
imagesc(mean_pv_corr_across_mice)
hold on
for line = 1:3
    plot(xlim,[30.5 30.5] +30*(line-1),'color',newmap3(1,:),'linewidth',1.25)
    plot([30.5 30.5] +30*(line-1),ylim,'color',newmap3(1,:),'linewidth',1.25)
end
set(gca,'xtick',15:30:120,'xticklabel',{'1st half','2nd half','1st half','2nd half'},...
    'ytick',15:30:120,'yticklabel',{'1st half','2nd half','1st half','2nd half'})
ytickangle(90)
xlabel('Block A                                           Block B')
ylabel('Block B                                  Block A')

title('Calcium imaging - Natural movie 3')
colormap(newmap3)
cb = colorbar;
cb.Label.String = 'PV correlation';
cb.FontSize = 12;

%% Figure S3B - Ensemble rate correlation within and between blocks of natural movie 3 (calcium imaging)
nat_movie = 2;  % natural movie 3
cell_cutoff = 20; % minimun number of cells recorded 
rate_corr_between_blocks_across_areas = {};
for area = 1:6
    valid_mice = calcium_excitatory_cell_count{area}(:,1) >= cell_cutoff; % including only mice with atleast 20 cells in session A
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie);  % subset mice that passed requiremnets of 'valid_mice'
    
    rate_corr_across_blocks = [];
    for mouse = 1:length(current_area) % loop over mice
        
        clc;
        disp(['Calculating ensemble rate correlation between blocks:'])
        disp(['Dataset: Calcium imaging | Stimulus: Natural movie 3 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse};  % subset neuronal activity of a single mouse (#cells x 30 time bins x 10 movie repeats)
        current_mouse_blockA =  current_mouse(:,:,1:5); % neuronal activity in block A (repeats 1-5)
        current_mouse_blockA_half1 = mean(mean(current_mouse_blockA(:,:,1:2),3,'omitnan'),2,'omitnan'); % average across time bins and movie repeats of first half of block A (repeats 1-2)
        current_mouse_blockA_half2 = mean(mean(current_mouse_blockA(:,:,3:5),3,'omitnan'),2,'omitnan'); % average across time bins and movie repeats of second half of block A (repeats 3-5)
        
        current_mouse_blockB =  current_mouse(:,:,6:10); % neuronal activity in block B (repeats 6-10)
        current_mouse_blockB_half1 = mean(mean(current_mouse_blockB(:,:,1:2),3,'omitnan'),2,'omitnan'); % average across time bins and movie repeats of first half of block B (repeats 6-7)
        current_mouse_blockB_half2 = mean(mean(current_mouse_blockB(:,:,3:5),3,'omitnan'),2,'omitnan'); % average across time bins and movie repeats of second half of block B (repeats 8-10)
        
        rate_corr = corr([current_mouse_blockA_half1,current_mouse_blockA_half2,...
            current_mouse_blockB_half1,current_mouse_blockB_half2]); % calculate the ensemble rate correlation between the four blocks halves
        
        rate_corr_across_blocks(mouse,1) = mean([rate_corr(1,2),rate_corr(3,4)],'omitnan'); % within block similarity
        rate_corr_across_blocks(mouse,2) = mean(mean(rate_corr(1:2,3:4),'omitnan'),'omitnan'); % between blocks similarity
        
    end
    rate_corr_between_blocks_across_areas{area} = rate_corr_across_blocks;
    
end


pvalue = [];
zvalue = [];
plt = [];
figure('Units','Normalized','Position',[0.2 0.4 0.2 0.325]) % visualizing ensemble rate correlation between blocks
for area = 1:6 % loop over areas
    current_area = rate_corr_between_blocks_across_areas{area}; % subset values of a single area
    mean_stability = mean(current_area,'omitnan'); % average across mice
    std_stability = std(current_area,'omitnan'); % standard deviation across mice
    ste_stability = std_stability./sqrt(size(current_area,1)); % standard error across mice
    
    [pvalue(area),~,stats] = signrank(current_area(:,1),current_area(:,2)); % perform two-sided Wilcoxon signed-rank test for difference within and between blocks
    zvalue(area) = stats.zval;
    
    
    hold on
    plt(area) = errorbar(mean_stability,ste_stability,'o','color',colors(area,:),...
        'markerfacecolor',colors(area,:),'capsize',0,'linestyle','-','linewidth',3,'markersize',5);
    
end
xlim([0.5 2.5])
ylim([0.625 0.85])
lgd = legend(plt,brain_areas(1:6));
legend('boxoff')
lgd.Position = [0.2 0.225 0.15 0.15];
set(gca,'xtick',[1,2],'xticklabel',{'Within block','Between blocks'},'ytick',0.65:0.05:0.8)
ylabel('Ensemble rate correlation')

% correct for multiple comparisons using bonferroni-holm method
corrected_pvalue = bonf_holm(pvalue);

% define statistics summary table
VarNames = {'area','zvalue','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas(1:6)'),zvalue(:),pvalue(:),corrected_pvalue(:),'VariableNames',VarNames);

% display statistics summary table
clc;
disp(['Ensemble rate correlation within blocks compared to between blocks'])
disp(['Two-tailed Wilcoxon signed-rank tests with Holm�Bonferroni correction:'])
disp(statistics)

%% Figure S3C - Tuning curve correlation within and between blocks of natural movie 3 (calcium imaging)
nat_movie = 2;  % natural movie 3
cell_cutoff = 20; % minimun number of cells recorded 
tuning_corr_between_blocks_across_areas = {};
for area = 1:6 % loop over areas
    valid_mice = calcium_excitatory_cell_count{area}(:,1) >= cell_cutoff; % including only mice with atleast 20 cells in session A
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie); % subset mice that passed requiremnets of 'valid_mice'
    
    tuning_corr_across_blocks = [];
    for mouse = 1:length(current_area) % loop over mice
        
        clc;
        disp(['Calculating tuning curve correlation between blocks:'])
        disp(['Dataset: Calcium imaging | Stimulus: Natural movie 3 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse}; % subset neuronal activity of a single mouse (#cells x 30 time bins x 10 movie repeats)
        current_mouse_halves = [];
        current_mouse_halves(:,:,1) = mean(current_mouse(:,:,1:2),3,'omitnan'); % average across movie repeats of first half of block A (repeats 1-2)
        current_mouse_halves(:,:,2) = mean(current_mouse(:,:,3:5),3,'omitnan'); % average across movie repeats of second half of block A (repeats 3-5)
        current_mouse_halves(:,:,3) = mean(current_mouse(:,:,6:7),3,'omitnan'); % average across movie repeats of first half of block B (repeats 6-7)
        current_mouse_halves(:,:,4) = mean(current_mouse(:,:,8:10),3,'omitnan'); % average across movie repeats of second half of block B (repeats 8-10)

        tuning_corr = [];
        for half1 = 1:4
            for half2 = 1:4
                tuning_corr(half1,half2) = median(diag(corr(current_mouse_halves(:,:,half1)',current_mouse_halves(:,:,half2)')),'omitnan'); % calculate the tuning curve corrlation across blocks halves and calculate the median across corresponding cells
            end
        end
        
        tuning_corr_across_blocks(mouse,1) = mean([tuning_corr(1,2),tuning_corr(3,4)],'omitnan'); % within block stability
        tuning_corr_across_blocks(mouse,2) = mean(mean(tuning_corr(1:2,3:4),'omitnan'),'omitnan'); % between blocks stability
        
    end
    tuning_corr_between_blocks_across_areas{area} = tuning_corr_across_blocks;
    
end


pvalue = [];
zvalue = [];
plt = [];
figure('Units','Normalized','Position',[0.2 0.4 0.2 0.325]) % visualize tuning curve correlation across blocks
for area = 1:6
    current_area = tuning_corr_between_blocks_across_areas{area}; % subset values of a single area
    mean_stability = mean(current_area,'omitnan'); % average across mice
    std_stability = std(current_area,'omitnan'); % standard deviation across mice
    ste_stability = std_stability./sqrt(size(current_area,1)); % standard error across mice
    
    [pvalue(area),~,stats] = signrank(current_area(:,1),current_area(:,2)); % perform two-sided Wilcoxon signed-rank test for difference within and between blocks
    zvalue(area) = stats.zval;
    
    
    hold on
    plt(area) = errorbar(mean_stability,ste_stability,'o','color',colors(area,:),...
        'markerfacecolor',colors(area,:),'capsize',0,'linestyle','-','linewidth',3,'markersize',5);
    
end
xlim([0.5 2.5])
ylim([0 0.5])
lgd = legend(plt,brain_areas(1:6));
legend('boxoff')
lgd.Position = [0.2 0.225 0.15 0.15];
set(gca,'xtick',[1,2],'xticklabel',{'Within block','Between blocks'},'ytick',0.65:0.05:0.8)
ylabel('Tuning curve correlation')

% correct for multiple comparisons using bonferroni-holm method
corrected_pvalue = bonf_holm(pvalue);

% define statistics summary table
VarNames = {'area','zvalue','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas(1:6)'),zvalue(:),pvalue(:),corrected_pvalue(:),'VariableNames',VarNames);

% display statistics summary table
clc;
disp(['Tuning curve correlation within blocks compared to between blocks'])
disp(['Two-tailed Wilcoxon signed-rank tests with Holm�Bonferroni correction:'])
disp(statistics)

%% Figure S3D - Ensmble rate and tuning curve correlation difference between blocks (calcium imaging)
nat_movie = 2;  % natural movie 3
cell_cutoff = 20; % minimun number of cells recorded 
rate_tuning_corr_diff_across_areas = {};
for area = 1:6 % loop over areas
    valid_mice = calcium_excitatory_cell_count{area}(:,1) >= cell_cutoff; % including only mice with atleast 20 cells in session A
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie); % subset mice that passed requiremnets of 'valid_mice'
    
    rate_corr_diff = [];
    tuning_corr_diff = [];
    for mouse = 1:length(current_area) % loop over mice
        
        clc;
        disp(['Calculating ensemble rate and tuning curve correlation between blocks:'])
        disp(['Dataset: Calcium imaging | Stimulus: Natural movie 3 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
           current_mouse = current_area{mouse}; % subset neuronal activity of a single mouse (#cells x 30 time bins x 10 movie repeats)
        current_mouse_halves = [];
        current_mouse_halves(:,:,1) = mean(current_mouse(:,:,1:2),3,'omitnan'); % average across movie repeats of first half of block A (repeats 1-2)
        current_mouse_halves(:,:,2) = mean(current_mouse(:,:,3:5),3,'omitnan'); % average across movie repeats of second half of block A (repeats 3-5)
        current_mouse_halves(:,:,3) = mean(current_mouse(:,:,6:7),3,'omitnan'); % average across movie repeats of first half of block B (repeats 6-7)
        current_mouse_halves(:,:,4) = mean(current_mouse(:,:,8:10),3,'omitnan'); % average across movie repeats of second half of block B (repeats 8-10)

        rate_corr = [];
        tuning_corr = [];
        for half1 = 1:4
            for half2 = 1:4
                rate_corr(half1,half2) = corr(mean(current_mouse_halves(:,:,half1),2,'omitnan'),mean(current_mouse_halves(:,:,half2),2,'omitnan')); % average activity across time bins and calculate the rate correlation between blocks halves
                tuning_corr(half1,half2) = median(diag(corr(current_mouse_halves(:,:,half1)',current_mouse_halves(:,:,half2)')),'omitnan'); % calculate the tuning curve corrlation across blocks halves and calculate the median across corresponding cells
            end
        end
        
        
        tuning_corr_across_blocks = [];
        tuning_corr_across_blocks(1) = mean([tuning_corr(1,2),tuning_corr(3,4)],'omitnan'); % within block tuning stability
        tuning_corr_across_blocks(2) = mean(mean(tuning_corr(1:2,3:4),'omitnan'),'omitnan'); % between blocks tuning stability
        
        rate_corr_across_blocks = [];
        rate_corr_across_blocks(1) = mean([rate_corr(1,2),rate_corr(3,4)],'omitnan'); % within block rate stability
        rate_corr_across_blocks(2) = mean(mean(rate_corr(1:2,3:4),'omitnan'),'omitnan'); % between blocks rate stability
        
        rate_corr_diff(mouse) = rate_corr_across_blocks(1) - rate_corr_across_blocks(2); % rate corr difference within and between blocks
        tuning_corr_diff(mouse) = tuning_corr_across_blocks(1) - tuning_corr_across_blocks(2); % tuning corr difference within and between blocks
    end 
    
    rate_tuning_corr_diff_across_areas{area} = [rate_corr_diff',tuning_corr_diff'];
    
end


pvalue = [];
zvalue = [];
mean_stability = [];
ste_stability = [];
for area = 1:6 % loop over areas
    current_area = rate_tuning_corr_diff_across_areas{area}; % subset corr values of a single area
    mean_stability(area,:) = mean(current_area,'omitnan'); % average across mice
    std_stability = std(current_area,'omitnan'); % standard deviation across mice
    ste_stability(area,:) = std_stability./sqrt(size(current_area,1)); % standard error across mice
    
    [pvalue(area,1),~,stats] = signrank(current_area(:,1)); % two-sided wilcoxon signed-rank test for rate corr difference within and between blocks
    zvalue(area,1) = stats.zval;
    
    [pvalue(area,2),~,stats] = signrank(current_area(:,2)); % two-sided wilcoxon signed-rank test for tuning corr difference within and between blocks
    zvalue(area,2) = stats.zval;
end

plt = [];
figure('Units','Normalized','Position',[0.3 0.4 0.225 0.325]) % visualize rate and tuning corr stability across blocks
xlim([0 7])
plot(xlim,[0 0],'--','color',[0.2 0.2 0.2],'linewidth',1.5)
hold on
plt(1) = errorbar(mean_stability(:,1),ste_stability(:,1),'o','color',[0.7 0.7 0.7],...
    'markerfacecolor',[0.7 0.7 0.7],'capsize',0,'linestyle','none','linewidth',3,'markersize',5);
plt(2) = errorbar(mean_stability(:,2),ste_stability(:,2),'o','color',[0.5 0.5 0.5],...
    'markerfacecolor',[0.5 0.5 0.5],'capsize',0,'linestyle','none','linewidth',3,'markersize',5);

lgd = legend(plt,{'Ensemble rate correlation','Tuning curve correlation'},'Location','southwest');
legend('boxoff')
ylim([-0.04 0.12])
set(gca,'xtick',1:6,'xticklabel',brain_areas(1:6),'box','off')
ylabel({'Correlation difference';'(within block - between blocks)'})

% correct for multiple comparisons using bonferroni-holm method
corrected_pvalue = [];
corrected_pvalue(:,1) = bonf_holm(pvalue(:,1)); % for rate corr
corrected_pvalue(:,2) = bonf_holm(pvalue(:,2)); % for tuning corr

% define statistics summary table
VarNames = {'area','Rate_zvalue','Rate_pvalue','Rate_bonf_holm','Tuning_zvalue','Tuning_pvalue','Tuning_bonf_holm'};
statistics = table(cell2mat(brain_areas(1:6)'),zvalue(:,1),pvalue(:,1),corrected_pvalue(:,1),...
    zvalue(:,2),pvalue(:,2),corrected_pvalue(:,2),'VariableNames',VarNames);


% display statistics summary table
clc;
disp(['Ensemble rate and tuning curve correlation within blocks compared to between blocks'])
disp(['Two-tailed Wilcoxon signed-rank tests with Holm�Bonferroni correction:'])
disp(statistics)

%% Figure S3E - Average ensemble rate correlation across V1 animals in the Brain Observatory group
cell_cutoff = 15; % minimum number off units recorded

reliability_cutoff = 0.6; % treshold for tuning curve correlation between blocks
area = 1; % example are V1
valid_mice = neuropixels_cell_count(:,area) >= cell_cutoff & movie_repeats(:,1) == 10; % include only mice from the 'Brain observatory' group with atleast 15 units
current_area = squeeze(neuropixels_population_vectors(valid_mice,area,:)); % subset only mice that met the requirements of 'valid_mice'
across_blocks_ensembles_corr = [];
sub = 1;
for mouse = 1:length(current_area) % loop over mice
    current_mouse = current_area(mouse,:); % subset neuronal activity of both movies a single mouse (each matrix of the size #units x 30 time bins x # num repeats)
    
    movie1_blockA = current_mouse{1}(:,:,1:10); % activity during block A of natural movie 1 (repeats 1-10)
    movie1_blockB = current_mouse{1}(:,:,11:20); % activity during block B of natural movie 1 (repeats 11-20)
    
    movie3_blockA = current_mouse{2}(:,:,1:5); % activity during block a of natural movie 3 (repeats 1-5)
    movie3_blockB = current_mouse{2}(:,:,6:10); % activity during block B of natural movie 3 (repeats 6-10)
    
    mean_movie1_blockA = mean(movie1_blockA,3,'omitnan'); % average activity across block of  block A of natural movie 1
    mean_movie1_blockB = mean(movie1_blockB,3,'omitnan'); % average activity across block of  block B of natural movie 1
    
    mean_movie3_blockA = mean(movie3_blockA,3,'omitnan'); % average activity across block of  block A of natural movie 3
    mean_movie3_blockB = mean(movie3_blockB,3,'omitnan'); % average activity across block of  block B of natural movie 3
    
    
    movie1_blocks_stability = diag(corr(mean_movie1_blockA',mean_movie1_blockB')); % tuning curve corr of corresponding units between blocks of natural movie 1
    movie3_blocks_stability = diag(corr(mean_movie3_blockA',mean_movie3_blockB')); % tuning curve corr of corresponding units between blocks of natural movie 3

    % 'valid_units_both_blocks' - include cells with tuning curve corr
    % between blocks of atleast 0.6 for both natual movie 1 and natural movie 3
    valid_units_both_blocks = (movie1_blocks_stability >= reliability_cutoff) & (movie3_blocks_stability >= reliability_cutoff);
    if sum(valid_units_both_blocks) >= cell_cutoff % include only mice with atleast 15 valid units
        
        % subset the neuronal activity of valid units 
        valid_movie1_blockA_stability = movie1_blockA(valid_units_both_blocks,:,:); 
        valid_movie1_blockB_stability = movie1_blockB(valid_units_both_blocks,:,:);
        valid_movie3_blockA_stability = movie3_blockA(valid_units_both_blocks,:,:);
        valid_movie3_blockB_stability = movie3_blockB(valid_units_both_blocks,:,:);
        
        % average the neuronal activity across time bins for each unit
        mean_valid_movie1_blockA_stability = squeeze(mean(valid_movie1_blockA_stability,2,'omitnan'));
        mean_valid_movie1_blockB_stability = squeeze(mean(valid_movie1_blockB_stability,2,'omitnan'));
        mean_valid_movie3_blockA_stability = squeeze(mean(valid_movie3_blockA_stability,2,'omitnan'));
        mean_valid_movie3_blockB_stability = squeeze(mean(valid_movie3_blockB_stability,2,'omitnan'));
        
        
        mean_valid_movie1_half_blockA_stability = [mean(mean_valid_movie1_blockA_stability(:,1:5),2,'omitnan'),...
            mean(mean_valid_movie1_blockA_stability(:,6:10),2,'omitnan')]; % average the neruonal activity across movie repeats for each half of block A of natural movie 1
        mean_valid_movie1_half_blockB_stability = [mean(mean_valid_movie1_blockB_stability(:,1:5),2,'omitnan'),...
            mean(mean_valid_movie1_blockB_stability(:,6:10),2,'omitnan')]; % average the neruonal activity across movie repeats for each half of block B of natural movie 1
        
        mean_valid_movie3_half_blockA_stability = [mean(mean_valid_movie3_blockA_stability(:,1:2),2,'omitnan'),...
            mean(mean_valid_movie3_blockA_stability(:,3:5),2)]; % average the neruonal activity across movie repeats for each half of block A of natural movie 3
        mean_valid_movie3_half_blockB_stability = [mean(mean_valid_movie3_blockB_stability(:,1:2),2,'omitnan'),...
            mean(mean_valid_movie3_blockB_stability(:,3:5),2,'omitnan')]; % average the neruonal activity across movie repeats for each half of block B of natural movie 3

        across_blocks_ensembles = [mean_valid_movie3_half_blockA_stability,mean_valid_movie1_half_blockA_stability,...
            mean_valid_movie3_half_blockB_stability,mean_valid_movie1_half_blockB_stability]; % sort activity of block halves based on their presentation time
        
        across_blocks_ensembles_corr(:,:,sub) = corr(across_blocks_ensembles); % ensemble rate correlation between blocks halves of different natural movies 
        
        sub = sub + 1;
    end
end
% visualize the ensemble rate corr across blocks of different natural movie
% for example area
figure('units','normalized','position',[0.3 0.3 0.25 0.35])
imagesc(mean(across_blocks_ensembles_corr,3,'omitnan')) % average matrices across mice
hold on
lines = [2.5:2:7];
for line = 1:3
    plot(xlim,[lines(line) lines(line)],'color',newmap3(1,:),'linewidth',2)
    plot([lines(line) lines(line)],ylim,'color',newmap3(1,:),'linewidth',2)
end
set(gca,'xtick',1.5:2:8,'xticklabels',{'NM3','NM1','NM3','NM1'},...
    'ytick',1.5:2:8,'yticklabels',{'NM3','NM1','NM3','NM1'})
ytickangle(90)
cb = colorbar;
cb.Label.String = 'Ensemble rate correlation';
cb.FontSize = 12;
colormap(newmap3)
title({'Neuropixels dataset';"'Brain Observatory' group:"})


%% Figure S3F - Units included in the analysis in the 'Brain Observatory' group
cell_cutoff = 15; % minimum number off units recorded

reliability_cutoff = 0.6; % treshold for tuning curve correlation between blocks
area = 1; % example are V1
valid_mice = neuropixels_cell_count(:,area) >= cell_cutoff & movie_repeats(:,1) == 10;  % include only mice from the 'Brain observatory' group with atleast 15 units
current_area = squeeze(neuropixels_population_vectors(valid_mice,area,:)); % subset only mice that met the requirements of 'valid_mice'
movie1_blocks_stability_all_units = [];
movie3_blocks_stability_all_units = [];

for mouse = 1:length(current_area) % loop over mice
    current_mouse = current_area(mouse,:); % subset neuronal activity of both movies a single mouse (each matrix of the size #units x 30 time bins x # num repeats)
    
   
    movie1_blockA = current_mouse{1}(:,:,1:10); % activity during block A of natural movie 1 (repeats 1-10)
    movie1_blockB = current_mouse{1}(:,:,11:20); % activity during block B of natural movie 1 (repeats 11-20)
    
    movie3_blockA = current_mouse{2}(:,:,1:5); % activity during block a of natural movie 3 (repeats 1-5)
    movie3_blockB = current_mouse{2}(:,:,6:10); % activity during block B of natural movie 3 (repeats 6-10)
    
    mean_movie1_blockA = mean(movie1_blockA,3,'omitnan'); % average activity across block of  block A of natural movie 1
    mean_movie1_blockB = mean(movie1_blockB,3,'omitnan'); % average activity across block of  block B of natural movie 1
    
    mean_movie3_blockA = mean(movie3_blockA,3,'omitnan'); % average activity across block of  block A of natural movie 3
    mean_movie3_blockB = mean(movie3_blockB,3,'omitnan'); % average activity across block of  block B of natural movie 3
    
    
    
    movie1_blocks_stability = diag(corr(mean_movie1_blockA',mean_movie1_blockB')); % tuning curve corr of corresponding units between blocks of natural movie 1
    movie3_blocks_stability = diag(corr(mean_movie3_blockA',mean_movie3_blockB')); % tuning curve corr of corresponding units between blocks of natural movie 3

    movie1_blocks_stability_all_units = [movie1_blocks_stability_all_units;movie1_blocks_stability]; % store values for natural movie 1
    movie3_blocks_stability_all_units = [movie3_blocks_stability_all_units;movie3_blocks_stability]; % store values for natural movie 3
end

% visualize tuning curve stability across block of natural movie 1 and
% natural movie 3
figure('units','normalized','position',[0.3 0.3 0.3 0.4])
ylim([-1 1])
xlim([-1 1])
hold on
dscatter(movie1_blocks_stability_all_units,movie3_blocks_stability_all_units)
plot(xlim,[0.6 0.6],'--','color',[0.7 0.7 0.7],'linewidth',2.5)
plot([0.6 0.6],ylim,'--','color',[0.7 0.7 0.7],'linewidth',2.5)
h=text(0.75,0.025,['Threshold'],'Units','normalized','fontsize',12);
set(h,'Rotation',90);
text(0.025,0.85,['Threshold'],'Units','normalized','fontsize',12)
ylabel({'Tuning curve corrlation';"between 'Natural movie 3' blocks"})
xlabel({'Tuning curve corrlation';"between 'Natural movie 1' blocks"})
cb = colorbar;
set(cb,'xtick',[0.01 0.995],'xticklabel',{'min','max'})
cb.Label.String = 'Density';
cb.FontSize = 12;
colormap(magma_colormap)
%% Fiugre S3G - Ensemble rate correlation across same and different movie blocks in the 'Brain Observatory' group

cell_cutoff = 15;% minimum number off units recorded
reliability_cutoff = 0.6;% treshold for tuning curve correlation between blocks

diff_maps_areas = {};
same_maps_areas = {};
for area = 1:6 % loop over mice
    valid_mice = (movie_repeats(:,1) == 10) &(neuropixels_cell_count(:,area,1) >= cell_cutoff); % include only mice from the 'Brain observatory' group with atleast 15 units
    current_area = squeeze(neuropixels_population_vectors(valid_mice,area,:)); % subset only mice that met the requirements of 'valid_mice'
 
    diff_maps = [];
    same_maps = [];
    sub = 1;
    for mouse = 1:length(current_area) % loop over mice
        clc;
        disp(['Calculating ensemble rate correlation between blocks:'])
        disp(['Dataset: Neuropixels | Stimulus: Natural movie 1 & Natural movie 3 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        current_mouse = current_area(mouse,:); % subset neuronal activity of both movies a single mouse (each matrix of the size #units x 30 time bins x # num repeats)
    
    movie1_blockA = current_mouse{1}(:,:,1:10); % activity during block A of natural movie 1 (repeats 1-10)
    movie1_blockB = current_mouse{1}(:,:,11:20); % activity during block B of natural movie 1 (repeats 11-20)
    
    movie3_blockA = current_mouse{2}(:,:,1:5); % activity during block a of natural movie 3 (repeats 1-5)
    movie3_blockB = current_mouse{2}(:,:,6:10); % activity during block B of natural movie 3 (repeats 6-10)
    
    mean_movie1_blockA = mean(movie1_blockA,3,'omitnan'); % average activity across block of  block A of natural movie 1
    mean_movie1_blockB = mean(movie1_blockB,3,'omitnan'); % average activity across block of  block B of natural movie 1
    
    mean_movie3_blockA = mean(movie3_blockA,3,'omitnan'); % average activity across block of  block A of natural movie 3
    mean_movie3_blockB = mean(movie3_blockB,3,'omitnan'); % average activity across block of  block B of natural movie 3
    
    
    movie1_blocks_stability = diag(corr(mean_movie1_blockA',mean_movie1_blockB')); % tuning curve corr of corresponding units between blocks of natural movie 1
    movie3_blocks_stability = diag(corr(mean_movie3_blockA',mean_movie3_blockB')); % tuning curve corr of corresponding units between blocks of natural movie 3

    % 'valid_units_both_blocks' - include cells with tuning curve corr
    % between blocks of atleast 0.6 for both natual movie 1 and natural movie 3
        valid_units_both_blocks = (movie1_blocks_stability >= reliability_cutoff) & (movie3_blocks_stability >= reliability_cutoff);
        if sum(valid_units_both_blocks) >= cell_cutoff % include only mice with atleast 15 valid units
            
           % subset the neuronal activity of valid units 
        valid_movie1_blockA_stability = movie1_blockA(valid_units_both_blocks,:,:); 
        valid_movie1_blockB_stability = movie1_blockB(valid_units_both_blocks,:,:);
        valid_movie3_blockA_stability = movie3_blockA(valid_units_both_blocks,:,:);
        valid_movie3_blockB_stability = movie3_blockB(valid_units_both_blocks,:,:);
        
        % average the neuronal activity across time bins for each unit
        mean_valid_movie1_blockA_stability = squeeze(mean(valid_movie1_blockA_stability,2,'omitnan'));
        mean_valid_movie1_blockB_stability = squeeze(mean(valid_movie1_blockB_stability,2,'omitnan'));
        mean_valid_movie3_blockA_stability = squeeze(mean(valid_movie3_blockA_stability,2,'omitnan'));
        mean_valid_movie3_blockB_stability = squeeze(mean(valid_movie3_blockB_stability,2,'omitnan'));

        mean_valid_movie1_half_blockA_stability = [mean(mean_valid_movie1_blockA_stability(:,1:5),2,'omitnan'),...
            mean(mean_valid_movie1_blockA_stability(:,6:10),2,'omitnan')]; % average the neruonal activity across movie repeats for each half of block A of natural movie 1
        mean_valid_movie1_half_blockB_stability = [mean(mean_valid_movie1_blockB_stability(:,1:5),2,'omitnan'),...
            mean(mean_valid_movie1_blockB_stability(:,6:10),2,'omitnan')]; % average the neruonal activity across movie repeats for each half of block B of natural movie 1
        
        mean_valid_movie3_half_blockA_stability = [mean(mean_valid_movie3_blockA_stability(:,1:2),2,'omitnan'),...
            mean(mean_valid_movie3_blockA_stability(:,3:5),2)]; % average the neruonal activity across movie repeats for each half of block A of natural movie 3
        mean_valid_movie3_half_blockB_stability = [mean(mean_valid_movie3_blockB_stability(:,1:2),2,'omitnan'),...
            mean(mean_valid_movie3_blockB_stability(:,3:5),2,'omitnan')]; % average the neruonal activity across movie repeats for each half of block B of natural movie 3

            
            across_blocks_ensembles = [mean_valid_movie3_half_blockA_stability,mean_valid_movie1_half_blockA_stability,...
                mean_valid_movie3_half_blockB_stability,mean_valid_movie1_half_blockB_stability]; % sort activity of block halves based on their presentation time
            
            across_blocks_ensembles_corr = corr(across_blocks_ensembles); % ensemble rate correlation between blocks halves of different natural movies 
            
            % cauculate mean ensemble rate across blocks halves of different movie blocks
            mean_across_blocks_ensembles_corr = [];
            for half1 = 1:4 % loop over blocks
                rows = [1:2] +2*(half1-1);
                for half2 = 1:4 % loop over blocks
                    cols = [1:2] +2*(half2-1);
                    current_half = across_blocks_ensembles_corr(rows,cols); % subset rate corr between blocks halves
                    mean_across_blocks_ensembles_corr(half1,half2) = mean(current_half(:),'omitnan'); % average rate corr between blocks halves
                end
            end
           
            diff_maps(sub,:) = [diag(mean_across_blocks_ensembles_corr,1)',mean_across_blocks_ensembles_corr(1,4)]; % mean rate corr across blocks of different movies
            same_maps(sub,:) = [mean_across_blocks_ensembles_corr(1,3),mean_across_blocks_ensembles_corr(2,4)]; % mean rate corr across blocks of the same movie
            
            sub = sub +1;
            
        end
    end
    diff_maps_areas{area} = diff_maps;
    same_maps_areas{area} = same_maps;
end

same_map_ticks = [20 72]; % interval in minutes between blocks of the same movie
diff_map_ticks = [0 15 47 77]; % interval in minutes between blocks of different movies
ylims=[0.85 1;0.685 1;0.85 1;0.75 1;0.825 1;0.825 1];
plt = [];
figure('units','normalized','position',[0.3 0.3 0.4 0.45]) % visualize the ensemble rate corr across block of different movies as function of elapsed time
for area = 1:6 % loop over areas
    current_area_diff_maps = diff_maps_areas{area}; % subset the rate corr across blocks of different movies
    current_area_same_maps = same_maps_areas{area}; % subset the rate corr across blocks of the same movies
    
    mean_diff_maps = mean(current_area_diff_maps,'omitnan'); % average across mice
    ste_diff_maps = std(current_area_diff_maps,'omitnan')./sqrt(size(current_area_diff_maps,1)); % standard error across mice
    
    mean_same_maps = mean(current_area_same_maps,'omitnan'); % average across mice
    ste_same_maps = std(current_area_same_maps,'omitnan')./sqrt(size(current_area_same_maps,1)); % standard error across mice
    
    subplot(2,3,area)
    hold on
    
    plt(1)=errorbar(diff_map_ticks,mean_diff_maps,ste_diff_maps,'o','markerfacecolor',[0.6 0.6 0.6],'linestyle','-',...
        'linewidth',2.5,'color',[0.6 0.6 0.6],'CapSize',0);
    
    plt(2)=errorbar(same_map_ticks,mean_same_maps,ste_same_maps,'o','markerfacecolor',[0.3 0.3 0.3],'linestyle','-',...
        'linewidth',2.5,'color',[0.3 0.3 0.3],'CapSize',0);
    text(0.7,0.9,brain_areas{area},'Units','normalized','FontSize',15)
    text(0.65,0.8,['N=',num2str(size(current_area_same_maps,1))],'Units','normalized','FontSize',13)
    
    xlim([-9 86])
    ylim(ylims(area,:))
    if area == 4 || area ==1
        ylabel('Ensemble rate correlation')
    elseif area == 5
        xlabel('Time between blocks (min)')
    end
    if area <4
        set(gca,'xtick',[])
    end
    
    legend(plt,{'Different','Same'},'Location','southwest')
    legend('boxoff')
end

%% Figure S3H - Average ensemble rate correlation across V1 animals in the Functional connectivity group

cell_cutoff = 15; % minimum number off units recorded
reliability_cutoff = 0.5; % treshold for tuning curve correlation between blocks
area = 1; % example area V1
valid_mice = neuropixels_cell_count(:,area) >= cell_cutoff & movie_repeats(:,1) == 30; % include only mice from the 'Functional connectivity' group with atleast 15 units
current_area = squeeze(neuropixels_population_vectors(valid_mice,area,:)); % subset only mice that met the requirements of 'valid_mice'

across_blocks_ensembles_corr = [];
sub = 1;
for mouse = 1:length(current_area) % loop over mice
    current_mouse = current_area(mouse,:); % subset neuronal activity of both movies a single mouse (each matrix of the size #units x 30 time bins x # num repeats)
    nm1_blockA = current_mouse{1}(:,:,1:30); % activity during block A of natural movie 1 (repeats 1-30)
    nm1_blockB = current_mouse{1}(:,:,31:60); % activity during block B of natural movie 1 (repeats 31-60)
    
    snm1_blockA = current_mouse{2}(:,:,1:10); % activity during block A of shuffled natural movie 1 (repeats 1-10)
    snm1_blockB = current_mouse{2}(:,:,11:20); % activity during block B of shuffled natural movie 1 (repeats 11-20)
    
    mean_nm1_blockA = mean(nm1_blockA,3,'omitnan'); % average activity across block of  block A of  natural movie 1
    mean_nm1_blockB = mean(nm1_blockB,3,'omitnan'); % average activity across block of  block B of  natural movie 1
    
    mean_snm1_blockA = mean(snm1_blockA,3,'omitnan'); % average activity across block of  block A of shuffled natural movie 1
    mean_snm1_blockB = mean(snm1_blockB,3,'omitnan'); % average activity across block of  block A of shuffled natural movie 1
    
    
    nm1_blocks_stability = diag(corr(mean_nm1_blockA',mean_nm1_blockB')); % tuning curve corr of corresponding units between blocks of natural movie 1
    snm1_blocks_stability = diag(corr(mean_snm1_blockA',mean_snm1_blockB')); % tuning curve corr of corresponding units between blocks of shuffled natural movie 1
    
    % 'valid_units_both_blocks' - include cells with tuning curve corr
    % between blocks of atleast 0.5 for both natual movie 1 
    valid_units_both_blocks = (nm1_blocks_stability >= reliability_cutoff);
    if sum(valid_units_both_blocks) >= cell_cutoff % include only mice with atleast 15 valid units
        
        % subset the neuronal activity of valid units 
        valid_nm1_blockA_stability = nm1_blockA(valid_units_both_blocks,:,:);
        valid_nm1_blockB_stability = nm1_blockB(valid_units_both_blocks,:,:);
        valid_snm1_blockA_stability = snm1_blockA(valid_units_both_blocks,:,:);
        valid_snm1_blockB_stability = snm1_blockB(valid_units_both_blocks,:,:);
        
        % average the neuronal activity across time bins for each unit
        mean_valid_nm1_blockA_stability = squeeze(mean(valid_nm1_blockA_stability,2,'omitnan'));
        mean_valid_nm1_blockB_stability = squeeze(mean(valid_nm1_blockB_stability,2,'omitnan'));
        mean_valid_snm1_blockA_stability = squeeze(mean(valid_snm1_blockA_stability,2,'omitnan'));
        mean_valid_snm1_blockB_stability = squeeze(mean(valid_snm1_blockB_stability,2,'omitnan'));
        
         % taking the first 10 movie repeats in each block to control for
         % the different number of movie repeats between movies
        mean_valid_nm1_half_blockA_stability = [mean(mean_valid_nm1_blockA_stability(:,1:5),2,'omitnan'),...
            mean(mean_valid_nm1_blockA_stability(:,6:10),2,'omitnan')];  % average the neruonal activity across movie repeats for each half of block A of natural movie 1
        mean_valid_nm1_half_blockB_stability = [mean(mean_valid_nm1_blockB_stability(:,1:5),2,'omitnan'),...
            mean(mean_valid_nm1_blockB_stability(:,6:10),2,'omitnan')];  % average the neruonal activity across movie repeats for each half of block B of natural movie 1
        
        mean_valid_snm1_half_blockA_stability = [mean(mean_valid_snm1_blockA_stability(:,1:5),2,'omitnan'),...
            mean(mean_valid_snm1_blockA_stability(:,6:10),2,'omitnan')];  % average the neruonal activity across movie repeats for each half of block A of shuffled natural movie 1
        mean_valid_snm1_half_blockB_stability = [mean(mean_valid_snm1_blockB_stability(:,1:5),2,'omitnan'),...
            mean(mean_valid_snm1_blockB_stability(:,6:10),2,'omitnan')];  % average the neruonal activity across movie repeats for each half of block B of shuffled natural movie 1

        across_blocks_ensembles = [mean_valid_nm1_half_blockA_stability,mean_valid_snm1_half_blockA_stability...
            mean_valid_snm1_half_blockB_stability,mean_valid_nm1_half_blockB_stability]; % sort activity of block halves based on their presentation time
        
        
        across_blocks_ensembles_corr(:,:,sub) = corr(across_blocks_ensembles); % ensemble rate correlation between blocks halves of different natural movies 
        sub = sub + 1;
    end
end


% visualize the ensemble rate corr across blocks of different natural movie
% for example area
figure
imagesc(mean(across_blocks_ensembles_corr,3,'omitnan')) % average matrices across mice
hold on
lines = [2.5:2:7];
for line = 1:3
    plot(xlim,[lines(line) lines(line)],'color',newmap3(1,:),'linewidth',2)
    plot([lines(line) lines(line)],ylim,'color',newmap3(1,:),'linewidth',2)
end
set(gca,'xtick',1.5:2:8,'xticklabels',{'NM1','SNM1','SNM1','NM1'},...
    'ytick',1.5:2:8,'yticklabels',{'NM1','SNM1','SNM1','NM1'})
ytickangle(90)
cb = colorbar;
cb.Label.String = 'Ensemble rate correlation';
cb.FontSize = 12;
colormap(newmap3)
title({'Neuropixels dataset';"'Brain Observatory' group:"})

%% Figure S3I - Units included in the analysis in the 'Functional Connectivity' group
cell_cutoff = 15; % minimum number off units recorded
area = 1; % example area V1
valid_mice = neuropixels_cell_count(:,area,1) >= cell_cutoff & movie_repeats(:,1) == 30; % include only mice from the 'Functional connectivity' group with atleast 15 units
current_area = squeeze(neuropixels_population_vectors(valid_mice,area,:)); % subset only mice that met the requirements of 'valid_mice'

movie1_blocks_stability = [];
for mouse = 1:length(current_area) % loop over mice
    current_mouse = current_area(mouse,:); % subset neuronal activity of both movies a single mouse (each matrix of the size #units x 30 time bins x # num repeats)
    movie1_blockA = current_mouse{1}(:,:,1:30); % activity during block A of natural movie 1 (repeats 1-30)
    movie1_blockB = current_mouse{1}(:,:,31:60); % activity during block B of natural movie 1 (repeats 31-60)
    mean_movie1_blockA = mean(movie1_blockA,3,'omitnan'); % average activity across block of  block A of  natural movie 1
    mean_movie1_blockB = mean(movie1_blockB,3,'omitnan'); % average activity across block of  block B of  natural movie 1
    
    movie1_blocks_stability = [movie1_blocks_stability;diag(corr(mean_movie1_blockA',mean_movie1_blockB'))]; % tuning curve corr of corresponding units between blocks of natural movie 1
end

figure('units','normalized','position',[0.3 0.3 0.275 0.4]) % visualize tuning curve corr distribution for natural movie 1
hold on
histogram(movie1_blocks_stability,40,'facecolor',[0.8 0.8 0.8],'edgecolor',[0.7 0.7 0.7])
xlim([-1 1])
plot([0.505 0.505],ylim,'--','linewidth',1.5,'color',[0.2 0.2 0.2])
% One arrow from left to right with text on left side
x = [0.71 0.9];    % adjust length and location of arrow
y = [0.87 0.87];      % adjust hieght and width of arrow
annotation('textarrow',x,y,'FontSize',13,'Linewidth',2)
text(0.775,0.975,['Units included'],'Units','normalized')
h=text(0.725,0.3,['Threshold'],'Units','normalized');
set(h,'Rotation',90);
ylabel('Unit count')
xlabel({'Tuning curve correlation';"between 'Natural movie 1' blocks"})
ylim([0 270])

%% Figure S3J - Ensemble rate correlation across same and different movie blocks in the 'Functional connectivity' group
cell_cutoff = 15; % minimum number off units recorded
reliability_cutoff = 0.5; % treshold for tuning curve correlation between blocks

diff_maps_areas = {};
same_maps_areas = {};
for area = 1:6 % loop over areas
    valid_mice = neuropixels_cell_count(:,area,1) >= cell_cutoff & movie_repeats(:,1) == 30;  % include only mice from the 'Functional connectivity' group with atleast 15 units
    current_area = squeeze(neuropixels_population_vectors(valid_mice,area,:)); % subset only mice that met the requirements of 'valid_mice'
    
    sub = 1;
    diff_maps = [];
    same_maps = [];
    for mouse = 1:length(current_area) % loop over mice
        clc;
        disp(['Calculating ensemble rate correlation between blocks:'])
        disp(['Dataset: Neuropixels | Stimulus: Natural movie 1 & Shuffled natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area(mouse,:); % subset neuronal activity of both movies a single mouse (each matrix of the size #units x 30 time bins x # num repeats)
        nm1_blockA = current_mouse{1}(:,:,1:30); % activity during block A of natural movie 1 (repeats 1-30)
        nm1_blockB = current_mouse{1}(:,:,31:60); % activity during block B of natural movie 1 (repeats 31-60)
        
        snm1_blockA = current_mouse{2}(:,:,1:10); % activity during block A of shuffled natural movie 1 (repeats 1-10)
        snm1_blockB = current_mouse{2}(:,:,11:20); % activity during block B of shuffled natural movie 1 (repeats 11-20)
        
        mean_nm1_blockA = mean(nm1_blockA,3,'omitnan'); % average activity across block of  block A of  natural movie 1
        mean_nm1_blockB = mean(nm1_blockB,3,'omitnan'); % average activity across block of  block B of  natural movie 1
        
        nm1_blocks_stability = diag(corr(mean_nm1_blockA',mean_nm1_blockB')); % tuning curve corr of corresponding units between blocks of natural movie 1
      
          % 'valid_units_both_blocks' - include cells with tuning curve corr
            % between blocks of atleast 0.5 for natual movie 1
        valid_units_both_blocks = nm1_blocks_stability >= reliability_cutoff;
        if sum(valid_units_both_blocks) >= cell_cutoff % include only mice with atleast 15 valid units
           
            % subset the neuronal activity of valid units 
            valid_nm1_blockA_stability = nm1_blockA(valid_units_both_blocks,:,:);
            valid_nm1_blockB_stability = nm1_blockB(valid_units_both_blocks,:,:);
            valid_snm1_blockA_stability = snm1_blockA(valid_units_both_blocks,:,:);
            valid_snm1_blockB_stability = snm1_blockB(valid_units_both_blocks,:,:);
            
            % average the neuronal activity across time bins for each unit
            mean_valid_nm1_blockA_stability = squeeze(mean(valid_nm1_blockA_stability,2,'omitnan'));
            mean_valid_nm1_blockB_stability = squeeze(mean(valid_nm1_blockB_stability,2,'omitnan'));
            mean_valid_snm1_blockA_stability = squeeze(mean(valid_snm1_blockA_stability,2,'omitnan'));
            mean_valid_snm1_blockB_stability = squeeze(mean(valid_snm1_blockB_stability,2,'omitnan'));
            
            mean_valid_nm1_half_blockA_stability = [mean(mean_valid_nm1_blockA_stability(:,1:5),2,'omitnan'),...
               mean(mean_valid_nm1_blockA_stability(:,6:10),2,'omitnan')];  % average the neruonal activity across movie repeats for each half of block A of natural movie 1
            mean_valid_nm1_half_blockB_stability = [mean(mean_valid_nm1_blockB_stability(:,1:5),2,'omitnan'),...
                mean(mean_valid_nm1_blockB_stability(:,6:10),2,'omitnan')]; % average the neruonal activity across movie repeats for each half of block B of natural movie 1
            
            mean_valid_snm1_half_blockA_stability = [mean(mean_valid_snm1_blockA_stability(:,1:5),2,'omitnan'),...
                mean(mean_valid_snm1_blockA_stability(:,6:10),2,'omitnan')];  % average the neruonal activity across movie repeats for each half of block A of shuffled natural movie 1
            mean_valid_snm1_half_blockB_stability = [mean(mean_valid_snm1_blockB_stability(:,1:5),2,'omitnan'),...
                mean(mean_valid_snm1_blockB_stability(:,6:10),2,'omitnan')];  % average the neruonal activity across movie repeats for each half of block B of shuffled natural movie 1
            

            across_blocks_ensembles = [mean_valid_nm1_half_blockA_stability,mean_valid_snm1_half_blockA_stability...
                mean_valid_snm1_half_blockB_stability,mean_valid_nm1_half_blockB_stability]; % sort activity of block halves based on their presentation time
            
            across_blocks_ensembles_corr = corr(across_blocks_ensembles);  % ensemble rate correlation between blocks halves of different natural movies 
            
            % cauculate mean ensemble rate across blocks halves of different movie blocks
            mean_across_blocks_ensembles_corr = [];
            for half1 = 1:4 % loop over blocks
                rows = [1:2] +2*(half1-1);
                for half2 = 1:4 % loop over blocks
                    cols = [1:2] +2*(half2-1);
                    current_half = across_blocks_ensembles_corr(rows,cols); % subset rate corr between blocks halves
                    mean_across_blocks_ensembles_corr(half1,half2) = mean(current_half(:),'omitnan'); % average rate corr between blocks halves
                    
                end
            end
            diff_maps(sub,:) = [mean([mean_across_blocks_ensembles_corr(1,2),mean_across_blocks_ensembles_corr(3,4)],'omitnan'),...
                mean([mean_across_blocks_ensembles_corr(1,3),mean_across_blocks_ensembles_corr(2,4)],'omitnan')];  % mean rate corr across blocks of different movies
            same_maps(sub,:) = [mean_across_blocks_ensembles_corr(2,3),mean_across_blocks_ensembles_corr(1,4)];  % mean rate corr across blocks of the same movie
            
            sub = sub + 1;
        end
    end
    diff_maps_areas{area} = diff_maps;
    same_maps_areas{area} = same_maps;
end

same_map_ticks = [60 70]; % interval in minutes between blocks of the same movie
diff_map_ticks = [0 65]; % interval in minutes between blocks of different movies
plt = [];
ylims=[0.69 1;0.7 1;0.65 1;0.69 1;0.7 1;0.775 1];
figure('units','normalized','position',[0.3 0.3 0.4 0.45])  % visualize the ensemble rate corr across block of different movies as function of elapsed time
for area = 1:6
    current_area_diff_maps = diff_maps_areas{area}; % subset the rate corr across blocks of different movies
    current_area_same_maps = same_maps_areas{area}; % subset the rate corr across blocks of the same movies
    
    mean_diff_maps = mean(current_area_diff_maps,'omitnan'); % average across mice
    ste_diff_maps = std(current_area_diff_maps,'omitnan')./sqrt(size(current_area_diff_maps,1)); % standard error across mice
    
    mean_same_maps = mean(current_area_same_maps,'omitnan'); % average across mice
    ste_same_maps = std(current_area_same_maps,'omitnan')./sqrt(size(current_area_same_maps,1)); % standard error across mice
    
    subplot(2,3,area)
    hold on
    plt(1)=errorbar(diff_map_ticks,mean_diff_maps,ste_diff_maps,'o','markerfacecolor',[0.6 0.6 0.6],'linestyle','-',...
        'linewidth',2.5,'color',[0.6 0.6 0.6],'CapSize',0);
    plt(2)=errorbar(same_map_ticks,mean_same_maps,ste_same_maps,'o','markerfacecolor',[0.3 0.3 0.3],'linestyle','-',...
        'linewidth',2.5,'color',[0.3 0.3 0.3],'CapSize',0);
    text(0.2,0.9,brain_areas{area},'Units','normalized','FontSize',15)
    text(0.15,0.8,['N=',num2str(size(current_area_same_maps,1))],'Units','normalized','FontSize',13)
    
    xlim([-9 79])
    ylim(ylims(area,:))
    
    if area == 4 || area == 1
        ylabel('Ensemble rate correlation')
    elseif area == 5
        xlabel('Time between blocks (min)')
    end
    
    if area <4
        set(gca,'xtick',[])
    end
    
    legend(plt,{'Different','Same'},'Location','southwest')
    legend('boxoff')
end

%% Figure S4A - Drifting gratings -  single units examples from area V1

area = 1; % area V1
mouse = 4; % example mouse #4

current_mouse = neuropixels_drifting_gratings{mouse,area}./2; % subset neuronal activity of example mouse during drifting gratings stimuli
% each block was presented for 2 seconds therefore deviding the average
% activity by 2 will result in the average activity in spikes/sec

unit_list = [63,104,98]; % list of example units to be visualized
ylims = [0 20; 0 20; 0 3.25];

plt = [];
figure('units','normalized','position',[0.3 0.3 0.4 0.3]) % visualize neuronal responses to drifting gratings for three exmaple units
for unit = 1:length(unit_list) % loop over units
    current_cell_ori = [];
    for block = 1:3 % loop over blocks
        current_cell_ori(block,:) = mean(current_mouse(unit_list(unit),:,block,:),4,'omitnan'); % calculate the average across temporal frequancies for a single example unit
    end
    
    subplot(1,3,unit) 
    hold on
    plt(1) = plot(current_cell_ori(1,:),'color',[0.3 0.3 0.3],'linewidth',3); % tuning during block A
    plt(2) = plot(current_cell_ori(2,:),'color',[0.55 0.55 0.55],'linewidth',3); % tuning during block B
    plt(3) = plot(current_cell_ori(3,:),'color',[0.8 0.8 0.8],'linewidth',3);  % tuning during block C
    xlim([1 8])
    ylim(ylims(unit,:))
    title(['Unit #',num2str(unit_list(unit))])
    if unit == 1
        ylabel('Mean activity rate (Hz)')
    elseif unit== 2
        xlabel('Drifting gratings direction')
    elseif unit == 3
        legend(plt,{'Block A','Block B','Block C'},'Location','northeast')
        legend('boxoff')
    end
    set(gca,'xtick',1:8,'xticklabel',0:45:315)
    xtickangle(45)

end

%% Figure S4B - Drifting gratings - PV correlation across V1 mice
area = 1; % area V1
cell_cutoff = 20; % treshold for minimum number of cells recorded
valid_mice = calcium_excitatory_cell_count{area}(:,1)>=cell_cutoff; % including only mice with atleast 20 cells recorded during session A
current_area = calcium_excitatory_drifting_gratings{area}(valid_mice); % subset mice that passed the requirments of 'valid_mice'

pv_corr_across_mice = [];
for mouse = 1:length(current_area) % loop over mice
    current_mouse = current_area{mouse}; % subset neuronal activity of a single mouse (#cells by 8 directions by 3 blocks by 5 temporal frequamcies)
    
    % reshape neuronal activity matrix from 4D into 2D
    current_mouse_sorted = [];
    for block =1:3 % loop over blocks
        for ori =1:8 % loop over gratings directions
            for freq = 1:5 % loop over temporal frequancies
                current_mouse_sorted = [current_mouse_sorted,current_mouse(:,ori,block,freq)];
            end
        end
    end
    pv_corr_across_mice(:,:,mouse) = corr(current_mouse_sorted,'rows','pairwise');% calculate the PV correlation between combinations of direction and temporal frequancies across blocks
end


mean_pv_corr_across_mice = mean(pv_corr_across_mice,3,'omitnan'); % average PV corr matrices across mice
main_diag_ind = boolean(eye(size(mean_pv_corr_across_mice,1)));
mean_pv_corr_across_mice(main_diag_ind) = NaN; % for visualiztion - convert main diagonal (values of 1 by definition) into NaN values
mean_pv_corr_across_mice(main_diag_ind) = max(mean_pv_corr_across_mice(:)); % for visualiztion - convert main diagonal into maximal value

figure % main - visualize the PV corr between block of drifting gratings
subplot(1,2,1,'units','normalized','position',[0.1 0.3 0.47 0.6])
imagesc(mean_pv_corr_across_mice)
hold on
for line = 1:24
    plot(xlim,[5.5 5.5]+5*(line-1),'color',newmap3(1,:),'linewidth',0.5)
    plot([5.5 5.5]+5*(line-1),ylim,'color',newmap3(1,:),'linewidth',0.5)
end
plot(xlim,[40.5 40.5],'color',newmap3(1,:),'linewidth',2)
plot(xlim,[40.5 40.5]+40,'color',newmap3(1,:),'linewidth',2)
plot([40.5 40.5],ylim,'color',newmap3(1,:),'linewidth',2)
plot([40.5 40.5]+40,ylim,'color',newmap3(1,:),'linewidth',2)

set(gca,'xtick',20:40:120,'xticklabel',{'Block A','Block B','Block C'},...
    'ytick',20:40:120,'yticklabel',{'Block A','Block B','Block C'})
ytickangle(90)
title('Calcium imaging - area V1')
colormap(newmap3)
cb = colorbar;
set(cb,'position',[0.6 0.3 0.04 0.325])
cb.Label.String = 'PV correlation';
cb.FontSize = 12;

% inset - calculate the average pv corr across corresponding pairs of
% directions X temporal frequancies combinations
pv_ori = [];
for row = 1:16 % loop over directions X temporal frequancies combinations
    rows_ind = [41:45]+5*(row-1);
    pv_ori(:,:,row) = mean_pv_corr_across_mice(rows_ind-40,rows_ind);
    
    if row<9
        pv_ori(:,:,row+16) = mean_pv_corr_across_mice(rows_ind-40,rows_ind+40);
    end
    
end
subplot(1,2,2,'units','normalized','position',[0.6 0.7 0.15 0.2])
imagesc(mean(pv_ori,3,'omitnan'))
set(gca,'xtick',1:5,'xticklabel',[1,2,4,8,15],'ytick',[])
title({'Temporal';'frequancy (Hz)'})
colormap(newmap3)

%% Figure S4C - Drifting gratings - PV correlation across directions

cell_cutoff = 20; % treshold for minimal number of recorded cells
for area = 1:6 % loop over areas
    valid_mice = calcium_excitatory_cell_count{area}(:,1)>=cell_cutoff; % include only mice with atleast 20 recorded cells in session A
    current_area = calcium_excitatory_drifting_gratings{area}(valid_mice); % subset mice that met the requirements of 'valid_mice'
    
    mean_elapsed_pv_shuffled = [];
    mean_elapsed_pv = [];
    for mouse = 1:length(current_area) % loop over mice
        clc;
        disp(['Calculating population vector correlation between blocks:'])
        disp(['Dataset: Calcium imaging | Stimulus: Drifting gratings | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse}; % subset neuronal activity of a single mouse (#cells by 8 directions by 3 blocks by 5 temporal frequamcies)
        
        % reshape neuronal activity from 4D into 3D (#cells by 40 directions X temporal frequancy combinations by 3 blocks) 
        current_mouse_sorted = [];
        for ori =1:size(current_mouse,2) % loop over grating directions
            for freq = 1:size(current_mouse,4) % loop over temporal frequancies
                current_mouse_sorted = [current_mouse_sorted,current_mouse(:,ori,:,freq)];
            end
        end
        
        % cyclic temporal shuffling of neuronal activity as a control
        current_mouse_sorted_shuffled = current_mouse_sorted;
        for cell = 1:size(current_mouse_sorted,1) % loop over cells
            for block = 1:3 % loop over blocks
                current_mouse_sorted_shuffled(cell,:,block) = current_mouse_sorted(cell,circshift([1:40],randperm(40,1)),block); % random cyclic temporal shuffle of single cells tuning curves
            end
        end
        
        
        
        % calculate the mean PV correlation as function of the difference
        % in degrees between gratings directions
        elapsed_pv = [];
        elapsed_pv_shuffled = [];
        ori_list = [5:5:35];
        ori_list2 = [-35:5:-5];
        block = 1;
        for block1 = 1:size(current_mouse_sorted,3) % loop over blocks
            current_mouse_block1 = current_mouse_sorted(:,:,block1);
            current_mouse_block1_shuffled = current_mouse_sorted_shuffled(:,:,block1); % subset neuronal activity during a single block
            
            for block2 = 1:size(current_mouse_sorted,3) % loop over blocks
                current_mouse_block2 = current_mouse_sorted(:,:,block2);
                current_mouse_block2_shuffled = current_mouse_sorted_shuffled(:,:,block2); % subset neuronal activity during a single block
                
                structure = corr(current_mouse_block1,current_mouse_block2);
                structure_shuffled = corr(current_mouse_block1_shuffled,current_mouse_block2_shuffled);
                
                if block1 < block2
                    for diagonal = 1:length(ori_list)+1 % loop over direction differences
                        if diagonal == 1 % for corresponding combinations across blocks
                            elapsed_pv(block,diagonal) = mean(diag(structure),'omitnan'); % calculate the mean pv across combinations of directions and temporal frequancies for non-shuffled data
                            elapsed_pv_shuffled(block,diagonal) = mean(diag(structure_shuffled),'omitnan'); % calculate the mean pv across combinations of directions and temporal frequancies for shuffled data
                            
                        else % across different combinations across blocks
                            elapsed_pv(block,diagonal) = mean([diag(structure,ori_list(diagonal-1));diag(structure,ori_list2(diagonal-1))],'omitnan'); % calculate the mean pv across combinations of directions and temporal frequancies for non-shuffled data
                            elapsed_pv_shuffled(block,diagonal) = mean([diag(structure_shuffled,ori_list(diagonal-1));diag(structure_shuffled,ori_list2(diagonal-1))],'omitnan'); % calculate the mean pv across combinations of directions and temporal frequancies for shuffled data
                            
                        end
                    end
                    block = block + 1;
                end
            end
        end
       
        mean_elapsed_pv(:,:,mouse) = [mean(elapsed_pv([1,3],:,:),'omitnan');elapsed_pv(2,:,:)]; % store mean pv corr values as function of elapsed timefor non-shuffled data
        mean_elapsed_pv_shuffled(:,:,mouse) = [mean(elapsed_pv_shuffled([1,3],:,:),'omitnan');elapsed_pv_shuffled(2,:,:)];  % store mean pv corr values as function of elapsed time for shuffled data
    end
    
    mean_elapsed_pv_area{area} = mean_elapsed_pv; 
    mean_elapsed_pv_area_shuffled{area} = mean_elapsed_pv_shuffled;
end


figure('units','normalized','position',[0.3 0.3 0.3 0.3]) % visualize pv corr as function of direction difference for each visual area
for area = 1:6 % loop over areas
    subplot(2,3,area)
    
    current_area = mean([mean_elapsed_pv_area{area}(:,5:8,:),mean_elapsed_pv_area{area}(:,1:5,:)],'omitnan'); % average pv corr values across time and reshape according to direction difference for non-shuffled data
    mean_current_area_shuffled = mean([mean_elapsed_pv_area_shuffled{area}(:,5:8,:),mean_elapsed_pv_area_shuffled{area}(:,1:5,:)],'omitnan'); % average pv corr values across time and reshape according to direction difference for shuffled data
   
    mean_pv = mean(current_area,3,'omitnan'); % calculate mean across mice for non-shuffled data
    ste_pv = std(current_area,[],3,'omitnan')./sqrt(size(current_area,3)); % calculate standard error across mice for non-shuffled data
    
    mean_pv_shuffled = mean(mean_current_area_shuffled,3,'omitnan'); % calculate mean across mice for non-shuffled data
    ste_pv_shuffled = std(mean_current_area_shuffled,[],3,'omitnan')./sqrt(size(mean_current_area_shuffled,3)); % calculate standard error across mice for non-shuffled data
    
    hold on
    errorbar(mean_pv_shuffled,ste_pv_shuffled,'o','color',[0.7 0.7 0.7],...
        'markerfacecolor',[0.7 0.7 0.7],'capsize',0,'linestyle','-','linewidth',2,'markersize',2)
    
    errorbar(mean_pv,ste_pv,'o','color',colors(area,:),...
        'markerfacecolor',colors(area,:),'capsize',0,'linestyle','-','linewidth',2,'markersize',2)
    
    text(0.1,0.95,brain_areas{area},'Units','normalized','FontSize',12)
    text(0.05,0.85,['N=',num2str(sum(size(current_area,3),2))],'Units','normalized','FontSize',10)
    
    text(0.8,0.95,['Data'],'Units','normalized','FontSize',10,'Color',colors(area,:))
    text(0.75,0.85,['Shuffle'],'Units','normalized','FontSize',10,'Color',[0.7 0.7 0.7])
    
    xlim([0 10])
    ylim([0 0.7])
    
    set(gca,'xtick',1:2:9,'xticklabel',-180:90:180)
    
    if area == 4 || area == 1
        ylabel('PV correlation')
    elseif area == 5
        xlabel('Direction difference (degrees)')
    end
    
end

%% Figure S4D - Drifting gratings - Ensemble rate correlation across blocks (calcium imaging)

cell_cutoff = 20; % treshold for minimal number of recorded cells
elapsed_block_rate_corr_area = {};
for area = 1:6 % loop over areas
    valid_mice = calcium_excitatory_cell_count{area}(:,1)>=cell_cutoff; % include only mice with atleast 20 recorded cells in session A
    current_area = calcium_excitatory_drifting_gratings{area}(valid_mice); % subset mice that met the requirements of 'valid_mice'
    
    elapsed_block_rate_corr = [];
    for mouse = 1:length(current_area) % loop over mice
        clc;
        disp(['Calculating ensemble rate correlation between blocks:'])
        disp(['Dataset: Calcium imaging | Stimulus: Drifting gratings | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse};  % subset neuronal activity of a single mouse (#cells by 8 directions by 3 blocks by 5 temporal frequamcies)
         
        current_mouse_sorted = []; % reshape neuronal activity from 4D into 3D (#cells by 40 directions X temporal frequancy combinations by 3 blocks) 
        for freq = 1:size(current_mouse,4) 
            current_mouse_sorted = [current_mouse_sorted,current_mouse(:,:,:,freq)];
        end
        
        mean_activity_each_block = squeeze(mean(current_mouse_sorted,2,'omitnan')); % average the activity across direction X temporal freq combinations
        rate_corr = corr(mean_activity_each_block); % calculate ensemble rate corr btween blocks
        
        elapsed_block_rate_corr(mouse,:) = [mean(diag(rate_corr,1),'omitnan'),diag(rate_corr,2)]; % calculate ensemble rate as function of elapsed time
        
    end
    elapsed_block_rate_corr_area{area} = elapsed_block_rate_corr; 
end

pvalues = [];
zvalues = [];
plt = [];
figure('units','normalized','position',[0.35 0.35 0.2 0.3]) % visualize the ensemble rate corr as function of elapsed time
for area = 1:6 % loop over areas
    current_area = elapsed_block_rate_corr_area{area}; % subset rate corr values for a single area
    
    [pvalues(area),~,stats] = signrank(current_area(:,1),current_area(:,2)); % two-sided wilcoxon signed-rank test for difference between proximal and distal blocks
    zvalues(area) = stats.zval;
    
    mean_stability = mean(current_area,'omitnan'); % average across mice
    ste_stability = std(current_area,'omitnan')./sqrt(size(current_area,1)); % standard error across mice
    hold on
    plt(area) = errorbar(mean_stability,ste_stability,'o','color',colors(area,:),...
        'markerfacecolor',colors(area,:),'capsize',0,'linestyle','-','linewidth',3,'markersize',3);
end
legend(plt,brain_areas(1:6),'Location','southwest')
legend('boxoff')
ylabel('Ensemble rate correlation')
set(gca,'xtick',1:2,'xticklabels',[15,30])
xlabel('Elapsed time (min)')
xlim([0.5 2.5])

% correct for multiple comparisons using bonferroni-holm method
corrected_pval = bonf_holm(pvalues);

% define statistics summary table
VarNames = {'area','zvalue','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas(1:6)'),zvalues(:),pvalues(:),corrected_pval(:),'VariableNames',VarNames);

% display statistics summary table
clc;
disp(['Ensemble rate correlation between proximal blocks compared to between distal blocks'])
disp(['Two-tailed Wilcoxon signed-rank tests with Holm�Bonferroni correction:'])
disp(statistics)

%% Figure S4E - Drifting gratings - Tuning curve correlation across blocks (calcium imaging)
cell_cutoff = 20; % treshold for minimal number of recorded cells
elapsed_block_tuning_corr_area = {};
for area = 1:6 % loop over areas
    valid_mice = calcium_excitatory_cell_count{area}(:,1)>=cell_cutoff; % include only mice with atleast 20 recorded cells in session A
    current_area = calcium_excitatory_drifting_gratings{area}(valid_mice); % subset mice that met the requirements of 'valid_mice'
    
    elapsed_block_tuning_corr = [];
    for mouse = 1:length(current_area) % loop over mice
        clc;
        disp(['Calculating tuning curve correlation between blocks:'])
        disp(['Dataset: Calcium imaging | Stimulus: Drifting gratings | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse};  % subset neuronal activity of a single mouse (#cells by 8 directions by 3 blocks by 5 temporal frequamcies)
        
        current_mouse_sorted = []; % reshape neuronal activity from 4D into 3D (#cells by 40 directions X temporal frequancy combinations by 3 blocks) 
        for freq = 1:size(current_mouse,4)
            current_mouse_sorted = [current_mouse_sorted,current_mouse(:,:,:,freq)];
        end
       
        mean_tuning_corr = [];
        for block1 = 1:size(current_mouse_sorted,3) % loop over blocks
            current_mouse_block1 = current_mouse_sorted(:,:,block1); % subset neuronal activity during a single block
            for block2 = 1:size(current_mouse_sorted,3) % loop over blocks
                current_mouse_block2 = current_mouse_sorted(:,:,block2); % subset neuronal activity during a single block
                mean_tuning_corr(block1,block2) = median(diag(corr(current_mouse_block1',current_mouse_block2','rows','complete')),'omitnan'); % calculate the median tuning correlation of corresponding cells between block
            end
        end
        
        elapsed_block_tuning_corr(mouse,:) = [mean(diag(mean_tuning_corr,1),'omitnan'),diag(mean_tuning_corr,2)]; % calculate tuning curve as function of elapsed time
    end
    elapsed_block_tuning_corr_area{area} = elapsed_block_tuning_corr;
end

pvalues = [];
zvalues = [];
figure('units','normalized','position',[0.35 0.35 0.2 0.3]) % visualize the tuning curve corr as function of elapsed time
for area = 1:6
    current_area = elapsed_block_tuning_corr_area{area}; % subset tuning curve corr values for a single area
    
    [pvalues(area),~,stats] = signrank(current_area(:,1),current_area(:,2)); % two-sided wilcoxon signed-rank test for difference between proximal and distal blocks
    zvalues(area) = stats.zval;
    
    mean_stability = mean(current_area,'omitnan'); % average across mice
    ste_stability = std(current_area,'omitnan')./sqrt(size(current_area,1)); % standard error across mice
    hold on
    errorbar(mean_stability,ste_stability,'o','color',colors(area,:),...
        'markerfacecolor',colors(area,:),'capsize',0,'linestyle','-','linewidth',3,'markersize',3);
    
end

ylabel('Tuning curve correlation')
set(gca,'xtick',1:2,'xticklabels',[15,30])
xlabel('Elapsed time (min)')
xlim([0.5 2.5])

% correct for multiple comparisons using bonferroni-holm method
corrected_pval = bonf_holm(pvalues);

% define statistics summary table
VarNames = {'area','zvalue','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas(1:6)'),zvalues(:),pvalues(:),corrected_pval(:),'VariableNames',VarNames);

% display statistics summary table
clc;
disp(['Tuning curve correlation between proximal blocks compared to between distal blocks'])
disp(['Two-tailed Wilcoxon signed-rank tests with Holm�Bonferroni correction:'])
disp(statistics)


%% Figure S4F - Drifting gratings - Ensemble rate correlation across blocks (neuropixels)
cell_cutoff = 15; % treshold for minimum number of cells recorded
elapsed_block_rate_corr_area = {};
for area = 1:6 % loop over areas
    valid_mice = neuropixels_cell_count(:,area,1) >= cell_cutoff & movie_repeats(:,1) == 10; % include only mice from the 'Brain Observatory' group with atleast 15 recorded units
    current_area = neuropixels_drifting_gratings(valid_mice,area); % subset mice that met the requirements of 'valid_mice'
    
    elapsed_block_rate_corr = [];
    for mouse = 1:length(current_area) % loop over mice
        clc;
        disp(['Calculating ensemble rate correlation between blocks:'])
        disp(['Dataset: Neuropixels | Stimulus: Drifting gratings | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse}; % subset neuronal activity of a single mouse (#units by 8 directions by 3 blocks by 5 temporal frequamcies)
        
        current_mouse_sorted = []; % reshape neuronal activity from 4D into 3D (#cells by 40 directions X temporal frequancy combinations by 3 blocks) 
        for freq = 1:5 % loop over temporal frequancies
            current_mouse_sorted = [current_mouse_sorted,current_mouse(:,:,:,freq)];
        end
        
        mean_activity_each_block = squeeze(mean(current_mouse_sorted,2,'omitnan')); % average the activity across direction X temporal freq combinations
        rate_corr = corr(mean_activity_each_block); % calculate ensemble rate corr btween blocks
        
        elapsed_block_rate_corr(mouse,:) = [mean(diag(rate_corr,1),'omitnan'),diag(rate_corr,2)]; % calculate ensemble rate as function of elapsed time
        
    end
    elapsed_block_rate_corr_area{area} = elapsed_block_rate_corr;
end


pvalues = [];
zvalues = [];
plt = [];
figure('units','normalized','position',[0.35 0.35 0.2 0.3])  % visualize the ensemble rate corr as function of elapsed time
for area = 1:6 % loop over areas
    current_area = elapsed_block_rate_corr_area{area}; % subset rate corr values for a single area
    
    [pvalues(area),~,stats] = signrank(current_area(:,1),current_area(:,2)); % two-sided wilcoxon signed-rank test for difference between proximal and distal blocks
    zvalues(area) = stats.zval;
    
    mean_stability = mean(current_area,'omitnan'); % average across mice
    ste_stability = std(current_area,'omitnan')./sqrt(size(current_area,1)); % standard error across mice
    
    hold on
    plt(area) = errorbar(mean_stability,ste_stability,'o','color',colors(area,:),...
        'markerfacecolor',colors(area,:),'capsize',0,'linestyle','-','linewidth',3,'markersize',3);
    
end
legend(plt,brain_areas(1:6),'Location','southwest')
legend('boxoff')
ylabel('Ensemble rate correlation')
set(gca,'xtick',1:2,'xticklabels',[15,30])
xlabel('Elapsed time (min)')
xlim([0.5 2.5])

% correct for multiple comparisons using bonferroni-holm method
corrected_pval = bonf_holm(pvalues);

% define statistics summary table
VarNames = {'area','zvalue','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas(1:6)'),zvalues(:),pvalues(:),corrected_pval(:),'VariableNames',VarNames);

% display statistics summary table
clc;
disp(['Ensemble rate correlation between proximal blocks compared to between distal blocks'])
disp(['Two-tailed Wilcoxon signed-rank tests with Holm�Bonferroni correction:'])
disp(statistics)

%% Figure S4G - Drifting gratings - Tuning curve correlation across blocks (neuropixels)
cell_cutoff = 15; % treshold for minimum number of cells recorded
elapsed_block_tuning_corr_area = {};
for area = 1:6 % loop over areas
    valid_mice = neuropixels_cell_count(:,area,1) >= cell_cutoff & movie_repeats(:,1) == 10; % include only mice from the 'Brain Observatory' group with atleast 15 recorded units
    current_area = neuropixels_drifting_gratings(valid_mice,area); % subset mice that met the requirements of 'valid_mice'
    
    elapsed_block_tuning_corr = [];
    for mouse = 1:length(current_area) % loop over mice
        clc;
        disp(['Calculating tuning curve correlation between blocks:'])
        disp(['Dataset: Neuropixels | Stimulus: Drifting gratings | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse}; % subset neuronal activity of a single mouse (#units by 8 directions by 3 blocks by 5 temporal frequamcies)
        
        current_mouse_sorted = []; % reshape neuronal activity from 4D into 3D (#cells by 40 directions X temporal frequancy combinations by 3 blocks) 
        for freq = 1:5 % loop over temporal frequancies
            current_mouse_sorted = [current_mouse_sorted,current_mouse(:,:,:,freq)];
        end

        mean_tuning_corr = [];
        for block1 = 1:size(current_mouse_sorted,3) % loop over blocks
            current_mouse_block1 = current_mouse_sorted(:,:,block1); % subset neuronal activity during a single block
            for block2 = 1:size(current_mouse_sorted,3) % loop over blocks
                current_mouse_block2 = current_mouse_sorted(:,:,block2); % subset neuronal activity during a single block
                mean_tuning_corr(block1,block2) = median(diag(corr(current_mouse_block1',current_mouse_block2','rows','complete')),'omitnan'); % calculate the median tuning correlation of corresponding cells between block
            end
        end
        elapsed_block_tuning_corr(mouse,:) = [mean(diag(mean_tuning_corr,1),'omitnan'),diag(mean_tuning_corr,2)]; % calculate tuning curve as function of elapsed time
    end
    elapsed_block_tuning_corr_area{area} = elapsed_block_tuning_corr;
end

pvalues = [];
zvalues = [];
figure('units','normalized','position',[0.35 0.35 0.2 0.3]) % visualize the tuning curve corr as function of elapsed time
for area = 1:6
    current_area = elapsed_block_tuning_corr_area{area}; % subset tuning curve corr values for a single area
    
    [pvalues(area),~,stats] = signrank(current_area(:,1),current_area(:,2)); % two-sided wilcoxon signed-rank test for difference between proximal and distal blocks
    zvalues(area) = stats.zval;
    
     mean_stability = mean(current_area,'omitnan'); % average across mice
    ste_stability = std(current_area,'omitnan')./sqrt(size(current_area,1)); % standard error across mice
    hold on
    errorbar(mean_stability,ste_stability,'o','color',colors(area,:),...
        'markerfacecolor',colors(area,:),'capsize',0,'linestyle','-','linewidth',3,'markersize',3);
    
end

ylabel('Tuning curve correlation')
set(gca,'xtick',1:2,'xticklabels',[15,30])
xlabel('Elapsed time (min)')
xlim([0.5 2.5])

% correct for multiple comparisons using bonferroni-holm method
corrected_pval = bonf_holm(pvalues);

% define statistics summary table
VarNames = {'area','zvalue','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas(1:6)'),zvalues(:),pvalues(:),corrected_pval(:),'VariableNames',VarNames);

% display statistics summary table
clc;
disp(['Tuning curve correlation between proximal blocks compared to between distal blocks'])
disp(['Two-tailed Wilcoxon signed-rank tests with Holm�Bonferroni correction:'])
disp(statistics)

%% Figure S5A - Ensemble rate and tuning curve correlation difference between proximal and distal sessions

nat_movie = 1; % natural movie 1
cell_cutoff = 20; % cell count threshold of 20 cells

rate_tuning_corr_diff_across_areas = {};
for area = 1:6 % loop over visual areas
    % valid_mice - true only if the number of cells recorded is
    % atleast 20 in each of the recording sessions
    valid_mice = min(calcium_excitatory_cell_count{area},[],2)>= cell_cutoff;
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie); % mice that passed the requirments of 'valid_mice'
    
    elapsed_rate_corr = []; % define empty variable that will store the rate correlation values for all mice
    elapsed_tuning_corr = []; % define empty variable that will store the tuning correlation values for all mice
    for mouse = 1:length(current_area) % loop over mice
        clc;
        disp(['Calculating ensemble rate and tuning correlation between ssessions:'])
        disp(['Dataset: Calcium imaging | Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse}; % subset neuronal activity for a single mouse (#cells by 30 time bins by 30 movie repeats)
        
       
        %calculate the average activity rate across movie repeats for each session
        mean_activity_per_sess = []; % define an empty variable that will store the averaged neuronal responses of each session
        for sess = 1:3 % loop over session halves
            current_sess = [1:10] + 10*(sess-1); % define range of movie repeats for current session
            mean_activity_per_sess(:,:,sess) = mean(current_mouse(:,:,current_sess),3,'omitnan'); % average the neuronal activity across chosen movie repeats
        end
        
        % calculate rate and tuning correlation between sessions using only
        % cells that were active in both compared time points
        rate_corr = []; % define empty variable that will store rate corr values for all mice
        tuning_corr = []; % define empty variable that will store tuning corr values for all mice
        for sessA = 1:size(mean_activity_per_sess,3) % loop over sessions
            sessA_activity = mean_activity_per_sess(:,:,sessA); % subset neuronal activity for a single session
            
            for sessB = 1:size(mean_activity_per_sess,3) % loop over sessions
                sessB_activity = mean_activity_per_sess(:,:,sessB); % subset neuronal activity for a single session
                
                % valid_cells - cells that were active in both of compared time points
                valid_cells = [mean(sessA_activity,2,'omitnan') ~= 0] & [mean(sessB_activity,2,'omitnan') ~= 0]; % true if average activity is above zero in both sessions
                valid_sessA_activity = sessA_activity(valid_cells,:); % subset only the cells that met the requirments of 'valid_cells'
                valid_sessB_activity = sessB_activity(valid_cells,:); % subset only the cells that met the requirments of 'valid_cells'

                rate_corr(sessA,sessB) = corr(mean(valid_sessA_activity,2,'omitnan'),mean(valid_sessB_activity,2,'omitnan')); % calculate ensemble rate corr between sessions
                tuning_corr(sessA,sessB) = median(diag(corr(valid_sessA_activity',valid_sessB_activity')),'omitnan'); % calculate the median tuning curve correlation of corresponding cells between sessions
            end
        end
         proximal_sessions_rate = mean(diag(rate_corr,1),'omitnan'); % rate corr values between proximal sessions (sessions 1&2 and sessions 2&3)
        distal_sessions_rate = diag(rate_corr,2); % rate corr values between distal sessions (sessions 1&3)
        elapsed_rate_corr(mouse,:) = [proximal_sessions_rate,distal_sessions_rate];
        
          proximal_sessions_tuning = mean(diag(tuning_corr,1),'omitnan'); % rate corr values between proximal sessions (sessions 1&2 and sessions 2&3)
        distal_sessions_tuning = diag(tuning_corr,2); % rate corr values between distal sessions (sessions 1&3)
        elapsed_tuning_corr(mouse,:) = [proximal_sessions_tuning,distal_sessions_tuning];
    end
    
    rate_corr_diff = elapsed_rate_corr(:,1) - elapsed_rate_corr(:,2); % calculate the rate corr difference between proximal and distal sessions
    tuning_corr_diff = elapsed_tuning_corr(:,1) - elapsed_tuning_corr(:,2); % calculate the tuning corr difference between proximal and distal sessions
    rate_tuning_corr_diff_across_areas{area} = [rate_corr_diff,tuning_corr_diff];
    
end

pvalue = [];
zvalue = [];
mean_stability = [];
ste_stability = [];
for area = 1:6 % loop over areas
    current_area = rate_tuning_corr_diff_across_areas{area}; % subset corr differnce values for a single area
    mean_stability(area,:) = mean(current_area,'omitnan');  % calculate the average across mice
    std_stability = std(current_area,'omitnan'); % calculate the standard deviation across mice
    ste_stability(area,:) = std_stability./sqrt(size(current_area,1)); % calculate standard error across mice
    
    [pvalue(area,1),~,stats] = signrank(current_area(:,1)); % perform one-sided wilcoxon signed-rank test for difference in rate corr between proximal and distal sessions
    zvalue(area,1) = stats.zval;
    
    [pvalue(area,2),~,stats] = signrank(current_area(:,2)); % perform one-sided wilcoxon signed-rank test for difference in tuning corr between proximal and distal sessions
    zvalue(area,2) = stats.zval;
end

plt = [];
figure('Units','Normalized','Position',[0.3 0.4 0.225 0.325]) % visualization of rate and tuning corr between sessions across areas and mice
xlim([0 7])
plot(xlim,[0 0],'--','color',[0.2 0.2 0.2],'linewidth',1.5)
hold on
plt(1) = errorbar(mean_stability(:,1),ste_stability(:,1),'o','color',[0.7 0.7 0.7],...
    'markerfacecolor',[0.7 0.7 0.7],'capsize',0,'linestyle','none','linewidth',3,'markersize',5);
plt(2) = errorbar(mean_stability(:,2),ste_stability(:,2),'o','color',[0.5 0.5 0.5],...
    'markerfacecolor',[0.5 0.5 0.5],'capsize',0,'linestyle','none','linewidth',3,'markersize',5);

lgd = legend(plt,{'Ensemble rate correlation','Tuning curve correlation'},'Location','southwest');
legend('boxoff')
ylim([-0.04 0.11])
set(gca,'xtick',1:6,'xticklabel',brain_areas(1:6),'box','off')
ylabel({'Correlation difference';'(within block - between blocks)'})

% correct for multiple comparisons using bonferroni-holm method
pvalue = pvalue./2; % one tail
corrected_pvalue = [];
corrected_pvalue(:,1) = bonf_holm(pvalue(:,1)); % correct for rate corr
corrected_pvalue(:,2) = bonf_holm(pvalue(:,2)); % correct for tuning corr

% define statistics summary table
VarNames = {'area','Rate_zvalue','Rate_pvalue','Rate_bonf_holm','Tuning_zvalue','Tuning_pvalue','Tuning_bonf_holm'};
statistics = table(cell2mat(brain_areas(1:6)'),zvalue(:,1),pvalue(:,1),corrected_pvalue(:,1),...
    zvalue(:,2),pvalue(:,2),corrected_pvalue(:,2),'VariableNames',VarNames);

% display statistics summary table
clc;
disp(['Ensemble rate and tuning curve correlation within blocks compared to between blocks'])
disp(['One-tailed Wilcoxon signed-rank tests with Holm�Bonferroni correction:'])
disp(statistics)

%% Figure S5B - PV correlation across sessions with cells turnover
nat_movie = 1; % natural movie 1
cell_cutoff = 20; % cell count threshold of 20 cells

pv_corr_areas = {}; % define empty variable that will store the pv corr values of within and between sessions across mice and visual areas
for area = 1:6 % loop over visual areas
    % valid_mice - true only if the number of cells recorded is
    % atleast 20 in each of the recording sessions
    valid_mice = min(calcium_excitatory_cell_count{area},[],2)>= cell_cutoff;
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie);  % mice that passed the requirments of 'valid_mice'
    
    for subset = 1:2 % loop over subset of neuronal data
        within_between_session_stability = []; % define empty variable that will store the pv correlation values for all mice
        for mouse = 1:length(current_area) % loop over mice
            clc;
            disp(['Calculating population vector correlation between ssessions:'])
            disp(['Dataset: Calcium imaging | Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
            
            current_mouse = current_area{mouse}; % subset neuronal activity for a single mouse (#cells by 30 time bins by 30 movie repeats)
            
            %calculate the average activity rate for session halves across movie repeats
            mean_activity_per_half = []; % define an empty variable that will store the averaged neuronal responses of each session half
            for half = 1:6 % loop over session halves
                current_half = [1:5] + 5*(half-1); % define range of movie repeats for current half
                mean_activity_per_half(:,:,half) = mean(current_mouse(:,:,current_half),3,'omitnan');  % average the neuronal activity across chosen movie repeats
            end

            mean_pv_corr = []; % define empty variable that will store pv corr values for all mice
            for halfA = 1:size(mean_activity_per_half,3) % loop over session halves
                halfA_activity = mean_activity_per_half(:,:,halfA); % subset neuronal activity of a single session half
                for halfB = 1:size(mean_activity_per_half,3) % loop over session halves
                    halfB_activity = mean_activity_per_half(:,:,halfB); % subset neuronal activity of a single session half
                    
                    if subset == 1 % cells that were active in both of the compared time points
                        valid_cells = [mean(halfA_activity,2,'omitnan') ~= 0] & [mean(halfB_activity,2,'omitnan') ~= 0];
                    elseif subset == 2 % cells that were active in at least one of the compared time points
                        valid_cells = [mean(halfA_activity,2,'omitnan') ~= 0] | [mean(halfB_activity,2,'omitnan') ~= 0];
                    end
                    
                    valid_halfA_activity = halfA_activity(valid_cells,:); % subset only the cells that met the requirments of 'valid_cells'
                    valid_halfB_activity = halfB_activity(valid_cells,:); % subset only the cells that met the requirments of 'valid_cells'
                    
                    pv_corr = corr(valid_halfA_activity,valid_halfB_activity); % calculate the pv correlation between time bins of different halves
                    mean_pv_corr(halfA,halfB) = mean(diag(pv_corr),'omitnan'); % calculate the mean pv corr across corresponding time bins
                    
                end
            end
            
             within_session_values = [mean_pv_corr(1,2),mean_pv_corr(3,4),mean_pv_corr(5,6)]; % pv corr values between halves of the same session (within session)
            proximal_sessions_values = [mean_pv_corr(1:2,3:4),mean_pv_corr(3:4,5:6)]; % pv corr values between halves of proximal sessions (sessions 1&2 and sessions 2&3)
            distal_session_values = mean_pv_corr(1:2,5:6); % pv corr values between halves of distal sessions (sessions 1&3)
            
            within_between_session_stability(mouse,1) = mean(within_session_values(:),'omitnan'); % average pv corr values for within session
            within_between_session_stability(mouse,2) = mean(proximal_sessions_values(:),'omitnan'); % average pv corr values for proximal sessions
            within_between_session_stability(mouse,3) = mean(distal_session_values(:),'omitnan'); % average pv corr values for distal sessions
        
        end
        
        pv_corr_areas{area,subset} =  within_between_session_stability; % store mean pv corr values for all mice of a given area 
    end
    
end

% define empty variables that will store the results of Wilcoxon signed-rank tests
pvalue = []; % define empty variable for the pvalues
zvalue = []; % define empty variable for the z values

figure('units','normalized','position',[0.3 0.3 0.4 0.4])
for area = 1:6 % loop over visual areas
    for subset = 1:2 % loop over dataset subsets
        current_area = pv_corr_areas{area,subset}; % subset pv corr values for a single area
        
        [pvalues(subset,area),~,stats] = signrank(current_area(:,2),current_area(:,3),'tail','right'); % perform one-sided wilcoxon signed-rank test for difference between proximal and distal sessions
        zvalues(subset,area) = stats.zval; % store z values
        
        mean_stability = mean(current_area,'omitnan'); % calculate mean pv across mice
        std_stability = std(current_area,'omitnan'); % calculate standard deviation across mice
        ste_stability = std_stability./sqrt(size(current_area,1)); % calculate standard error across mice
        
        subplot(2,3,area)
        hold on
        if subset == 1
            plt(subset) = errorbar(mean_stability,ste_stability,'o','color',[0.7 0.7 0.7],...
                'markerfacecolor',[0.7 0.7 0.7],'capsize',0,'linestyle','-','linewidth',3,'markersize',5);
        else
            plt(subset) = errorbar(mean_stability,ste_stability,'o','color',colors(area,:),...
                'markerfacecolor',colors(area,:),'capsize',0,'linestyle','-','linewidth',3,'markersize',5);
        end
        xlim([0.5 3.5])
        if area > 3
            set(gca,'xtick',1:3,'xticklabel',{'Within session','Proximal sessions','Distal sessions'})
            xtickangle(15)
        else
            set(gca,'xtick',[],'xticklabel',[])
        end
        
        if area == 4 || area == 1
            ylabel('PV correlation')
        end
        
    end
    legend(plt,{'Both','Active>1'},'Location','best')
    legend('boxoff')
end


% correct for multiple comparisons using bonferroni-holm method
corrected_pval = bonf_holm(pvalues(2,:));

% define statistics summary table
VarNames = {'area','zvalue','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas(1:6)'),zvalues(2,:)',pvalues(2,:)',corrected_pval(:),'VariableNames',VarNames);

% display statistics summary table
clc;
disp(['PV correlation  of proximal sessions compared to distal sessions when including only cells active in atleast one of the compred time points'])
disp(['One-tailed Wilcoxon signed-rank tests with Holm�Bonferroni correction:'])
disp(statistics)


%% Figure S5C - Ensemble rate correlation across sessions for spontaneous activity
nat_movie = 1; % natural movie 1
cell_cutoff = 20; % cell count threshold of 20 cells

rate_corr_areas = {}; % define empty variable that will store the rate corr values of within and between sessions across mice and visual areas
for area = 1:6 % loop over visual areas
    % valid_mice - true only if the number of cells recorded is
    % atleast 20 in each of the recording sessions
    valid_mice = min(calcium_excitatory_cell_count{area},[],2)>= cell_cutoff;
    
    current_area_movie = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie); % mice that passed the requirments of 'valid_mice' during natural movies
    current_area_spont = calcium_excitatory_spont_population_vectors{area}(valid_mice,nat_movie); % mice that passed the requirments of 'valid_mice' during spont. activity
    
    elapsed_sess_rate_corr_spont = []; % define empty variable that will store the ensemble rate correlation values during natural movies
    elapsed_sess_rate_corr_movie = []; % define empty variable that will store the ensemble rate correlation values during spont. activity
    for mouse = 1:length(current_area_movie) % loop over mice
        clc;
        disp(['Calculating ensemble rate correlation between ssessions:'])
        disp(['Dataset: Calcium imaging | Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area_movie))])
        
        current_mouse_movie = current_area_movie{mouse}; % subset neuronal activity for a single mouse during natural movies (#cells by 30 time bins by 30 movie repeats)
        current_mouse_spont = current_area_spont{mouse}; % subset neuronal activity for a single mouse during spont. activity (#cells by 30 time bins by 30 movie repeats)
        

        rate_corr_movie = [];
        rate_corr_spont= [];
        for half1 = 1:6 % loop over session halves
            current_half1_movie = mean(current_mouse_movie(:,:,[1:5]+5*(half1-1)),3,'omitnan'); % calculate the average neuronal activity across movie repeats of a single session half during natural movies
            current_half1_spont = mean(current_mouse_spont(:,:,[1:5]+5*(half1-1)),3,'omitnan'); % calculate the average neuronal activity across movie repeats of a single session half during spont. activity
            for  half2 = 1:6 % loop over session halves
                current_half2_movie = mean(current_mouse_movie(:,:,[1:5]+5*(half2-1)),3,'omitnan'); % calculate the average neuronal activity across movie repeats of a single session half during natural movies
                current_half2_spont = mean(current_mouse_spont(:,:,[1:5]+5*(half2-1)),3,'omitnan'); % calculate the average neuronal activity across movie repeats of a single session half during spont. activity
                
                % 'valid_cells' - cells that were active in both of the
                % compared time points during movie presentation
                valid_cells = [mean(current_half1_movie,2,'omitnan') > 0] & [mean(current_half2_movie,2,'omitnan') > 0];
                valid_half1_activity_movie = mean(current_half1_movie(valid_cells,:),2,'omitnan'); % calculate the average activity across time bins only the cells that met the requirments of 'valid_cells'
                valid_half2_activity_movie = mean(current_half2_movie(valid_cells,:),2,'omitnan'); % calculate the average activity across time bins only the cells that met the requirments of 'valid_cells'
                
                % 'valid_cells' - cells that were active in both of the
                % compared time points during spont. activity
                valid_cells = [mean(current_half1_spont,2,'omitnan') > 0] & [mean(current_half2_spont,2,'omitnan') > 0];
                valid_half1_activity_spont = mean(current_half1_spont(valid_cells,:),2,'omitnan'); % calculate the average activity across time bins only the cells that met the requirments of 'valid_cells'
                valid_half2_activity_spont = mean(current_half2_spont(valid_cells,:),2,'omitnan'); % calculate the average activity across time bins only the cells that met the requirments of 'valid_cells'
                
                 rate_corr_movie(half1,half2) = corr(valid_half1_activity_movie,valid_half2_activity_movie); % calculate the ensemble rate corr between sessions halves of natural movies
                rate_corr_spont(half1,half2) = corr(valid_half1_activity_spont,valid_half2_activity_spont); % calculate the ensemble rate corr between sessions halves of spont. activity
            
            end
        end

        % during presentation of natural movies
           within_session_values_movie = [rate_corr_movie(1,2),rate_corr_movie(3,4),rate_corr_movie(5,6)]; % rate corr values between halves of the same session (within session)
            proximal_sessions_values_movie = [rate_corr_movie(1:2,3:4),rate_corr_movie(3:4,5:6)]; % rate corr values between halves of proximal sessions (sessions 1&2 and sessions 2&3)
            distal_session_values_movie = rate_corr_movie(1:2,5:6); % rate corr values between halves of distal sessions (sessions 1&3)
            
        elapsed_sess_rate_corr_movie(mouse,1) = mean(within_session_values_movie(:),'omitnan'); % average rate corr values for within session
        elapsed_sess_rate_corr_movie(mouse,2) = mean(proximal_sessions_values_movie(:),'omitnan'); % average rate corr values for proximal sessions
        elapsed_sess_rate_corr_movie(mouse,3) = mean(distal_session_values_movie(:),'omitnan'); % average rate corr values for distal sessions
        
        
               % during blocks of spont. activity
           within_session_values_spont = [rate_corr_spont(1,2),rate_corr_spont(3,4),rate_corr_spont(5,6)]; % rate corr values between halves of the same session (within session)
            proximal_sessions_values_spont = [rate_corr_spont(1:2,3:4),rate_corr_spont(3:4,5:6)]; % rate corr values between halves of proximal sessions (sessions 1&2 and sessions 2&3)
            distal_session_values_spont = rate_corr_spont(1:2,5:6); % rate corr values between halves of distal sessions (sessions 1&3)
            
        elapsed_sess_rate_corr_spont(mouse,1) = mean(within_session_values_spont(:),'omitnan'); % average rate corr values for within session
        elapsed_sess_rate_corr_spont(mouse,2) = mean(proximal_sessions_values_spont(:),'omitnan'); % average rate corr values for proximal sessions
        elapsed_sess_rate_corr_spont(mouse,3) = mean(distal_session_values_spont(:),'omitnan'); % average rate corr values for distal sessions
        
    end
    
    % store ensemble rate corr values for all mice of a given area 
    rate_corr_areas{1,area} = elapsed_sess_rate_corr_movie;
    rate_corr_areas{2,area} = elapsed_sess_rate_corr_spont;
    
end

% define empty variables that will store the results of Wilcoxon signed-rank tests
pvalues = []; % define empty variable for the pvalues
zvalues = []; % define empty variable for the z values
plt = [];

figure('units','normalized','position',[0.3 0.3 0.35 0.35])
for area = 1:6 % loop over visual areas
    subplot(2,3,area)
    % during presentation of natural movies
    mean_rate_corr_movie = mean(rate_corr_areas{1,area},'omitnan'); % calculate mean across mice
    ste_rate_corr_movie = std(rate_corr_areas{1,area},'omitnan')./sqrt(size(rate_corr_areas{1,area},1)); % calculate standard error across mice
    
    % during blocks of spont. activity
    mean_rate_corr_spont = mean(rate_corr_areas{2,area},'omitnan'); % calculate mean across mice
    ste_rate_corr_spont = std(rate_corr_areas{2,area},'omitnan')./sqrt(size(rate_corr_areas{2,area},1)); % calculate standard error across mice
    
    
    [pvalues(area),~,stats] = signrank(rate_corr_areas{2,area}(:,2),rate_corr_areas{2,area}(:,3),'tail','right'); % perform one-sided wilcoxon signed-rank test for difference between proximal and distal sessions of spont. activity
    zvalues(area) = stats.zval;
    subplot(2,3,area)
    hold on
    plt(1) = errorbar(mean_rate_corr_movie,ste_rate_corr_movie,'o','color',[0.7 0.7 0.7],...
        'markerfacecolor',[0.7 0.7 0.7],'capsize',0,'linestyle','-','linewidth',2);
    
    plt(2) = errorbar(mean_rate_corr_spont,ste_rate_corr_spont,'o','color',colors(area,:),...
        'markerfacecolor',colors(area,:),'capsize',0,'linestyle','-','linewidth',2);
    
    text(0.75,0.85,brain_areas{area},'Units','normalized','FontSize',15)
    xlim([0.5 3.5])
    if area >3
        xtickangle(15)
        set(gca,'xtick',1:3,'xticklabel',{'Within session','Proximal sessions','Distal sessions'})
    else
        set(gca,'xtick',[])
    end
    if area == 4 || area == 1
        ylabel('Ensemble rate correlation')
    end
    
    legend(plt,{'Nat. Mov.','Spontaneous'},'Location','southwest','fontsize',8)
    legend('boxoff')
    ylim([0.2 max(ylim)])
end

% correct for multiple comparisons using bonferroni-holm method
corrected_pvalues = bonf_holm(pvalues);

% define statistics summary table
VarNames = {'area','zvalue','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas(1:6)'),zvalues',pvalues',corrected_pvalues','VariableNames',VarNames);

% display statistics summary table
clc;
disp(['Difference in the ensemble rate correlation values between proximal'])
disp(['and distal sessions during blocks of spontaneous activity'])
disp(['One-tailed Wilcoxon signed-rank tests with Holm�Bonferroni correction:'])
disp(statistics)


%% Figure S5D - Mean activity rate across sessions

cell_cutoff = 20; % cell count threshold of 20 cells
nat_movie = 1; % natural movie 1

mean_activity_all_mice_area = {};
for area = 1:6 % loop over visual areas
     % valid_mice - true only if the number of cells recorded is
    % atleast 20 in each of the recording sessions
    valid_mice = min(calcium_excitatory_cell_count{area},[],2)>= cell_cutoff;
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie); % subset mice that met the requirements of 'valid_mice'
    
    mean_activity_all_mice = [];
    for mouse = 1:length(current_area) % loop over mice
        current_mouse_activity = current_area{mouse}*30; % subset neuronal activity for a single mouse (#cells by 30 time bins by 30 movie repeats)
        
        % subset neuronal activity only of cells active in session 1
        sess1 = current_mouse_activity(:,:,1:10);
        valid_sess1 = mean(mean(sess1,2,'omitnan'),3,'omitnan')>0;
        
        % subset neuronal activity only of cells active in session 2
        sess2 = current_mouse_activity(:,:,11:20);
        valid_sess2 = mean(mean(sess2,2,'omitnan'),3,'omitnan')>0;
        
        % subset neuronal activity only of cells active in session 3
        sess3 = current_mouse_activity(:,:,21:30);
        valid_sess3 = mean(mean(sess3,2,'omitnan'),3,'omitnan')>0;
        
        mean_sess1 = mean(mean(sess1(valid_sess1,:,:),3,'omitnan'),2,'omitnan'); % calculate the average neuronal activity across movie repeats and time bins for session 1
        mean_sess2 = mean(mean(sess2(valid_sess2,:,:),3,'omitnan'),2,'omitnan'); % calculate the average neuronal activity across movie repeats and time bins for session 2
        mean_sess3 = mean(mean(sess3(valid_sess3,:,:),3,'omitnan'),2,'omitnan'); % calculate the average neuronal activity across movie repeats and time bins for session 3

        mean_activity_all_mice(mouse,:) = [mean(mean_sess1,'omitnan'),mean(mean_sess2,'omitnan'),mean(mean_sess3,'omitnan')]; % calculate average across cells for each session
    end
    mean_activity_all_mice_area{area} = mean_activity_all_mice;
end

pvalues = [];
zvalues = [];
figure('units','normalized','position',[0.3 0.3 0.315 0.35]) % visualize the distribution of average activity rates for each session
for area = 1:6 % loop over areas
    current_area = mean_activity_all_mice_area{area}; % subset activity rates for a single area
    subplot(2,3,area)
    figure_boxplot(current_area);
    title([brain_areas{area},'              ',['N=',num2str(size(current_area,1))]])
    set(gca,'ActivePositionProperty','position')
    xlim([0.1 3.9])
    if area == 4 || area == 1
        ylabel({'Mean activity rate';'(events/sec)'})
    elseif area ==5
        xlabel('Session')
    end
    if area <3
        set(gca,'xtick',[])
    end
    
    % test for differences between sessions
    sig_mat = [];
    z_mat = [];
    for day1 = 1:3 % loop over sessions
        for day2 = 1:3 % loop over sessions
            if day1 ~= day2
                [sig_mat(day1,day2),~,stat] = signrank(current_area(:,day1),current_area(:,day2)); % two-sided Mann-Whitney test for difference between sessions
                z_mat(day1,day2) = stat.zval;
            end
        end
    end
    pvalues(area,:) = [sig_mat(1,2),sig_mat(2,3),sig_mat(1,3)];
    zvalues(area,:) = [z_mat(1,2),z_mat(2,3),z_mat(1,3)];
end

% define statistics summary table
VarNames = {'area','zvalue_sess1_sess2','pvalue_sess1_sess2','zvalue_sess2_sess3','pvalue_sess2_sess3','zvalue_sess1_sess3','pvalue_sess1_sess3'};
statistics = table(cell2mat(brain_areas(1:6)'),zvalues(:,1),pvalues(:,1),zvalues(:,2),pvalues(:,2),zvalues(:,3),pvalues(:,3),'VariableNames',VarNames);

% display statistics summary table
clc;
disp(['Difference in activity rates across sessions'])
disp(['Two-tailed Wilcoxon signed-rank tests with Holm�Bonferroni correction:'])
disp(statistics)

%% Figure S5E - Number of active cells across sessions

cell_cutoff = 20; % cell count threshold of 20 cells
nat_movie = 1; % natural movie 1

mean_num_active_cells_area = {};
for area = 1:6 % loop over visual areas
     % valid_mice - true only if the number of cells recorded is
    % atleast 20 in each of the recording sessions
    valid_mice = min(calcium_excitatory_cell_count{area},[],2)>= cell_cutoff;
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie); % subset mice that met the requirements of 'valid_mice'
    
    mean_num_active_cells = [];
    for mouse = 1:length(current_area) % loop over mice
        current_mouse_activity = current_area{mouse}; % subset neuronal activity for a single mouse (#cells by 30 time bins by 30 movie repeats)
        
        % subset neuronal activity only of cells active in session 1
        sess1 = current_mouse_activity(:,:,1:10);
        valid_sess1 = mean(mean(sess1,2,'omitnan'),3,'omitnan')>0;
        
        % subset neuronal activity only of cells active in session 2
        sess2 = current_mouse_activity(:,:,11:20);
        valid_sess2 = mean(mean(sess2,2,'omitnan'),3,'omitnan')>0;
        
        % subset neuronal activity only of cells active in session 3
        sess3 = current_mouse_activity(:,:,21:30);
        valid_sess3 = mean(mean(sess3,2,'omitnan'),3,'omitnan')>0;
        
        mean_num_active_cells(mouse,:) = [sum(valid_sess1),sum(valid_sess2),sum(valid_sess3)]; % calculate the number of cells found active in each session
    end
    mean_num_active_cells_area{area} = mean_num_active_cells;
end

pvalues = [];
zvalues = [];
figure('units','normalized','position',[0.3 0.3 0.315 0.35]) % visualize the distribution of active cells for each session
for area = 1:6 % loop over areas
    current_area = mean_num_active_cells_area{area};
    subplot(2,3,area)
    figure_boxplot(current_area);
    title([brain_areas{area},'              ',['N=',num2str(size(current_area,1))]])
    set(gca,'ActivePositionProperty','position')
    xlim([0.1 3.9])
    if area == 4 || area == 1
        ylabel({'Number of active cells'})
    elseif area ==5
        xlabel('Session')
    end
    if area <3
        set(gca,'xtick',[])
    end
    
    % test for differences between sessions
    sig_mat = [];
    for day1 = 1:3 % loop over sessions
        for day2 = 1:3 % loop over sessions
            if day1 ~= day2
                [sig_mat(day1,day2),~,stat] =signrank(current_area(:,day1),current_area(:,day2)); % two-sided Mann-Whitney test for difference between sessions
                z_mat(day1,day2) = stat.zval;
            end
        end
    end
    pvalues(area,:) = [sig_mat(1,2),sig_mat(2,3),sig_mat(1,3)];
    zvalues(area,:) = [z_mat(1,2),z_mat(2,3),z_mat(1,3)];
end

% define statistics summary table
VarNames = {'area','zvalue_sess1_sess2','pvalue_sess1_sess2','zvalue_sess2_sess3','pvalue_sess2_sess3','zvalue_sess1_sess3','pvalue_sess1_sess3'};
statistics = table(cell2mat(brain_areas(1:6)'),zvalues(:,1),pvalues(:,1),zvalues(:,2),pvalues(:,2),zvalues(:,3),pvalues(:,3),'VariableNames',VarNames);

% display statistics summary table
clc;
disp(['Difference in number of active cells across sessions'])
disp(['Two-tailed Wilcoxon signed-rank tests with Holm�Bonferroni correction:'])
disp(statistics)


%% Figure S5F - Mean running speed across sessions

cell_cutoff = 20; % cell count threshold of 20 cells
nat_movie = 1; % natural movie 1

mean_running_speed_area = {};
for area = 1:6 % loop over visual areas
     % valid_mice - true only if the number of cells recorded is
    % atleast 20 in each of the recording sessions
    valid_mice = min(calcium_excitatory_cell_count{area},[],2)>= cell_cutoff;
    current_area = calcium_excitatory_running_speed{area}(valid_mice,nat_movie); % subset mice that met the requirements of 'valid_mice'
    
    mean_running_speed = [];
    for mouse = 1:length(current_area ) % loop over mice
        
        current_mouse = current_area{mouse}; % subset running speed of a single mouse (1 by 30 time bins by 30 movie repeats)
       
        sess1 = mean(mean(current_mouse(:,:,1:10),2,'omitnan'),3,'omitnan'); % calculate mean running speed across time bins and movie repeats for session 1
        sess2 = mean(mean(current_mouse(:,:,11:20),2,'omitnan'),3,'omitnan'); % calculate mean running speed across time bins and movie repeats for session 2
        sess3 = mean(mean(current_mouse(:,:,21:30),2,'omitnan'),3,'omitnan'); % calculate mean running speed across time bins and movie repeats for session 3
        
        mean_running_speed(mouse,:) = [sess1,sess2,sess3]; % store mean running speed for all three sessions
    end
    mean_running_speed_area{area,1} = mean_running_speed;
end

mean_running_speed_pooled = cell2mat(mean_running_speed_area); % pool values across all mice and visual areas

figure('units','normalized','position',[0.35 0.35 0.2 0.3]) % visualize the distribution of mean running speed for each session
figure_boxplot(mean_running_speed_pooled);
xlim([0.1 3.9])
ylabel({'Mean running speed (cm/sec)'})
xlabel('Session')

% test for differences between sessions
sig_mat = [];
for day1 = 1:3 % loop over sessions
    for day2 = 1:3 % loop over sessions
        if day1 ~= day2
            [sig_mat(day1,day2),~,stat] =signrank(mean_running_speed_pooled(:,day1),mean_running_speed_pooled(:,day2)); % two-sided Mann-Whitney test for difference between sessions
            z_mat(day1,day2) = stat.zval;
        end
    end
end
pvalues = [sig_mat(1,2),sig_mat(2,3),sig_mat(1,3)];
zvalues = [z_mat(1,2),z_mat(2,3),z_mat(1,3)];
stats = [pvalues;zvalues];

% define statistics summary table
VarNames = {'statistic','sess1_sess2','sess2_sess3','sess1_sess3'};
statistics = table({'pvalue';'zvalue'},stats(:,1),stats(:,2),stats(:,3),'VariableNames',VarNames);

% display statistics summary table
clc;
disp(['Difference in mean running speed across sessions'])
disp(['Two-tailed Wilcoxon signed-rank tests with Holm�Bonferroni correction:'])
disp(statistics)

%% Figure S5G - Mean pupil size across sessions

cell_cutoff = 20; % cell count threshold of 20 cells
nat_movie = 1; % natural movie 1

mean_pupil_size_area = {};
for area = 1:6 % loop over visual areas
     % valid_mice - true only if the number of cells recorded is
    % atleast 20 in each of the recording sessions
    valid_mice = min(calcium_excitatory_cell_count{area},[],2)>= cell_cutoff;
    current_area = calcium_excitatory_pupil_size{area}(valid_mice,nat_movie); % subset mice that met the requirements of 'valid_mice'
    
    mean_pupil_size = [];
    for mouse = 1:length(current_area)  % loop over mice
        current_mouse = current_area{mouse}; % subset pupil area of a single mouse (1 by 30 time bins by 30 movie repeats)
        
        sess1 = mean(mean(current_mouse(:,:,1:10),2,'omitnan'),3,'omitnan'); % calculate mean pupil area across time bins and movie repeats for session 1
        sess2 = mean(mean(current_mouse(:,:,11:20),2,'omitnan'),3,'omitnan'); % calculate mean pupil area across time bins and movie repeats for session 2
        sess3 = mean(mean(current_mouse(:,:,21:30),2,'omitnan'),3,'omitnan'); % calculate mean pupil area across time bins and movie repeats for session 3

        mean_pupil_size(mouse,:) = [sess1,sess2,sess3]; % store mean pupil area for all three sessions
    end
    mean_pupil_size_area{area,1} = mean_pupil_size;
end

mean_pupil_size_pooled = cell2mat(mean_pupil_size_area); % pool values across all mice and visual areas

figure('units','normalized','position',[0.35 0.35 0.2 0.3]) % visualize the distribution of mean pupil area for each session
figure_boxplot(mean_pupil_size_pooled);
xlim([0.1 3.9])
ylabel({'Pupil size (a.u.)'})
xlabel('Session')

% test for differences between sessions
sig_mat = [];
for day1 = 1:3 % loop over sessions
    for day2 = 1:3 % loop over sessions
        if day1 ~= day2
            [sig_mat(day1,day2),~,stat] =signrank(mean_pupil_size_pooled(:,day1),mean_pupil_size_pooled(:,day2)); % two-sided Mann-Whitney test for difference between sessions
            z_mat(day1,day2) = stat.zval;
        end
    end
end
pvalues = [sig_mat(1,2),sig_mat(2,3),sig_mat(1,3)];
zvalues = [z_mat(1,2),z_mat(2,3),z_mat(1,3)];
stats = [pvalues;zvalues];

% define statistics summary table
VarNames = {'statistic','sess1_sess2','sess2_sess3','sess1_sess3'};
statistics = table({'pvalue';'zvalue'},stats(:,1),stats(:,2),stats(:,3),'VariableNames',VarNames);

% display statistics summary table
clc;
disp(['Difference in mean pupil size across sessions'])
disp(['Two-tailed Wilcoxon signed-rank tests with Holm�Bonferroni correction:'])
disp(statistics)

%% Figure S5H - Within day decoder across sessions

cell_cutoff = 20; % cell count threshold of 20 cells
nat_movie = 1; % natural movie 1

within_day_decoder_areas = {};
for area = 1:6 % loop over visual areas
     % valid_mice - true only if the number of cells recorded is
    % atleast 20 in each of the recording sessions
    valid_mice = min(calcium_excitatory_cell_count{area},[],2)>= cell_cutoff;
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie); % subset mice that met the requirements of 'valid_mice'
    
    within_day_decoder = [];
    for mouse = 1:length(current_area) % loop over mice
        clc;
        disp(['Performing time-lapse decoding within sessions:'])
        disp(['Dataset: Calcium imaging | Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse}; % subset neuronal activity for a single mouse (#cells by 30 time bins by 30 movie repeats)
        
        % within session decoding - use the average activity of the first
        % half of movie repeats (repeats 1-5) to predict the time in movie
        % based on the average neuronal activity of the second half of
        % movie repeats (repeats 6-10)
        for day = 1:3 % loop over sessions
            current_day = current_mouse(:,:,[1:10]+10*(day-1)); % subset neuronal activity of a single session
            valid_cells = mean(mean(current_day,3,'omitnan'),2,'omitnan')>0; % find the cells that were active in that session
            
            train_data = mean(current_day(valid_cells,:,1:5),3,'omitnan'); % subset only active cells and calculate the average activity over the first half of movie repeats and use it as train data
            test_data = mean(current_day(valid_cells ,:,6:10),3,'omitnan'); % subset only active cells and calculate the average activity over the second half of movie repeats and use it as test data
            
            mdl = fitcknn(train_data',[1:30]); % fit KNN decoder to predict the time bin based on the activity pattern of each time bin in the train data
            ypred = predict(mdl,test_data'); % predict the time bin identity based on the activity patterns in the test data

            within_day_decoder(mouse,day) = (sum(ypred == [1:30]')./30)*100; % calculate the percentage of correct classifications
        end
    end
    within_day_decoder_areas{area} = within_day_decoder;
end

pvalues = [];
zvalues = [];
figure('units','normalized','position',[0.3 0.3 0.315 0.35]) % visualize the distribution of within session decoder performance
for area = 1:6 % loop over areas
    current_area = within_day_decoder_areas{area}; % subset decoder performance of a single area
    subplot(2,3,area)
    figure_boxplot(current_area);
    xlim([0.1 3.9])
    ylim([0 100])
    hold on
    plot(xlim,[100 100]./30,'--','color',[0.2 0.2 0.2],'linewidth',1)
    title([brain_areas{area},'              ',['N=',num2str(size(current_area,1))]])
    set(gca,'ActivePositionProperty','position')
    
    if area == 4 || area == 1
        ylabel({'Successful';'classifications (%)'})
    elseif area ==5
        xlabel('Session')
    end
    if area <3
        set(gca,'xtick',[])
    end
    
    % test for differences between sessions
    sig_mat = [];
    for day1 = 1:3 % loop over sessions
        for day2 = 1:3 % loop over sessions
            if day1 ~= day2
                [sig_mat(day1,day2),~,stat] =signrank(current_area(:,day1),current_area(:,day2)); % two-sided Mann-Whitney test for difference between sessions
                z_mat(day1,day2) = stat.zval;
            end
        end
    end
    pvalues(area,:) = [sig_mat(1,2),sig_mat(2,3),sig_mat(1,3)];
    zvalues(area,:) = [z_mat(1,2),z_mat(2,3),z_mat(1,3)];
end

% define statistics summary table
VarNames = {'area','zvalue_sess1_sess2','pvalue_sess1_sess2','zvalue_sess2_sess3','pvalue_sess2_sess3','zvalue_sess1_sess3','pvalue_sess1_sess3'};
statistics = table(cell2mat(brain_areas(1:6)'),zvalues(:,1),pvalues(:,1),zvalues(:,2),pvalues(:,2),zvalues(:,3),pvalues(:,3),'VariableNames',VarNames);

% display statistics summary table
clc;
disp(['Difference in within day decoder performance across sessions'])
disp(['Two-tailed Wilcoxon signed-rank tests with Holm�Bonferroni correction:'])
disp(statistics)

%% Figure S5I - With day PV correlation across sessions

cell_cutoff = 20; % cell count threshold of 20 cells
nat_movie = 1; % natural movie 1

within_day_pv_areas = {};
for area = 1:6 % loop over visual areas
    % valid_mice - true only if the number of cells recorded is
    % atleast 20 in each of the recording sessions
    valid_mice = min(calcium_excitatory_cell_count{area},[],2)>= cell_cutoff;
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie); % subset mice that met the requirements of 'valid_mice'
    
    within_day_pv = [];
    for mouse = 1:length(current_area) % loop over mice
        clc;
        disp(['Calculation population vector correlation within sessions:'])
        disp(['Dataset: Calcium imaging | Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse}; % subset neuronal activity for a single mouse (#cells by 30 time bins by 30 movie repeats)
        for day = 1:3 % loop over sessions
            current_day = current_mouse(:,:,[1:10]+10*(day-1)); % subset neuronal activity of a single session
            valid_cells = mean(mean(current_day,3,'omitnan'),2,'omitnan')>0; % find the cells that were active in that session
            
            activity_first_half = mean(current_day(valid_cells,:,1:5),3,'omitnan'); % subset only active cells and calculate the average activity over the first half of movie repeats and use it as train data
            activity_second_half = mean(current_day(valid_cells ,:,6:10),3,'omitnan'); % subset only active cells and calculate the average activity over the second half of movie repeats and use it as test data
             within_day_pv(mouse,day) = mean(diag(corr(activity_first_half,activity_second_half)),'omitnan'); % calculate the mean pv correlation across corresponding time bins between first and second half of movie repeats
        end
    end
    within_day_pv_areas{area} = within_day_pv;
end

pvalues = [];
zvalues = [];
figure('units','normalized','position',[0.3 0.3 0.315 0.35]) % visualize the distribution of within session decoder performance
for area = 1:6 % loop over areas
    current_area = within_day_pv_areas{area}; % subset pv corr values of a single area
    subplot(2,3,area)
    figure_boxplot(current_area);
    xlim([0.1 3.9])
    
    title([brain_areas{area},'              ',['N=',num2str(size(current_area,1))]])
    set(gca,'ActivePositionProperty','position')
    
    if area == 4 || area == 1
        ylabel({'Within session';'PV correlation'})
    elseif area ==5
        xlabel('Session')
    end
    if area <3
        set(gca,'xtick',[])
    end
    
    % test for differences between sessions
    sig_mat = [];
    for day1 = 1:3 % loop over sessions
        for day2 = 1:3 % loop over sessions
            if day1 ~= day2
                [sig_mat(day1,day2),~,stat] =signrank(current_area(:,day1),current_area(:,day2)); % two-sided Mann-Whitney test for difference between sessions
                z_mat(day1,day2) = stat.zval;
            end
        end
    end
    pvalues(area,:) = [sig_mat(1,2),sig_mat(2,3),sig_mat(1,3)];
    zvalues(area,:) = [z_mat(1,2),z_mat(2,3),z_mat(1,3)];
end

% define statistics summary table
VarNames = {'area','zvalue_sess1_sess2','pvalue_sess1_sess2','zvalue_sess2_sess3','pvalue_sess2_sess3','zvalue_sess1_sess3','pvalue_sess1_sess3'};
statistics = table(cell2mat(brain_areas(1:6)'),zvalues(:,1),pvalues(:,1),zvalues(:,2),pvalues(:,2),zvalues(:,3),pvalues(:,3),'VariableNames',VarNames);

% display statistics summary table
clc;
disp(['Difference in within day PV correlation across sessions'])
disp(['Two-tailed Wilcoxon signed-rank tests with Holm�Bonferroni correction:'])
disp(statistics)

%% Figure S5J - PV correlation difference between pairs of subsequent sessions

cell_cutoff = 20; % cell count threshold of 20 cells
nat_movie = 1; % natural movie 1

pvalue = [];
zvalue = [];
pv_corr_diff_across_areas = nan(100,6);
for area = 1:6 % loop over visual areas
    % valid_mice - true only if the number of cells recorded is
    % atleast 20 in each of the recording sessions
    valid_mice = min(calcium_excitatory_cell_count{area},[],2)>= cell_cutoff;
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie); % subset mice that met the requirements of 'valid_mice'
    
    elapsed_pv_corr = [];
    for mouse = 1:length(current_area) % loop over mice
        clc;
        disp(['Calculating population vector correlation between ssessions:'])
        disp(['Dataset: Calcium imaging | Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse}; % subset neuronal activity for a single mouse (#cells by 30 time bins by 30 movie repeats)
        
        % average the neuronal activity across movie repeats for each session
        mean_activity_per_sess = [];
        for sess = 1:3
            current_sess = [1:10] + 10*(sess-1); % range of movie repeats for a given session
            mean_activity_per_sess(:,:,sess) = mean(current_mouse(:,:,current_sess),3,'omitnan'); % average the neuronal activity across range of movie repeats
        end
        
        pv_corr = [];
        for sessA = 1:size(mean_activity_per_sess,3) % loop over sessions
            sessA_activity = mean_activity_per_sess(:,:,sessA); % subset neuronal activity of a single session
            for sessB = 1:size(mean_activity_per_sess,3) % loop over sessions
                sessB_activity = mean_activity_per_sess(:,:,sessB); % subset neuronal activity of a single session
                
                valid_cells = [mean(sessA_activity,2,'omitnan') ~= 0] & [mean(sessB_activity,2,'omitnan') ~= 0]; % include only cells that are active in both compared sessions
                valid_sessA_activity = sessA_activity(valid_cells,:); % subset neuronal activity only for cells that area active in both compared sessions
                valid_sessB_activity = sessB_activity(valid_cells,:); % subset neuronal activity only for cells that area active in both compared sessions
    
                pv_corr(sessA,sessB) = mean(diag(corr(valid_sessA_activity,valid_sessB_activity)),'omitnan'); % calculate the mean pv corr across corresponding time bins between pair of sessions
            end
        end
        
        elapsed_pv_corr(mouse,:) = [diag(pv_corr,1)]; % subset the pv corr values of subsequent pairs of sessions
    end
  
    pv_corr_diff = elapsed_pv_corr(:,1) - elapsed_pv_corr(:,2); % calculate the pv corr difference between different pairs of subsequent sessions
    pv_corr_diff_across_areas(1:mouse,area) = pv_corr_diff; % store pv corr difference values
    [pvalue(area),~,stat] = signrank(elapsed_pv_corr(:,1),elapsed_pv_corr(:,2)); % perform two-sided Wilcoxon signed-rank test for pv corr difference between different pairs of subsequent sessions
    zvalue(area) = stat.zval;
end

figure('units','normalized','position',[0.3 0.3 0.2 0.3]) % visualize the distribution of pv correlation difference between different pairs of subsequent sessions
xlim([0 7])
ylim([-0.6 0.6])
hold on
figure_boxplot(pv_corr_diff_across_areas);
ylabel('PV correlation difference')
set(gca,'xtick',1:6,'xticklabel',brain_areas(1:6))

% define statistics summary table
VarNames = {'area','zvalue','pvalue'};
statistics = table(cell2mat(brain_areas(1:6)'),zvalue',pvalue','VariableNames',VarNames);

% display statistics summary table
clc;
disp(['Difference in the PV correlation values between pairs of subsequent sessions'])
disp(['Two-tailed Wilcoxon signed-rank tests with Holm�Bonferroni correction:'])
disp(statistics)


%% Figure S5K - PV correlation across sessions using dF/F traces

nat_movie = 1; % natural movie 1
cell_cutoff = 20; % cell count threshold of 20 cells

pv_corr_areas = {}; % define empty variable that will store the pv corr values of within and between sessions across mice and visual areas
for area = 1:6% loop over visual areas
    % valid_mice - true only if the number of cells recorded is
    % atleast 20 in each of the recording sessions
    valid_mice = min(calcium_excitatory_cell_count{area},[],2)>= cell_cutoff;
    
    % subset mice that passed the requirments of 'valid_mice'
    current_area_events = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie); % mice that passed the requirments of 'valid_mice'
    current_area_raw = calcium_excitatory_population_vectors_raw{area}(valid_mice);
    
    elapsed_sess_pv_corr_raw = []; % define empty variable that will store the pv correlation values for all mice using raw df/f0
    elapsed_sess_pv_corr_events = []; % define empty variable that will store the pv correlation values for all mice using event detection
    for mouse = 1:length(current_area_events) % loop over mice
        clc;
        disp(['Calculating population vector correlation between ssessions:'])
        disp(['Dataset: Calcium imaging | Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area_events))])
        
        current_mouse_events = current_area_events{mouse}; % subset neuronal activity (events) for a single mouse (#cells by 30 time bins by 30 movie repeats)
        current_mouse_raw = current_area_raw{mouse}; % subset neuronal activity (df/f0) for a single mouse (#cells by 30 time bins by 30 movie repeats)
        
        pv_corr_events = []; % define empty variable that will store pv corr values for event detection
        pv_corr_raw= []; % define empty variable that will store pv corr values for df/f0
        for half1 = 1:6 % loop over session halves
            current_half1_events = mean(current_mouse_events(:,:,[1:5]+5*(half1-1)),3,'omitnan'); % average the neuronal activity (events) over movie repeats for a single session half
            current_half1_raw = mean(current_mouse_raw(:,:,[1:5]+5*(half1-1)),3,'omitnan'); % average the neuronal activity (df/f0) over movie repeats for a single session half
            for  half2 = 1:6 % loop over session halves
                current_half2_events = mean(current_mouse_events(:,:,[1:5]+5*(half2-1)),3,'omitnan'); % average the neuronal activity (events) over movie repeats for a single session half
                current_half2_raw = mean(current_mouse_raw(:,:,[1:5]+5*(half2-1)),3,'omitnan'); % average the neuronal activity (df/f0) over movie repeats for a single session half

                % include only cells that were found active in both compared time points
                valid_cells = [nanmean(current_half1_events,2) > 0] & [nanmean(current_half2_events,2) > 0]; 
                
                valid_half1_activity_events = current_half1_events(valid_cells,:); 
                valid_half2_activity_events = current_half2_events(valid_cells,:); 
                
                valid_half1_activity_raw = current_half1_raw(valid_cells,:);
                valid_half2_activity_raw = current_half2_raw(valid_cells,:);
                
                % calculate the mean PV corr across corresponding time bins
                % between sessions hlaves
                pv_corr_events(half1,half2) = mean(diag(corr(valid_half1_activity_events,valid_half2_activity_events)),'omitnan'); % using event detection
                pv_corr_raw(half1,half2) = mean(diag(corr(valid_half1_activity_raw,valid_half2_activity_raw)),'omitnan'); % using raw df/f0
            end
        end
        
        % for event detection data
        within_session_values_events = [pv_corr_events(1,2),pv_corr_events(3,4),pv_corr_events(5,6)]; % pv corr values between halves of the same session (within session)
        proximal_sessions_values_events = [pv_corr_events(1:2,3:4),pv_corr_events(3:4,5:6)]; % pv corr values between halves of proximal sessions (sessions 1&2 and sessions 2&3)
        distal_session_values_events = pv_corr_events(1:2,5:6); % pv corr values between halves of distal sessions (sessions 1&3)
        
        elapsed_sess_pv_corr_events(mouse,1) = mean(within_session_values_events(:),'omitnan'); % average pv corr values for within session
        elapsed_sess_pv_corr_events(mouse,2) = mean(proximal_sessions_values_events(:),'omitnan'); % average pv corr values for proximal sessions
        elapsed_sess_pv_corr_events(mouse,3) = mean(distal_session_values_events(:),'omitnan'); % average pv corr values for distal sessions
        
     
        % for raw df/f0 data
        within_session_values_raw = [pv_corr_raw(1,2),pv_corr_raw(3,4),pv_corr_raw(5,6)]; % pv corr values between halves of the same session (within session)
        proximal_sessions_values_raw = [pv_corr_raw(1:2,3:4),pv_corr_raw(3:4,5:6)]; % pv corr values between halves of proximal sessions (sessions 1&2 and sessions 2&3)
        distal_session_values_raw = pv_corr_raw(1:2,5:6); % pv corr values between halves of distal sessions (sessions 1&3)
        
        elapsed_sess_pv_corr_raw(mouse,1) = mean(within_session_values_raw(:),'omitnan'); % average pv corr values for within session
        elapsed_sess_pv_corr_raw(mouse,2) = mean(proximal_sessions_values_raw(:),'omitnan'); % average pv corr values for proximal sessions
       elapsed_sess_pv_corr_raw(mouse,3) = mean(distal_session_values_raw(:),'omitnan'); % average pv corr values for distal sessions

        
    end

    pv_corr_areas{1,area} = elapsed_sess_pv_corr_events;
    pv_corr_areas{2,area} = elapsed_sess_pv_corr_raw;
    
end

figure('units','normalized','position',[0.3 0.3 0.35 0.35])
pvalues = [];
zvalues = [];
plt = [];
for area = 1:6
    subplot(2,3,area)
    % for event detection data
    mean_pv_corr_events = mean(pv_corr_areas{1,area},'omitnan'); % calculate mean pv across mice
    std_pv_corr_events = std(pv_corr_areas{1,area})./sqrt(size(pv_corr_areas{1,area},1)); % calculate standard error across mice
    
    % for raw df/f0 data
    mean_pv_corr_raw = mean(pv_corr_areas{2,area},'omitnan'); % calculate mean pv across mice
    std_pv_corr_raw = std(pv_corr_areas{2,area})./sqrt(size(pv_corr_areas{2,area},1)); % calculate standard error across mice
    
    
    [pvalues(area),~,stats] = signrank(pv_corr_areas{2,area}(:,2),pv_corr_areas{2,area}(:,3),'tail','right'); % perform one-sided wilcoxon signed-rank test for difference between proximal and distal sessions for df/f0
    zvalues(area) = stats.zval;
    subplot(2,3,area)
    hold on
    plt(1) = errorbar(mean_pv_corr_raw,std_pv_corr_raw,'o','color',[0.7 0.7 0.7],...
        'markerfacecolor',[0.7 0.7 0.7],'capsize',0,'linestyle','-','linewidth',2);
    
    plt(2) = errorbar(mean_pv_corr_events,std_pv_corr_events,'o','color',colors(area,:),...
        'markerfacecolor',colors(area,:),'capsize',0,'linestyle','-','linewidth',2);
    
    text(0.1,0.15,brain_areas{area},'Units','normalized','FontSize',15)
    xlim([0.5 3.5])
    if area >3
        xtickangle(15)
        set(gca,'xtick',1:3,'xticklabel',{'Within session','Proximal sessions','Distal sessions'})
    else
        set(gca,'xtick',[])
    end
    if area == 4 || area == 1
        ylabel('PV correlation')
    end
    
    legend(plt,{'dF(t)/F0','Events'},'Location','best','fontsize',8)
    legend('boxoff')
end

% define statistics summary table
VarNames = {'area','zvalue','pvalue'};
statistics = table(cell2mat(brain_areas(1:6)'),zvalues',pvalues','VariableNames',VarNames);

% display statistics summary table
clc;
disp(['Difference in the PV correlation values between proximal'])
disp(['and distal sessions when using dF(t)/F0 traces'])
disp(['One-tailed Wilcoxon signed-rank tests:'])
disp(statistics)


%% Figure S5L - PV correlation across sessions using dF/F traces (scatter plot)

nat_movie = 1; % natural movie 1
cell_cutoff = 20; % cell count threshold of 20 cells

pv_corr_areas = {}; % define empty variable that will store the pv corr values of within and between sessions across mice and visual areas
for area = 1:6% loop over visual areas
    % valid_mice - true only if the number of cells recorded is
    % atleast 20 in each of the recording sessions
    valid_mice = min(calcium_excitatory_cell_count{area},[],2)>= cell_cutoff;
    
    % subset mice that passed the requirments of 'valid_mice'
    current_area_events = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie); % mice that passed the requirments of 'valid_mice'
    current_area_raw = calcium_excitatory_population_vectors_raw{area}(valid_mice);
    
    elapsed_sess_pv_corr_raw = []; % define empty variable that will store the pv correlation values for all mice using raw df/f0
    elapsed_sess_pv_corr_events = []; % define empty variable that will store the pv correlation values for all mice using event detection
    for mouse = 1:length(current_area_events) % loop over mice
        clc;
        disp(['Calculating population vector correlation between ssessions:'])
        disp(['Dataset: Calcium imaging | Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area_events))])
        
        current_mouse_events = current_area_events{mouse}; % subset neuronal activity (events) for a single mouse (#cells by 30 time bins by 30 movie repeats)
        current_mouse_raw = current_area_raw{mouse}; % subset neuronal activity (df/f0) for a single mouse (#cells by 30 time bins by 30 movie repeats)
        
        pv_corr_events = []; % define empty variable that will store pv corr values for event detection
        pv_corr_raw= []; % define empty variable that will store pv corr values for df/f0
        for half1 = 1:6 % loop over session halves
            current_half1_events = mean(current_mouse_events(:,:,[1:5]+5*(half1-1)),3,'omitnan'); % average the neuronal activity (events) over movie repeats for a single session half
            current_half1_raw = mean(current_mouse_raw(:,:,[1:5]+5*(half1-1)),3,'omitnan'); % average the neuronal activity (df/f0) over movie repeats for a single session half
            for  half2 = 1:6 % loop over session halves
                current_half2_events = mean(current_mouse_events(:,:,[1:5]+5*(half2-1)),3,'omitnan'); % average the neuronal activity (events) over movie repeats for a single session half
                current_half2_raw = mean(current_mouse_raw(:,:,[1:5]+5*(half2-1)),3,'omitnan'); % average the neuronal activity (df/f0) over movie repeats for a single session half

                % include only cells that were found active in both compared time points
                valid_cells = [nanmean(current_half1_events,2) > 0] & [nanmean(current_half2_events,2) > 0]; 
                
                valid_half1_activity_events = current_half1_events(valid_cells,:); 
                valid_half2_activity_events = current_half2_events(valid_cells,:); 
                
                valid_half1_activity_raw = current_half1_raw(valid_cells,:);
                valid_half2_activity_raw = current_half2_raw(valid_cells,:);
                
                % calculate the mean PV corr across corresponding time bins
                % between sessions hlaves
                pv_corr_events(half1,half2) = mean(diag(corr(valid_half1_activity_events,valid_half2_activity_events)),'omitnan'); % using event detection
                pv_corr_raw(half1,half2) = mean(diag(corr(valid_half1_activity_raw,valid_half2_activity_raw)),'omitnan'); % using raw df/f0
            end
        end
        
        % for event detection data
        within_session_values_events = [pv_corr_events(1,2),pv_corr_events(3,4),pv_corr_events(5,6)]; % pv corr values between halves of the same session (within session)
        proximal_sessions_values_events = [pv_corr_events(1:2,3:4),pv_corr_events(3:4,5:6)]; % pv corr values between halves of proximal sessions (sessions 1&2 and sessions 2&3)
        distal_session_values_events = pv_corr_events(1:2,5:6); % pv corr values between halves of distal sessions (sessions 1&3)
        
        elapsed_sess_pv_corr_events(mouse,1) = mean(within_session_values_events(:),'omitnan'); % average pv corr values for within session
        elapsed_sess_pv_corr_events(mouse,2) = mean(proximal_sessions_values_events(:),'omitnan'); % average pv corr values for proximal sessions
        elapsed_sess_pv_corr_events(mouse,3) = mean(distal_session_values_events(:),'omitnan'); % average pv corr values for distal sessions
        
     
        % for raw df/f0 data
        within_session_values_raw = [pv_corr_raw(1,2),pv_corr_raw(3,4),pv_corr_raw(5,6)]; % pv corr values between halves of the same session (within session)
        proximal_sessions_values_raw = [pv_corr_raw(1:2,3:4),pv_corr_raw(3:4,5:6)]; % pv corr values between halves of proximal sessions (sessions 1&2 and sessions 2&3)
        distal_session_values_raw = pv_corr_raw(1:2,5:6); % pv corr values between halves of distal sessions (sessions 1&3)
        
        elapsed_sess_pv_corr_raw(mouse,1) = mean(within_session_values_raw(:),'omitnan'); % average pv corr values for within session
        elapsed_sess_pv_corr_raw(mouse,2) = mean(proximal_sessions_values_raw(:),'omitnan'); % average pv corr values for proximal sessions
       elapsed_sess_pv_corr_raw(mouse,3) = mean(distal_session_values_raw(:),'omitnan'); % average pv corr values for distal sessions

        
    end

    pv_corr_areas{1,area} = elapsed_sess_pv_corr_events;
    pv_corr_areas{2,area} = elapsed_sess_pv_corr_raw;
    
end

figure('units','normalized','position',[0.3 0.3 0.35 0.35])
for area = 1:6 % loop over areas
    subplot(2,3,area)
    xlim([-0.2 1])
    ylim([-0.2 1])
    hold on
    scatter(pv_corr_areas{1,area}(:),pv_corr_areas{2,area}(:),10,[0.7 0.7 0.7],'filled','MarkerfaceAlpha',0.8)
    
    plot([-0.2 1],[-0.2 1],'--','color',[0.2 0.2 0.2],'linewidth',1)
    
    [r,p] = corr(pv_corr_areas{1,area}(:),pv_corr_areas{2,area}(:));
    text(0.1, 0.85,['r = ',num2str(round(r,2))],'Units','normalized','color',[0.25 0.25 0.25],'Fontsize',10)
    text(0.1, 0.95,[brain_areas{area}],'Units','normalized','color',[0.25 0.25 0.25],'Fontsize',12)
    
    if area == 5
        xlabel('PV correlation using event detection')
    elseif area == 4
        ylabel('PV correlation using raw dF/F')
    end
    
end

%% Figure S6A-B - ROIs across sessions for a single LM representative mouse

mouse = 14; % example mouse #14

% load neuronal and registration data for example mouse from area LM
results_path = ['E:\daniel_master\AllenBrainObservatory\calcium_imaging\results_files\excitatory4\VISl\'];
mat_list = dir([results_path,'*.mat']);
mat_list = {mat_list.name};
registration_path = ['D:\daniel-master\AllenBrainObservatory\Advanced\experimental_data\VISl\'];

load([registration_path,mat_list{mouse}(1:end-4),'\registration\aligned_data_struct.mat']) 

reg_file_name = dir([registration_path,mat_list{mouse}(1:end-4),'\registration\cellRegistered*.mat']);
quality_file_name = dir([registration_path,mat_list{mouse}(1:end-4),'\registration\','*.txt']);

load([registration_path,mat_list{mouse}(1:end-4),'\registration\aligned_data_struct.mat'])
load([registration_path,mat_list{mouse}(1:end-4),'\registration\',reg_file_name(2).name])
load([results_path,mat_list{mouse}],'mouse_age','cell_registration')

cellreg = cell_registered_struct.cell_to_index_map;
[~,exp_order] = sort(mouse_age);

sorted_spatial_footprints_across_days = [];
for session = 1:3
    experiment = exp_order(session);
    exp_list = dir([registration_path,mat_list{mouse}(1:end-4)]);
    exp_list(ismember({exp_list.name}, {'.', '..','registration'})) = [];
    myDir = find(vertcat(exp_list.isdir));
    exp_list = {exp_list(myDir).name};
    
    
    ROI_file_names  = dir([registration_path,mat_list{mouse}(1:end-4),'\',exp_list{experiment},'\ROIs\*.csv']);
    ROI_file_names = {ROI_file_names.name};
    currented_labels_sess = cell2mat(cellfun(@str2num,cellfun(@(x) x(1:end-4),ROI_file_names,'UniformOutput',false),'UniformOutput',false)');
    [~,ai,bi] = intersect(currented_labels_sess,cell_registration{experiment});
    
    current_sess_centeroid = aligned_data_struct.centroid_locations_corrected{1,experiment};
    current_sess_foorprints = aligned_data_struct.spatial_footprints_corrected{1,experiment};
    
    sorted_reg = cell_registration{experiment}(bi);
    for cell = 1:size(cellreg,1)
        clc;[mouse,session,cell]
        ordered_cell = find(cellreg(:,experiment) == cell);
        if ~isempty(ordered_cell)
            
            current_foot_print = squeeze(current_sess_foorprints(cell,:,:));
            sorted_spatial_footprints_across_days(:,:,ordered_cell,session) = current_foot_print;
        end
    end
    
end

ROI_overlay = [];
ROI_overlay(:,:,1) =  sum(sorted_spatial_footprints_across_days(:,:,:,1),3)>0;
ROI_overlay(:,:,2) =  sum(sorted_spatial_footprints_across_days(:,:,:,2),3)>0;
ROI_overlay(:,:,3) =  sum(sorted_spatial_footprints_across_days(:,:,:,3),3)>0;


figure('Units','normalized','position',[0.1 0.3 0.85 0.3])
subplot(1,4,1)
session1_FOV = ROI_overlay;
session1_FOV(:,:,2:3) = 0;
imagesc(session1_FOV)
title(['Session 1 (Day ',num2str(mouse_age(exp_order(1))),') - ',num2str(length(cell_registration{exp_order(1)})),' cells'])
set(gca,'xtick',[],'ytick',[])

subplot(1,4,2)
session2_FOV = ROI_overlay;
session2_FOV(:,:,[1,3]) = 0;
imagesc(session2_FOV)
title(['Session 2 (Day ',num2str(mouse_age(exp_order(2))),') - ',num2str(length(cell_registration{exp_order(2)})),' cells'])
set(gca,'xtick',[],'ytick',[])

subplot(1,4,3)
session3_FOV = ROI_overlay;
session3_FOV(:,:,1:2) = 0;
imagesc(session3_FOV)
title(['Session 3 (Day ',num2str(mouse_age(exp_order(3))),') - ',num2str(length(cell_registration{exp_order(3)})),' cells'])
set(gca,'xtick',[],'ytick',[])

subplot(1,4,4)
imagesc(ROI_overlay)
title(['RGB overlay across sessions:'])
set(gca,'xtick',[],'ytick',[])

%% Figure S6C - Responsiveness of three selected neurons across sessions
nat_movie = 1;
area = 2;
mouse = 14;

results_path = ['E:\daniel_master\AllenBrainObservatory\calcium_imaging\results_files\excitatory4\VISl\'];
mat_list = dir([results_path,'*.mat']);
mat_list = {mat_list.name};
registration_path = ['D:\daniel-master\AllenBrainObservatory\Advanced\experimental_data\VISl\'];

load([registration_path,mat_list{mouse}(1:end-4),'\registration\aligned_data_struct.mat'])

reg_file_name = dir([registration_path,mat_list{mouse}(1:end-4),'\registration\cellRegistered*.mat']);
quality_file_name = dir([registration_path,mat_list{mouse}(1:end-4),'\registration\','*.txt']);

load([registration_path,mat_list{mouse}(1:end-4),'\registration\aligned_data_struct.mat'])
load([registration_path,mat_list{mouse}(1:end-4),'\registration\',reg_file_name(2).name])
load([results_path,mat_list{mouse}],'mouse_age','cell_registration','filtered_traces_days_events')

current_mouse = filtered_traces_days_events(:,1);
cellreg = cell_registered_struct.cell_to_index_map;
[~,exp_order] = sort(mouse_age);

sorted_traces_across_days = zeros(size(cellreg,1),9000,3);
for session = 1:3
    experiment = exp_order(session);
    exp_list = dir([registration_path,mat_list{mouse}(1:end-4)]);
    exp_list(ismember({exp_list.name}, {'.', '..','registration'})) = [];
    myDir = find(vertcat(exp_list.isdir));
    exp_list = {exp_list(myDir).name};
    
    
    ROI_file_names  = dir([registration_path,mat_list{mouse}(1:end-4),'\',exp_list{experiment},'\ROIs\*.csv']);
    ROI_file_names = {ROI_file_names.name};
    currented_labels_sess = cell2mat(cellfun(@str2num,cellfun(@(x) x(1:end-4),ROI_file_names,'UniformOutput',false),'UniformOutput',false)');
    [~,ai,bi] = intersect(currented_labels_sess,cell_registration{experiment});
    current_sess = current_mouse{experiment}(bi,:);
    current_sess_centeroid = aligned_data_struct.centroid_locations_corrected{1,experiment};
    current_sess_foorprints = aligned_data_struct.spatial_footprints_corrected{1,experiment};
    
    sorted_reg = cell_registration{experiment}(bi);
    for cell = 1:size(currented_labels_sess,1)
        clc;
        disp(['Registration of neuronal responses across sessions using Sheintuch et al. (2017) method:'])
        disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area}, ' | Session: ',num2str(session),'/3 | Cell: ',num2str(cell),'/',num2str(size(currented_labels_sess,1))])
        ordered_cell = find(cellreg(:,experiment) == cell);
        if ~isempty(ordered_cell)
            sorted_traces_across_days(ordered_cell,:,session) = current_sess(cell,:);
            sorted_centeroids_across_days(ordered_cell,:,session) = current_sess_centeroid(cell,:);
            current_foot_print = squeeze(current_sess_foorprints(cell,:,:));
            sorted_spatial_footprints_across_days(:,:,ordered_cell,session) = current_foot_print;
            cell_reg_allen_reg(ordered_cell) = sorted_reg(cell);
        end
    end
    
    
end

frames_rate = 30;
repeats_movie1 = 10;
movie1_length = 30;
movie1_frames = frames_rate * movie1_length;
movie1_bin_size = 30;
binned_movie1 = ones(movie1_frames,1);
movie1_bin_edges = 1:movie1_bin_size:movie1_frames;
for bin = 1:length(movie1_bin_edges)
    binned_movie1(movie1_bin_edges(bin):movie1_bin_edges(bin)+movie1_bin_size-1) = bin;
end
binned_movie_repeated1 = repmat(binned_movie1,[repeats_movie1,1]);

% binization for natural movies
pop_vector_info_trials = [];
sub  = 1;
for session = 1:3
    for repeat = 1:10
        frames_temp = [1:movie1_frames] + (movie1_frames*(repeat-1));
        current_repeat = sorted_traces_across_days(:,frames_temp,session);
        
        for bin = 1:length(movie1_bin_edges)
            pop_vector_info_trials(:,bin,sub) = nanmean(current_repeat(:,binned_movie_repeated1(frames_temp) == bin),2);
            
        end
        sub = sub + 1;
    end
end

cell_list = [57,20,81];
chosen_cells_roi = nanmean(sorted_spatial_footprints_across_days(:,:,cell_list,1),3);
chosen_cells_roi(chosen_cells_roi>0) = 1;

figure
imagesc(chosen_cells_roi)
colormap(gray)

figure
for cell = 1:3
    current_cell = squeeze(pop_vector_info_trials(cell_list(cell),:,:))';
    smooth_current_cell = [];
    for repeat = 1:size(current_cell,1)
        smooth_current_cell(repeat,:) = imgaussfilt(current_cell(repeat,:),2);
    end
    
    norm_current_cell = smooth_current_cell ./ max(smooth_current_cell,[],2);
    [row,col] = find(norm_current_cell == 1);
    [B,I] = sort(row);
    
    mean_current_cell = [];
    mean_current_cell(1,:) = nanmean(current_cell(1:10,:));
    mean_current_cell(2,:) = nanmean(current_cell(11:20,:));
    mean_current_cell(3,:) = nanmean(current_cell(21:30,:));
    
    
    norm_current_cell(1:10,:) = norm_current_cell(1:10,:) *  max(mean_current_cell(1,:));
    norm_current_cell(11:20,:) = norm_current_cell(11:20,:) *  max(mean_current_cell(2,:));
    norm_current_cell(21:30,:) = norm_current_cell(21:30,:) *  max(mean_current_cell(3,:));
    
    
    subplot(2,3,cell)
    imagesc(norm_current_cell)
    hold on
    plot(xlim, [10 10]+0.5,'--','linewidth',2,'color','w')
    plot(xlim, [20 20]+0.5,'--','linewidth',2,'color','w')
    if cell == 1
        ylabel('Movie repeat')
    end
    colormap(newmap3)
    title(['Cell #',num2str(cell)])
    
    subplot(2,3,cell+3)
    
    hold on
    plot(mean_current_cell(1,:)*30,'-','markersize',20,'linewidth',1.5,'color',[0.7 0 0])
    plot(mean_current_cell(2,:)*30,'-','markersize',20,'linewidth',1.5,'color',[0 0.7 0])
    plot(mean_current_cell(3,:)*30,'-','markersize',20,'linewidth',1.5,'color',[0 0 0.7])
    
    xlim([1 30])
    if cell == 1
        ylabel({'Mean activity rate';'(events/sec)'})
        
    elseif cell == 2
        xlabel('Time in movie (sec)')
    elseif cell == 3
        legend({'Session 1','Session 2','Session 3'},'location','best','Fontsize',8)
        legend('boxoff')
    end
    
end

%% Cell registration contol - loading registration data and calculating the different measurments
% IMPORTANT! You must run this section before attempting to run Figures S6E-J!

area = 2;
nat_movie = 1;
cell_cutoff = 20;
valid_mice = min(calcium_excitatory_cell_count{area},[],2) >= cell_cutoff;
current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie);

results_path = ['E:\daniel_master\AllenBrainObservatory\calcium_imaging\results_files\excitatory4\VISl\'];
mat_list = dir([results_path,'*.mat']);
mat_list = {mat_list.name};
mat_list = mat_list(valid_mice);

registration_path = ['D:\daniel-master\AllenBrainObservatory\Advanced\experimental_data\VISl\'];

allen_elapsed_session = [];
cellreg_elapsed_session = [];

num_cells_across_registrations =[];
sorted_elapsed_time_all_mice = [];
quality_across_mice = [];

all_active_cells_across_mice = [];
active_both_cells_across_mice = [];
cellreg_all_active_cells_across_mice = [];
cellreg_active_both_cells_across_mice = [];

allen_pv_elapsed_full_sess = [];
cellreg_pv_elapsed_full_sess=[];

pairs_footprint_same_cells_across_mice = {};
pairs_footprint_diff_cells_across_mice = {};
sorted_elapsed_time_per_cell_across_mice = {};

for mouse  = 1:length(current_area)
    clc;
    disp(['Registration of neuronal responses across sessions using Sheintuch et al. (2017) method:'])
    disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
    
    current_mouse = current_area{mouse};
    mean_current_mouse_sessions = [];
    mean_current_mouse_sessions(:,:,1) =  nanmean(current_mouse(:,:,1:5),3);
    mean_current_mouse_sessions(:,:,2)  =  nanmean(current_mouse(:,:,6:10),3);
    mean_current_mouse_sessions(:,:,3)  =  nanmean(current_mouse(:,:,11:15),3);
    mean_current_mouse_sessions(:,:,4) =  nanmean(current_mouse(:,:,16:20),3);
    mean_current_mouse_sessions(:,:,5)  =  nanmean(current_mouse(:,:,21:25),3);
    mean_current_mouse_sessions(:,:,6)  =  nanmean(current_mouse(:,:,26:30),3);
    
    mean_pv_corr_pairwise_strict = [];
    for sessionA = 1:6
        sessionA_activity = mean_current_mouse_sessions(:,:,sessionA);
        for sessionB = 1:6
            sessionB_activity = mean_current_mouse_sessions(:,:,sessionB);
            
            valid_cells = [nanmean(sessionA_activity,2) ~= 0] & [nanmean(sessionB_activity,2) ~= 0];
            valid_sessionA_activity = sessionA_activity(valid_cells,:);
            valid_sessionB_activity = sessionB_activity(valid_cells,:);
            
            mean_pv_corr_pairwise_strict(sessionA,sessionB) = nanmean(diag(corr(valid_sessionA_activity,valid_sessionB_activity)));
            
        end
    end
    
    allen_elapsed_session(mouse,1) = nanmean([mean_pv_corr_pairwise_strict(1,2),mean_pv_corr_pairwise_strict(3,4),mean_pv_corr_pairwise_strict(5,6)]);
    allen_elapsed_session(mouse,2) = nanmean([nanmean(nanmean(mean_pv_corr_pairwise_strict(1:2,3:4))),nanmean(nanmean(mean_pv_corr_pairwise_strict(3:4,5:6)))]);
    allen_elapsed_session(mouse,3) = nanmean(nanmean(mean_pv_corr_pairwise_strict(1:2,5:6)));
    
    mean_current_mouse_sessions = [];
    mean_current_mouse_sessions(:,:,1) =  nanmean(current_mouse(:,:,1:10),3);
    mean_current_mouse_sessions(:,:,2)  =  nanmean(current_mouse(:,:,11:20),3);
    mean_current_mouse_sessions(:,:,3)  =  nanmean(current_mouse(:,:,21:30),3);
    num_cells_across_registrations(mouse,1) = size(mean_current_mouse_sessions,1);
    
    
    mean_pv_corr_pairwise_strict = [];
    all_active_cells_across_sess = [];
    active_both_cells_across_sess = [];
    for sessionA = 1:3
        sessionA_activity = mean_current_mouse_sessions(:,:,sessionA);
        for sessionB = 1:3
            sessionB_activity = mean_current_mouse_sessions(:,:,sessionB);
            
            valid_cells = [nanmean(sessionA_activity,2) ~= 0] & [nanmean(sessionB_activity,2) ~= 0];
            all_active_cells_across_sess(sessionA,sessionB) = sum([nanmean(sessionA_activity,2) ~= 0] | [nanmean(sessionB_activity,2) ~= 0]);
            active_both_cells_across_sess(sessionA,sessionB) = sum([nanmean(sessionA_activity,2) ~= 0] & [nanmean(sessionB_activity,2) ~= 0]);
            
            valid_sessionA_activity = sessionA_activity(valid_cells,:);
            valid_sessionB_activity = sessionB_activity(valid_cells,:);
            
            mean_pv_corr_pairwise_strict(sessionA,sessionB) = nanmean(diag(corr(valid_sessionA_activity,valid_sessionB_activity)));
        end
    end
    
    all_active_cells_across_mice(mouse,:) = [all_active_cells_across_sess(1,2),all_active_cells_across_sess(1,3),all_active_cells_across_sess(2,3)];
    active_both_cells_across_mice(mouse,:) = [active_both_cells_across_sess(1,2),active_both_cells_across_sess(1,3),active_both_cells_across_sess(2,3)];
    allen_pv_elapsed_full_sess(mouse,:) = [mean_pv_corr_pairwise_strict(1,2),mean_pv_corr_pairwise_strict(1,3),mean_pv_corr_pairwise_strict(2,3)];
    
    
    reg_file_name = dir([registration_path,mat_list{mouse}(1:end-4),'\registration\cellRegistered*.mat']);
    quality_file_name = dir([registration_path,mat_list{mouse}(1:end-4),'\registration\','*.txt']);
    
    load([results_path,mat_list{mouse}],'mouse_age','filtered_traces_days_events','cell_registration','sorted_elapsed_time')
    current_mouse = filtered_traces_days_events(:,1);
    [~,exp_order] = sort(mouse_age);
    sorted_elapsed_time_all_mice(mouse,:) = sorted_elapsed_time;
    load([registration_path,mat_list{mouse}(1:end-4),'\registration\aligned_data_struct.mat'])
    for psame = 1:3
        load([registration_path,mat_list{mouse}(1:end-4),'\registration\',reg_file_name(psame).name])
        
        quality_data = readtable([registration_path,mat_list{mouse}(1:end-4),'\registration\',quality_file_name(psame).name]);
        quality_across_mice(mouse,:,psame) = cellfun(@str2num,cellfun(@(x) x(31:end),table2cell(quality_data(30:31,1)),'UniformOutput',false)); % flase positive ; false negative
        
        cellreg = cell_registered_struct.cell_to_index_map;
        sorted_traces_across_days = zeros(size(cellreg,1),9000,3);
        cell_reg_allen_reg = [];
        sorted_centeroids_across_days = nan(size(cellreg,1),2,3);
        temp = aligned_data_struct.spatial_footprints_corrected{1};
        sorted_spatial_footprints_across_days =  nan(size(temp,2)*size(temp,3),size(cellreg,1),3);
        for session = 1:3
            experiment = exp_order(session);
            exp_list = dir([registration_path,mat_list{mouse}(1:end-4)]);
            exp_list(ismember({exp_list.name}, {'.', '..','registration'})) = [];
            myDir = find(vertcat(exp_list.isdir));
            exp_list = {exp_list(myDir).name};
            
            ROI_file_names  = dir([registration_path,mat_list{mouse}(1:end-4),'\',exp_list{experiment},'\ROIs\*.csv']);
            ROI_file_names = {ROI_file_names.name};
            currented_labels_sess = cell2mat(cellfun(@str2num,cellfun(@(x) x(1:end-4),ROI_file_names,'UniformOutput',false),'UniformOutput',false)');
            [~,ai,bi] = intersect(currented_labels_sess,cell_registration{experiment});
            current_sess = current_mouse{experiment}(bi,:);
            current_sess_foorprints = aligned_data_struct.spatial_footprints_corrected{1,experiment};
            
            sorted_reg = cell_registration{experiment}(bi);
            for cell = 1:size(current_sess,1)
                ordered_cell = find(cellreg(:,experiment) == cell);
                if ~isempty(ordered_cell)
                    sorted_traces_across_days(ordered_cell,:,session) = current_sess(cell,:);
                    current_foot_print = squeeze(current_sess_foorprints(cell,:,:));
                    sorted_spatial_footprints_across_days(:,ordered_cell,session) = current_foot_print(:);
                    cell_reg_allen_reg(ordered_cell) = sorted_reg(cell);
                end
            end
        end
        
        frames_rate = 30;
        repeats_movie1 = 10;
        movie1_length = 30;
        movie1_frames = frames_rate * movie1_length;
        movie1_bin_size = 30;
        binned_movie1 = ones(movie1_frames,1);
        movie1_bin_edges = 1:movie1_bin_size:movie1_frames;
        for bin = 1:length(movie1_bin_edges)
            binned_movie1(movie1_bin_edges(bin):movie1_bin_edges(bin)+movie1_bin_size-1) = bin;
        end
        binned_movie_repeated1 = repmat(binned_movie1,[repeats_movie1,1]);
        
        % binization for natural movies
        pop_vector_info_trials = [];
        sub  = 1;
        for session = 1:3
            for repeat = 1:10
                frames_temp = [1:movie1_frames] + (movie1_frames*(repeat-1));
                current_repeat = sorted_traces_across_days(:,frames_temp,session);
                
                for bin = 1:length(movie1_bin_edges)
                    pop_vector_info_trials(:,bin,sub) = nanmean(current_repeat(:,binned_movie_repeated1(frames_temp) == bin),2);
                    
                end
                sub = sub + 1;
            end
        end
        
        
        mean_current_mouse_sessions = [];
        mean_current_mouse_sessions(:,:,1) =  nanmean(pop_vector_info_trials(:,:,1:5),3);
        mean_current_mouse_sessions(:,:,2)  =  nanmean(pop_vector_info_trials(:,:,6:10),3);
        mean_current_mouse_sessions(:,:,3)  =  nanmean(pop_vector_info_trials(:,:,11:15),3);
        mean_current_mouse_sessions(:,:,4) =  nanmean(pop_vector_info_trials(:,:,16:20),3);
        mean_current_mouse_sessions(:,:,5)  =  nanmean(pop_vector_info_trials(:,:,21:25),3);
        mean_current_mouse_sessions(:,:,6)  =  nanmean(pop_vector_info_trials(:,:,26:30),3);
        
        mean_pv_corr_pairwise_strict = [];
        for sessionA = 1:6
            sessionA_activity = mean_current_mouse_sessions(:,:,sessionA);
            for sessionB = 1:6
                sessionB_activity = mean_current_mouse_sessions(:,:,sessionB);
                
                valid_cells = [nanmean(sessionA_activity,2) ~= 0] & [nanmean(sessionB_activity,2) ~= 0];
                valid_sessionA_activity = sessionA_activity(valid_cells,:);
                valid_sessionB_activity = sessionB_activity(valid_cells,:);
                
                mean_pv_corr_pairwise_strict(sessionA,sessionB) = nanmean(diag(corr(valid_sessionA_activity,valid_sessionB_activity)));
            end
        end
        
        cellreg_elapsed_session(mouse,1,psame) = nanmean([mean_pv_corr_pairwise_strict(1,2),mean_pv_corr_pairwise_strict(3,4),mean_pv_corr_pairwise_strict(5,6)]);
        cellreg_elapsed_session(mouse,2,psame) = nanmean([nanmean(nanmean(mean_pv_corr_pairwise_strict(1:2,3:4))),nanmean(nanmean(mean_pv_corr_pairwise_strict(3:4,5:6)))]);
        cellreg_elapsed_session(mouse,3,psame) = nanmean(nanmean(mean_pv_corr_pairwise_strict(1:2,5:6)));
        
        mean_current_mouse_sessions = [];
        mean_current_mouse_sessions(:,:,1) =  nanmean(pop_vector_info_trials(:,:,1:10),3);
        mean_current_mouse_sessions(:,:,2)  =  nanmean(pop_vector_info_trials(:,:,11:20),3);
        mean_current_mouse_sessions(:,:,3)  =  nanmean(pop_vector_info_trials(:,:,21:30),3);
        num_cells_across_registrations(mouse,psame+1) = size(mean_current_mouse_sessions,1);
        
        mean_pv_corr_pairwise_strict = [];
        
        active_both_cells_across_sess = [];
        all_active_cells_across_sess = [];
        
        pairs_footprint = {};
        pairs_footprint_diff = {};
        for sessionA = 1:3
            sessionA_activity = mean_current_mouse_sessions(:,:,sessionA);
            for sessionB = 1:3
                sessionB_activity = mean_current_mouse_sessions(:,:,sessionB);
                
                valid_cells = [nanmean(sessionA_activity,2) ~= 0] & [nanmean(sessionB_activity,2) ~= 0];
                valid_sessionA_activity = sessionA_activity(valid_cells,:);
                valid_sessionB_activity = sessionB_activity(valid_cells,:);
                
                all_active_cells_across_sess(sessionA,sessionB) = sum([nanmean(sessionA_activity,2) ~= 0] | [nanmean(sessionB_activity,2) ~= 0]);
                active_both_cells_across_sess(sessionA,sessionB) = sum([nanmean(sessionA_activity,2) ~= 0] & [nanmean(sessionB_activity,2) ~= 0]);
                mean_pv_corr_pairwise_strict(sessionA,sessionB) = nanmean(diag(corr(valid_sessionA_activity,valid_sessionB_activity)));
                
                if sessionA < sessionB
                    footprint_corr_per_sess = corr(sorted_spatial_footprints_across_days(:,:,sessionA),sorted_spatial_footprints_across_days(:,:,sessionB));
                    pairs_footprint{sessionA,sessionB} = diag(footprint_corr_per_sess);
                    temp = footprint_corr_per_sess;
                    temp(boolean(eye(size(temp,1)))) = NaN;
                    pairs_footprint_diff{sessionA,sessionB} = max(temp,[],2,'omitnan');
                end
            end
        end
        
        if psame == 2
            cellreg_all_active_cells_across_mice(mouse,:) = [all_active_cells_across_sess(1,2),all_active_cells_across_sess(1,3),all_active_cells_across_sess(2,3)];
            cellreg_active_both_cells_across_mice(mouse,:) = [active_both_cells_across_sess(1,2),active_both_cells_across_sess(1,3),active_both_cells_across_sess(2,3)];
            
            cellreg_pv_elapsed_full_sess(mouse,:) = [mean_pv_corr_pairwise_strict(1,2),mean_pv_corr_pairwise_strict(1,3),mean_pv_corr_pairwise_strict(2,3)];
            
            pairs_footprint_same_cells_across_mice(mouse) = {[pairs_footprint{1,2},pairs_footprint{2,3},pairs_footprint{1,3}]};
            pairs_footprint_diff_cells_across_mice(mouse) = {[pairs_footprint_diff{1,2},pairs_footprint_diff{2,3},pairs_footprint_diff{1,3}]};
            sorted_elapsed_time_per_cell_across_mice(mouse) = {repmat(sorted_elapsed_time([1,3,2]),[length(pairs_footprint{1,2}),1])};
        end
    end
end


%% Figure S6E - Percentage of registration errors

figure('units','normalized','position',[0.35 0.35 0.2 0.3])
figure_boxplot(quality_across_mice(:,:,2)*100)
set(gca,'xtick',[1,2],'xticklabel',{'False positives','False negatives'})
ylabel('Percentage of registration errors')
ylim([-1 30])

%% Figure S6F - Total number of unique cells registered across all sessions as function of Psame

figure('units','normalized','position',[0.35 0.35 0.2 0.3])
figure_boxplot(num_cells_across_registrations)
ylabel({'Total number of cells';'registered across all sessions'})
set(gca,'xtick',1:4,'xticklabel',{'ABO default','Psame = 0.05','Psame = 0.50','Psame = 0.95'})
xtickangle(15)

%% Figure S6G - Fraction of active cells using Sheituch or ABO registration

overlap_allen = active_both_cells_across_mice./all_active_cells_across_mice;
overlap_cellreg = cellreg_active_both_cells_across_mice./cellreg_all_active_cells_across_mice;

temp_allen = overlap_allen(:);
temp_cellreg = overlap_cellreg(:);
[~,i] = sort(temp_allen);
x = temp_allen(i);
y = temp_cellreg(i);
[p,s] = polyfit(x,y,1);
[yfit,dy] = polyconf(p,x,s,'predopt','curve');

[r,p] = corr(temp_allen,temp_cellreg);

figure('units','normalized','position',[0.35 0.35 0.2 0.3])
hold on
scatter(temp_allen,temp_cellreg,[],[0.8 0.8 0.8],'filled','MarkerfaceAlpha',0.6)
h = fill([x;flipud(x)],[yfit-dy;flipud(yfit+dy)],newmap3(1,:),'linestyle','none');
line(x,yfit,'color',newmap3(1,:),'linewidth',2)
plot([0 1],[0 1],'--','color',[0.5 0.5 0.5],'linewidth',0.5)

text(0.05, 0.95,['R^2 = ',num2str(round(r.^2,4))],'Units','normalized','fontsize',10)
text(0.05, 0.885,['p = ',num2str(p)],'Units','normalized','fontsize',10)
text(0.05, 0.82,['N = ',num2str(size(x,1)),' pairs of sessions'],'Units','normalized','fontsize',10)

set(h,'facealpha',.25)
ylabel({'Fraction of active cells using';'Sheituch et al. (2017) registration'})
xlabel({'Fraction of active cells using';'ABO default registration'})
ylim([0 1])
xlim([0 1])



%% Figure S6H - PV correlation across sessions using Sheituch or ABO registration

temp_allen = allen_pv_elapsed_full_sess(:);
temp_cellreg = cellreg_pv_elapsed_full_sess(:);
[~,i] = sort(temp_allen);
x = temp_allen(i);
y = temp_cellreg(i);
[p,s] = polyfit(x,y,1);
[yfit,dy] = polyconf(p,x,s,'predopt','curve');

[r,p] = corr(temp_allen,temp_cellreg);

figure('units','normalized','position',[0.35 0.35 0.2 0.3])
hold on
scatter(temp_allen,temp_cellreg,[],[0.8 0.8 0.8],'filled','MarkerfaceAlpha',0.6)
h = fill([x;flipud(x)],[yfit-dy;flipud(yfit+dy)],newmap3(1,:),'linestyle','none');
line(x,yfit,'color',newmap3(1,:),'linewidth',2)
plot([0 1],[0 1],'--','color',[0.5 0.5 0.5],'linewidth',0.5)

text(0.05, 0.95,['R^2 = ',num2str(round(r.^2,4))],'Units','normalized','fontsize',10)
text(0.05, 0.885,['p = ',num2str(p)],'Units','normalized','fontsize',10)
text(0.05, 0.82,['N = ',num2str(size(x,1)),' pairs of sessions'],'Units','normalized','fontsize',10)

set(h,'facealpha',.25)
xlabel({'PV correlation using';'Sheituch et al. (2017) registration'})
ylabel({'PV correlation using';'ABO default registration'})
ylim([0 1])
xlim([0 1])



%% Figure S6I - PV correlation across session as a function of Psame

figure('units','normalized','position',[0.35 0.35 0.225 0.35])
mean_stability = nanmean(allen_elapsed_session);
std_stability = nanstd(allen_elapsed_session)./sqrt(size(allen_elapsed_session,1));
errorbar([0.5 3.5 6.5],mean_stability,std_stability,'o','color',colors(1,:),...
    'markerfacecolor',colors(1,:),'capsize',0,'linestyle','none','linewidth',2)
hold on
for psame = 1:3
    mean_stability = nanmean(cellreg_elapsed_session(:,:,psame));
    std_stability = nanstd(cellreg_elapsed_session(:,:,psame))./sqrt(size(cellreg_elapsed_session,1));
    errorbar([0.5 3.5 6.5] + 0.25*(psame),mean_stability,std_stability,'o','color',[0.8 0.8 0.8]-0.3*(psame-1),...
        'markerfacecolor',[0.8 0.8 0.8]-0.3*(psame-1),'capsize',0,'linestyle','none','linewidth',2)
end
legend({'ABO default','Psame = 0.05','Psame = 0.50','Psame = 0.95'})
legend('boxoff')
xlim([0 8])
set(gca,'xtick',0.85:3:6.85,'xticklabel',{'Within session','Proximal sessions','Distal sessions'},'box','off')
ylabel('PV correlation')

%% Figure S6J - Spatial footprint correlation within and across cells as a function of elapsed time

temp1 = cell2mat(pairs_footprint_same_cells_across_mice');
temp1 = temp1(:);
temp2 = cell2mat(pairs_footprint_diff_cells_across_mice');
temp2 = temp2(:);
temp3 = cell2mat(sorted_elapsed_time_per_cell_across_mice');
temp3 = temp3(:);
day_vals = unique(temp3);
mean_similarity = [];
std_similarity = [];
for day = 1:length(day_vals)
    mean_similarity(1,day) = nanmean(temp1(temp3==day));
    mean_similarity(2,day) =  nanmean(temp2(temp3==day));
    std_similarity(1,day) = nanstd(temp1(temp3==day));
    std_similarity(2,day) = nanstd(temp2(temp3==day));
end

figure('units','normalized','position',[0.3 0.3 0.225 0.35])
hold on
xlim([0 21])
ylim([-0.2 1.19])
errorbar(mean_similarity(1,:),std_similarity(1,:),'o','color',colors(1,:),'linewidth',2,'capsize',0,'MarkerFaceColor',colors(1,:))
errorbar(mean_similarity(2,:),std_similarity(2,:),'o','color',[0.7 0.7 0.7],'linewidth',2,'capsize',0,'MarkerFaceColor',[0.7 0.7 0.7])
text(0.7, 0.95,['Within cells'],'Units','normalized','color',colors(1,:),'fontsize',12)
text(0.7, 0.45,['Across cells'],'Units','normalized','color',[0.7 0.7 0.7],'fontsize',12)
text(0.685, 0.385,['(most similar)'],'Units','normalized','color',[0.7 0.7 0.7],'fontsize',12)

ylabel('Spatial footprint correlation')
xlabel('Days between recording sessions')

%% Figure S7A - Internal structure in the reduced space using tSNE as a function of cell count - takes ~20 min to run

% perform tsne analysis for neuropixels data
nat_movie = 1; % natural movie 1
area = 1; % area V1

valid_mice = movie_repeats(:,nat_movie) == 30; % include only mice from the neuropixels 'Functional connectivity' group
V1_pseudo = cell2mat(neuropixels_population_vectors_tsne(valid_mice,area,nat_movie)); % pool all V1 cells across mice to create a single pseudo mouse
V1_pseudo = V1_pseudo(:,:,31:60); % subset the neuronal activity during block B (repeats 31-60)


load('FigureS7A_neuropixels.mat','state') % load seed for tSNE visualization
rng(state) % set seed
cells_rand_ind = randperm(size(V1_pseudo,1)); % random permutatuion on cells IDs that will later be subsampled from the entire population
cell_num_list = [75,100,125,150,200,250,500,1000,1500]; % list of number of cells included in the analysis

neuropixels_dim_reduceV1_cell_count = {}; % will store the reduced space comp for each tsne analysis
for cell_included = 1:length(cell_num_list) % loop over cell count treshold
    clc;
    disp('Calculating structures:')
    disp(['Dataset: Neuropixels | Cells: ',num2str(cell_num_list(cell_included)), ' | ',...
        num2str(cell_included),'/',num2str(length(cell_num_list)),...
        ' | ', num2str([cell_included./length(cell_num_list)]*100),'%'])
    
    cell_num = cells_rand_ind(1:cell_num_list(cell_included)); % subsampled cells ID
    subset_V1_pseudo = V1_pseudo(cell_num,:,:); % subsample neuronal activity
    V1_pseudo_strcut = reshape(subset_V1_pseudo,[size(subset_V1_pseudo,1),90*30]); % reshape the neuronal activity from 3D into 2D (#units by 90 time bins x 30 movie repeats)
   
    dim_reduceV1 = tsne(V1_pseudo_strcut','Algorithm','exact','Distance','cosine',...
        'NumDimensions',2,'NumPCAComponents',20,'Perplexity',200); % perform tSNE dim reduction on cells
    neuropixels_dim_reduceV1_cell_count{cell_included} = dim_reduceV1; % store tsne components
end

% visualize the internal structure as function of number of cells included in the analysis
xy_list = [2,1;1 2;1 2;2 1;1 2;1 2;1 2;1 2;1 2];
direction_list = [1 1;1 1;-1 1;1 -1; -1 -1;1 -1;-1 -1;1 -1; -1 -1];
figure('units','normalized','position',[0.3 0.2 0.3 0.525]) % neuropixels
for cell_included = 1:length(cell_num_list)
    current_dir = direction_list(cell_included,:);
    current_xy = xy_list(cell_included,:);
    dim_reduceV1 = neuropixels_dim_reduceV1_cell_count{cell_included};
    subplot(3,3,cell_included)
    scatter(current_dir(1).*dim_reduceV1(:,current_xy(1)),current_dir(2).*dim_reduceV1(:,current_xy(2)),15,repmat([1:90],[1 30]),'filled')
    xlim([-25 25])
    ylim([-25 25])
    text(0.65,0.9,['N=',num2str(cell_num_list(cell_included))],'Units','normalized','FontSize',10)
    set(gca,'xtick',[],'ytick',[])
    if cell_included == 7
        ylabel('Comp. 2')
        xlabel('Comp. 1')
    end
end
suptitle('tSNE | Neuropixels | Natural movie 1 | Area V1:')
colormap(new_jet_colormap)
cb = colorbar;
set(cb,'position',[0.925 0.105 0.04 0.77],'Ticks',[])
cb.Label.Position(1) = 0.005;
cb.Label.String = 'Start                     Time in movie                     End';
cb.FontSize = 14;

% perform tsne analysis for calcium imaging data
area = 1; % area V1
V1_pseudo = cell2mat(calcium_excitatory_population_vectors{area}(:,3)); % pool all V1 cells across mice to create a single pseudo mouse

active_sess1 = mean(mean(V1_pseudo(:,:,1:10),3,'omitnan'),2,'omitnan')>0; % active cells in session 1 (repeats 1-10)
active_sess2 = mean(mean(V1_pseudo(:,:,11:20),3,'omitnan'),2,'omitnan')>0; % active cells in session 1 (repeats 11-20)
active_sess3 = mean(mean(V1_pseudo(:,:,21:30),3,'omitnan'),2,'omitnan')>0; % active cells in session 1 (repeats 21-30)
valid_cells = find([active_sess1 & active_sess2 & active_sess3]); % active in all three sessions

V1_pseudo = V1_pseudo(valid_cells,:,:); % include only cells that were active in all three sessions
V1_pseudo(V1_pseudo==0) = 0.000001; % convert zero values into epsilon due to sparssness of data (required for tsne)

load('FigureS7A_calcium.mat','state') % load seed for tSNE visualization
rng(state) % set seed

cells_rand_ind = randperm(size(V1_pseudo,1)); % random permutatuion on cells IDs that will later be subsampled from the entire population
cell_num_list = [75,100,125,150,200,250,500,1000,1500]; % list of number of cells included in the analysis


calcium_dim_reduceV1_cell_count = {}; % will store the reduced space comp for each tsne analysis
for cell_included = 1:length(cell_num_list)  % loop over cell count treshold
    clc;
    disp('Calculating structures:')
    disp(['Dataset: Calcium imaging | Cells: ',num2str(cell_num_list(cell_included)), ' | ',...
        num2str(cell_included),'/',num2str(length(cell_num_list)),...
        ' | ', num2str([cell_included./length(cell_num_list)]*100),'%'])
    
    cell_num = cells_rand_ind(1:cell_num_list(cell_included)); % subsampled cells ID
    subset_V1_pseudo = V1_pseudo(cell_num,:,:); % subsample neuronal activity
    V1_pseudo_strcut = reshape(subset_V1_pseudo,[size(subset_V1_pseudo,1),90*30]); % reshape the neuronal activity from 3D into 2D (#cells by 90 time bins x 30 movie repeats)

    dim_reduceV1 = tsne(V1_pseudo_strcut','Algorithm','exact','Distance','cosine',...
        'NumDimensions',2,'NumPCAComponents',20,'Perplexity',200); % perform tSNE dim reduction on cells
    calcium_dim_reduceV1_cell_count{cell_included} = dim_reduceV1; % store tsne components
end

% visualize the internal structure as function of number of cells included in the analysis
xy_list = [2 1;1 2;1 2;2 1;2 1;2 1;2 1;1 2;1 2];
direction_list = [-1 -1;1 -1;-1 1;1 -1; 1 -1;-1 1;1 1;1 -1; 1 1];
figure('units','normalized','position',[0.3 0.2 0.3 0.525]) % calcium imaging
for cell_included = 1:length(cell_num_list)
    current_dir = direction_list(cell_included,:);
    current_xy = xy_list(cell_included,:);
    dim_reduceV1 = calcium_dim_reduceV1_cell_count{cell_included};
    subplot(3,3,cell_included)
    scatter(current_dir(1).*dim_reduceV1(:,current_xy(1)),current_dir(2).*dim_reduceV1(:,current_xy(2)),15,repmat([1:90],[1 30]),'filled')
    xlim([-20 20])
    ylim([-20 20])
    text(0.65,0.9,['N=',num2str(cell_num_list(cell_included))],'Units','normalized','FontSize',10)
    set(gca,'xtick',[],'ytick',[])
    if cell_included == 7
        ylabel('Comp. 2')
        xlabel('Comp. 1')
    end
end
suptitle('tSNE | Calcium imaging | Natural movie 1 | Area V1:')
colormap(new_jet_colormap)
cb = colorbar;
set(cb,'position',[0.925 0.105 0.04 0.77],'Ticks',[])
cb.Label.Position(1) = 0.005;
cb.Label.String = 'Start                     Time in movie                     End';
cb.FontSize = 14;


%% Figure S7B - Internal structure in the reduced space using PCA as a function of cell count

% perform tsne analysis for neuropixels data
nat_movie = 1; % natural movie 1
area = 1; % area V1

valid_mice = movie_repeats(:,nat_movie) == 30; % include only mice from the neuropixels 'Functional connectivity' group
V1_pseudo = cell2mat(neuropixels_population_vectors_tsne(valid_mice,area,nat_movie)); % pool all V1 cells across mice to create a single pseudo mouse
V1_pseudo = V1_pseudo(:,:,31:60); % subset the neuronal activity during block B (repeats 31-60)


load('FigureS7A_neuropixels.mat','state')% load seed for PCA visualization
rng(state) % set seed
cells_rand_ind = randperm(size(V1_pseudo,1)); % random permutatuion on cells IDs that will later be subsampled from the entire population
cell_num_list = [75,100,125,150,200,250,500,1000,1500]; % list of number of cells included in the analysis

neuropixels_dim_reduceV1_cell_count = {}; % will store the reduced space comp for each tsne analysis
for cell_included = 1:length(cell_num_list)
   clc;
    disp('Calculating structures:')
    disp(['Dataset: Neuropixels | Cells: ',num2str(cell_num_list(cell_included)), ' | ',...
        num2str(cell_included),'/',num2str(length(cell_num_list)),...
        ' | ', num2str([cell_included./length(cell_num_list)]*100),'%'])
    
    cell_num = cells_rand_ind(1:cell_num_list(cell_included)); % subsampled cells ID
    subset_V1_pseudo = V1_pseudo(cell_num,:,:); % subsample neuronal activity
    V1_pseudo_strcut = reshape(subset_V1_pseudo,[size(subset_V1_pseudo,1),90*30]); % reshape the neuronal activity from 3D into 2D (#units by 90 time bins x 30 movie repeats)
    
    [~,dim_reduceV1] = pca(V1_pseudo_strcut'); % perform PCA dim reduction on cells
    neuropixels_dim_reduceV1_cell_count{cell_included} = dim_reduceV1; % store tsne components
end

% visualize the internal structure as function of number of cells included in the analysis
xy_list = [1,2;1 2;1 2;2 1;1 2;1 2;1 2;1 2;1 2];
direction_list = [1 1;1 1;-1 1;1 -1; -1 -1;1 -1;-1 -1;1 -1; -1 -1];
figure('units','normalized','position',[0.3 0.2 0.3 0.525]) % neuropixels
for cell_included = 1:length(cell_num_list)
    current_dir = direction_list(cell_included,:);
    current_xy = xy_list(cell_included,:);
    dim_reduceV1 = neuropixels_dim_reduceV1_cell_count{cell_included};
    subplot(3,3,cell_included)
    scatter(dim_reduceV1(:,1),dim_reduceV1(:,2),15,repmat([1:90],[1 30]),'filled')
    
    xlim([-10 10])
    ylim([-10 10])
    text(0.65,0.9,['N=',num2str(cell_num_list(cell_included))],'Units','normalized','FontSize',10)
    set(gca,'xtick',[],'ytick',[])
    if cell_included == 7
        ylabel('PC 2')
        xlabel('Comp. 1')
    end
end
suptitle('PCA | Neuropixels | Natural movie 1 | Area V1:')
colormap(new_jet_colormap)
cb = colorbar;
set(cb,'position',[0.925 0.105 0.04 0.77],'Ticks',[])
cb.Label.Position(1) = 0.005;
cb.Label.String = 'Start                     Time in movie                     End';
cb.FontSize = 14;


% perform tsne analysis for calcium imaging data
area = 1; % area V1
V1_pseudo = cell2mat(calcium_excitatory_population_vectors{area}(:,3));

active_sess1 = mean(mean(V1_pseudo(:,:,1:10),3,'omitnan'),2,'omitnan')>0; % active cells in session 1 (repeats 1-10)
active_sess2 = mean(mean(V1_pseudo(:,:,11:20),3,'omitnan'),2,'omitnan')>0; % active cells in session 1 (repeats 11-20)
active_sess3 = mean(mean(V1_pseudo(:,:,21:30),3,'omitnan'),2,'omitnan')>0; % active cells in session 1 (repeats 21-30)
valid_cells = find([active_sess1 & active_sess2 & active_sess3]); % active in all three sessions
V1_pseudo = V1_pseudo(valid_cells,:,:); % include only cells that were active in all three sessions
V1_pseudo(V1_pseudo==0) = 0.000001; % convert zero values into epsilon due to sparssness of data (required for tsne)

load('FigureS7A_calcium.mat','state') % load seed for PCA visualization
rng(state) % set seed

cells_rand_ind = randperm(size(V1_pseudo,1)); % random permutatuion on cells IDs that will later be subsampled from the entire population
cell_num_list = [75,100,125,150,200,250,500,1000,1500]; % list of number of cells included in the analysis

calcium_dim_reduceV1_cell_count = {}; % will store the reduced space comp for each tsne analysis
for cell_included = 1:length(cell_num_list)  % loop over cell count treshold
    clc;
    disp('Calculating structures:')
    disp(['Dataset: Calcium imaging | Cells: ',num2str(cell_num_list(cell_included)), ' | ',...
        num2str(cell_included),'/',num2str(length(cell_num_list)),...
        ' | ', num2str([cell_included./length(cell_num_list)]*100),'%'])
    
    cell_num = cells_rand_ind(1:cell_num_list(cell_included));% subsampled cells ID
    subset_V1_pseudo = V1_pseudo(cell_num,:,:); % subsample neuronal activity
    V1_pseudo_strcut = reshape(subset_V1_pseudo,[size(subset_V1_pseudo,1),90*30]); % reshape the neuronal activity from 3D into 2D (#celss by 90 time bins x 30 movie repeats)
    
    [~,dim_reduceV1] = pca(V1_pseudo_strcut'); % perform PCA dim reduction on cells
    calcium_dim_reduceV1_cell_count{cell_included} = dim_reduceV1; % store tsne components
end

% visualize the internal structure as function of number of cells included in the analysis
figure('units','normalized','position',[0.3 0.2 0.3 0.525])
for cell_included = 1:length(cell_num_list)
    current_dir = direction_list(cell_included,:);
    current_xy = xy_list(cell_included,:);
    dim_reduceV1 = calcium_dim_reduceV1_cell_count{cell_included};
    subplot(3,3,cell_included)
    scatter(dim_reduceV1(:,1),dim_reduceV1(:,2),15,repmat([1:90],[1 30]),'filled')
    
    xlim([-2 2])
    ylim([-2 2])
    text(0.65,0.9,['N=',num2str(cell_num_list(cell_included))],'Units','normalized','FontSize',10)
    set(gca,'xtick',[],'ytick',[])
    if cell_included == 7
        ylabel('PC 2')
        xlabel('Comp. 1')
    end
end
suptitle('PCA | Calcium imaging | Natural movie 1 | Area V1:')
colormap(new_jet_colormap)
cb = colorbar;
set(cb,'position',[0.925 0.105 0.04 0.77],'Ticks',[])
cb.Label.Position(1) = 0.005;
cb.Label.String = 'Start                     Time in movie                     End';
cb.FontSize = 14;

%% Figure S7C - Between pseudo-mice decoder - cell count (calcium)- takes hours!
nat_movie = 1; % natural movie 1

% subsample for each mouse in the neuropixels dataset the neuronal activity during the first 20 movie repeats
subset_population_vectors = {}; % define an empty variable that will store the subsmaple neuronal activity of indevidual mice
for area = 1:6 % loop over areas
    subset_population_vectors{area} = calcium_excitatory_population_vectors{area}(:,nat_movie); 
end

% decode the identity of each visual area based on the
% similarity between the internal structures of the two psedo-mice
cell_cutoff_list = [5,10,25,50,75,100,150,200,250,300,350,400,500,600,700,800,1000,1200,1400]; % define list of how many cells are included in the analysis
between_similarity_acc_all_cutoffs = []; % define an empty variable that will store the decoder results for the non-shuffled pseudo-mice
between_similarity_acc_shuffled_all_cutoffs = []; % define an empty variable that will store the decoder results for the suffled pseudo-mice

num_shuffles = 3000;% number of pseudo-mice realizations for each cell count threshold
min_cell_num_all_shuffles = []; % define empty variable to store number of cells sampled in leach realization of pseudo-mice
for cells_included = 1:length(cell_cutoff_list)+1 % loop over cell count thresholds
    
    between_similarity_acc = nan(num_shuffles,6); % define empty variable that store decoders predictions for non-shuffled pseudo-mice
    between_similarity_acc_shuffled = nan(num_shuffles,6); % define empty variable that store decoders predictions for shuffled pseudo-mice
    
    for shuffle = 1:num_shuffles % loop over realizations
    
        % determine the number of cells to be incuded in the analysis
        if cells_included <= length(cell_cutoff_list) % if below 1400 cells than use 'cell_cutoff_list'
            breaker = cell_cutoff_list(cells_included);
        else %used the entire dataset
            breaker = 1;
        end
        
        % create realizations of pseudo-mice until all pseudo-areas contain
        % more cells relative to 'cell_cutoff_list'
        pseudo_area_cell_num = zeros(2,6); % define a variable for of zeros that will store the number of units in each pseudo-area
        
        while   breaker > min(pseudo_area_cell_num(:)) % check if all pseudo-areas contain more cells relative to 'cell_cutoff_list'
            
            % for each visual area in each pseudo mouse, pool all units across mice
            % to create 12 pseudo areas (6 areas x 2 pseudo-mice)
            pseudo_mouseA = {}; % define an empty variable that will store the pooled units across mice for each pseudo-area of pseudo-mouse A
            pseudo_mouseB = {}; % define an empty variable that will store the pooled units across mice for each pseudo-area of pseudo-mouse B
            for area = 1:6
                % split the dataset into two independent group of mice
                pseudo_mouseA_ind = sort(randperm(size(subset_population_vectors{area},1),round(size(subset_population_vectors{area},1)./2)));
                pseudo_mouseB_ind = find(~ismember([1:size(subset_population_vectors{area},1)],pseudo_mouseA_ind));
                
                pseudo_mouseA{area} = cell2mat(subset_population_vectors{area}(pseudo_mouseA_ind)); % pool units across mice for a single visual area (pseudo-mouse A)
                pseudo_mouseB{area} = cell2mat(subset_population_vectors{area}(pseudo_mouseB_ind)); % pool units across mice for a single visual area (pseudo-mouse B)
                pseudo_area_cell_num(1,area) = size(pseudo_mouseA{area},1); % store number of cells in pseudo-area of pseudo-mouse A
                pseudo_area_cell_num(2,area) = size(pseudo_mouseB{area},1); % store number of cells in pseudo-area of pseudo-mouse B
            end
            % if the minimum value in "pseudo_area_cell_num" is higher than
            % the value in "breaker" than the loop breaks and the algorithm continues
        end
        
        % min_cell_num - the minimal number of cells across pseudo-areas and pseudo-mice.
        % this number will be used to randomly subsample the same number of cells for all pseudo-ares
        if cells_included <= length(cell_cutoff_list) % when below 1400
            min_cell_num = cell_cutoff_list(cells_included);
        else % when using the entire dataset
            min_cell_num = min(pseudo_area_cell_num(:));
            min_cell_num_all_shuffles(shuffle) = min_cell_num;
        end
        
         clc;
        disp(['Performing between calcium imaging pseudo-mice decoding. Cells included: ',num2str(min_cell_num),' | Realization: ',num2str(shuffle),'/',num2str(num_shuffles)])

        % subsampling randomly the same number of units for all pseudo-areas
        pseudo_mouseA_subset = {}; % define an empty variable that will store the subsampled units for each pseudo-area of pseudo-mouse A
        pseudo_mouseB_subset = {}; % define an empty variable that will store the subsampled units for each pseudo-area of pseudo-mouse B
        for area = 1:6 % loop over areas
            % random sampling #min_cell_num of units from each pseudo-area of pseudo-mouse A
            subset_cell_ids_mouseA = sort(randperm(pseudo_area_cell_num(1,area),min_cell_num));
            pseudo_mouseA_subset{area} = pseudo_mouseA{area}(subset_cell_ids_mouseA,:,:);
            
            % random sampling #min_cell_num of units from each pseudo-area of pseudo-mouse B
            subset_cell_ids_mouseB = sort(randperm(pseudo_area_cell_num(2,area),min_cell_num));
            pseudo_mouseB_subset{area} = pseudo_mouseB{area}(subset_cell_ids_mouseB,:,:);
        end
        
        % creating shuffled pseudo mice
        all_cells_pseudo_mouseA = cell2mat(pseudo_mouseA_subset'); % pooling all the cells across all pseudo-areas of pseudo-mouse A
        all_cells_pseudo_mouseB = cell2mat(pseudo_mouseB_subset'); % pooling all the cells across all pseudo-areas of pseudo-mouse B
        
        % random redistribution of cells across pseudo-areas within a given pseudo-mouse
        rand_cells_id_mouseA = randperm(size(all_cells_pseudo_mouseA,1)); % random permutation of cells indices for pseudo-mouse A
        rand_cells_id_mouseB = randperm(size(all_cells_pseudo_mouseB,1)); % random permutation of cells indices for pseudo-mouse B
        
        shuffled_pseudo_mouseA_subset = {}; % define an empty variable that will store the redistributed cells across areas for pseudo-mouse A
        shuffled_pseudo_mouseB_subset = {}; % define an empty variable that will store the redistributed cells across areas for pseudo-mouse B
        for area = 1:6 % loop over areas
            current_pseudo_area = [1:min_cell_num] + min_cell_num*(area-1);  % define range of cell indices for the current area
            
            shuffled_pseudo_mouseA_subset{area} = all_cells_pseudo_mouseA(rand_cells_id_mouseA(current_pseudo_area),:,:); % randomly subsample #min_cell_num cells to current pseud-area of pseudo-mouse A
            shuffled_pseudo_mouseB_subset{area} = all_cells_pseudo_mouseB(rand_cells_id_mouseB(current_pseudo_area),:,:); % randomly subsample #min_cell_num cells to current pseud-area of pseudo-mouse B
        end

        % calculating the internal structures for each movie repeat for all
        % pseudo-areas of both example pseudo-mice
        internal_structures_pseudoA = []; % define an empty variable the will store the internal structures of pseudo-mouse A
        internal_structures_pseudoB = []; % define an empty variable the will store the internal structures of pseudo-mouse B
        internal_structures_pseudoA_shuffle = []; % define an empty variable the will store the internal structures of shuffled pseudo-mouse A
        internal_structures_pseudoB_shuffle = []; % define an empty variable the will store the internal structures of shuffled pseudo-mouse B
        internal_structures_labels = []; % define an empty variable that will store the label of each internal structure
        triu_ind = boolean(triu(ones(30),1)); % define a boolian matrix with true values in the upper half of the matrix
        for area = 1:6 % loop over areas
            % for non-shuffled pseudo-mice
            current_structure_mouseA = [];
            current_structure_mouseA(:,:,1) = corr(mean(pseudo_mouseA_subset{area}(:,:,1:10),3,'omitnan')); % calculate the internal struture of session 1 for a single area in pseudo-mouse A
            current_structure_mouseA(:,:,2) = corr(mean(pseudo_mouseA_subset{area}(:,:,11:20),3,'omitnan')); % calculate the internal struture of session 2 for a single area in pseudo-mouse A
            current_structure_mouseA(:,:,3) = corr(mean(pseudo_mouseA_subset{area}(:,:,21:30),3,'omitnan')); % calculate the internal struture of session 3 for a single area in pseudo-mouse A
            current_structure_mouseA = mean(current_structure_mouseA,3,'omitnan'); % average internal structure across sessions
            internal_structures_pseudoA(:,area) = current_structure_mouseA(triu_ind);  % vectorize and store the internal structure
          
            current_structure_mouseB = [];
            current_structure_mouseB(:,:,1) = corr(mean(pseudo_mouseB_subset{area}(:,:,1:10),3,'omitnan')); % calculate the internal struture of session 1 for a single area in pseudo-mouse B
            current_structure_mouseB(:,:,2) = corr(mean(pseudo_mouseB_subset{area}(:,:,11:20),3,'omitnan')); % calculate the internal struture of session 2 for a single area in pseudo-mouse B
            current_structure_mouseB(:,:,3) = corr(mean(pseudo_mouseB_subset{area}(:,:,21:30),3,'omitnan')); % calculate the internal struture of session 3 for a single area in pseudo-mouse B
            current_structure_mouseB = mean(current_structure_mouseB,3,'omitnan'); % average internal structure across sessions
            internal_structures_pseudoB(:,area) = current_structure_mouseB(triu_ind); % vectorize and store the internal structure
           

            % for shuffled pseudo-mice
            current_structure_mouseA_shuffle = [];
            current_structure_mouseA_shuffle(:,:,1) = corr(mean(shuffled_pseudo_mouseA_subset{area}(:,:,1:10),3,'omitnan'));  % calculate the internal struture of session 1 for a single area in shuffled pseudo-mouse A
            current_structure_mouseA_shuffle(:,:,2) = corr(mean(shuffled_pseudo_mouseA_subset{area}(:,:,11:20),3,'omitnan')); % calculate the internal struture of session 2 for a single area in shuffled pseudo-mouse A
            current_structure_mouseA_shuffle(:,:,3) = corr(mean(shuffled_pseudo_mouseA_subset{area}(:,:,21:30),3,'omitnan')); % calculate the internal struture of session 3 for a single area in shuffled pseudo-mouse A
            current_structure_mouseA_shuffle = mean(current_structure_mouseA_shuffle ,3,'omitnan'); % average internal structure across sessions
            internal_structures_pseudoA_shuffle(:,area) = current_structure_mouseA_shuffle(triu_ind); % vectorize and store the internal structure
            
            
             current_structure_mouseB_shuffle = [];
            current_structure_mouseB_shuffle(:,:,1) = corr(mean(shuffled_pseudo_mouseB_subset{area}(:,:,1:10),3,'omitnan'));  % calculate the internal struture of session 1 for a single area in shuffled pseudo-mouse B
            current_structure_mouseB_shuffle(:,:,2) = corr(mean(shuffled_pseudo_mouseB_subset{area}(:,:,11:20),3,'omitnan')); % calculate the internal struture of session 2 for a single area in shuffled pseudo-mouse B
            current_structure_mouseB_shuffle(:,:,3) = corr(mean(shuffled_pseudo_mouseB_subset{area}(:,:,21:30),3,'omitnan')); % calculate the internal struture of session 3 for a single area in shuffled pseudo-mouse B
            current_structure_mouseB_shuffle = mean(current_structure_mouseB_shuffle ,3,'omitnan'); % average internal structure across sessions
            internal_structures_pseudoB_shuffle(:,area) = current_structure_mouseB_shuffle(triu_ind); % vectorize and store the internal structure
        
     
        end
        
        % calculating the similarity between the internal structures of a
        % reference mouse (pseudo-mouse A) and all 720 permutations of the test
        % mouse (pseudo-mouse B)
        similarity_sum = []; % define an empty variable that will store the total similarity (correlation sum) between non-shuffled pseudo-mice
        similarity_sum_shuffled = []; % define an empty variable that will store the total similarity (correlation sum) between shuffled pseudo-mice
        permutations = flipud(perms([1:6])); % define a matrix with all possible permutations
        for perm = 1:size(permutations,1) % loop over permutations
            current_perm = permutations(perm,:); % set the current permutations
            similarity_between_pseudomice = corr(internal_structures_pseudoA,internal_structures_pseudoB(:,current_perm)); % calculate  the correlation between the internal structures of reference and permutated non-shuffled pseudo-mice
            similarity_between_pseudomice_shuffled = corr(internal_structures_pseudoA_shuffle,internal_structures_pseudoB_shuffle(:,current_perm)); % calculate  the correlation between the internal structures of reference and permutated shuffled pseudo-mice
            
            similarity_sum(perm) = sum(diag(similarity_between_pseudomice)); % calculate the sum of correlations between corresponding pseudo-areas for non-shuffled pseudo-mice
            similarity_sum_shuffled(perm) = sum(diag(similarity_between_pseudomice_shuffled)); % calculate the sum of correlations between corresponding pseudo-areas for shuffled pseudo-mice
        end

        [B,I] = max(similarity_sum); % find the permutation with the highest similarity between non-shuffled pseudo-mice
        between_similarity_acc(shuffle,:) = permutations(I,:) == [1:6]; % assess decoder prediction for non-shuffled pseudo-mice
        
        [B,I] = max(similarity_sum_shuffled); % find the permutation with the highest similarity between shuffled pseudo-mice
        between_similarity_acc_shuffled(shuffle,:) = permutations(I,:) == [1:6]; % assess decoder prediction for shuffled pseudo-mice
        
    end
    between_similarity_acc_all_cutoffs(cells_included,:) = sum(between_similarity_acc)./num_shuffles; % store calculate decoder overall performance across all realizations of non-shuflled pseudo-mice
    between_similarity_acc_shuffled_all_cutoffs(cells_included,:) = sum(between_similarity_acc_shuffled)./num_shuffles; % store calculate decoder overall performance across all realizations of shuffled pseudo-mice
end


plt = []; % define empty variable to store plot information for decoder performance
xvalues = [cell_cutoff_list,round(mean(min_cell_num_all_shuffles,'omitnan'))]; % define x axis values
figure('units','normalized','position',[0.35 0.4 0.25 0.375]) % visualize decoder performance as function of the number of cells included
for area = 1:6 % loop over areas
    hold on
    plt(area) = plot(xvalues,between_similarity_acc_all_cutoffs(:,area)*100,'color',colors(area,:),'linewidth',3);
    plot(xvalues,between_similarity_acc_shuffled_all_cutoffs(:,area)*100,'color',[0.2 0.2 0.2]+0.1*(area-1),'linewidth',3)
end
text(0.55, 0.225,['Shuffled pseudo-mice'],'Units','normalized','color',[0.4 0.4 0.4],'fontsize',12)
set(gca,'xtick',[0:100:800],'xticklabels',[0:100:800])
ylim([0 100])
xlim([-10 820])
ylabel('Successful classifications (%)')
xlabel('Number of cells included')
legend(plt,brain_areas(1:6),'Location','best')
legend('boxoff')

%% Figure S7D - Between pseudo-mice permutation decoder - cell count - shuffle control

nat_movie = 1; % natural movie 1

% subsample for each mouse in the neuropixels dataset the neuronal activity during the first 20 movie repeats
subset_population_vectors = {}; % define an empty variable that will store the subsmaple neuronal activity of indevidual mice
for area = 1:6 % loop over areas
    for mouse = 1:size(neuropixels_population_vectors,1) % loop over mice
        if ~isempty(neuropixels_population_vectors{mouse,area,nat_movie}) % test if current mouse was recorded from current area
            subset_population_vectors{mouse,area} = neuropixels_population_vectors{mouse,area,nat_movie}(:,:,1:20); % subset the neuronal activity during the first 20 movie repeats
        end
    end
end
% decode the identity of each visual area based on the
% similarity between the internal structures of the two psedo-mice
cell_cutoff_list = [5,10,25,50,75,100,150,200,250,300,350,400,500,600,700]; % define list of how many cells are included in the analysis
between_similarity_acc_shuffled_all_cutoffs = []; % define an empty variable that will store the decoder results for the suffled internal structures

num_shuffles = 3000; % number of pseudo-mice realizations for each cell count threshold
min_cell_num_all_shuffles = []; % define empty variable to store number of cells sampled in leach realization of pseudo-mice
for cells_included = 1:length(cell_cutoff_list)+1 % loop over cell count thresholds
    
    between_similarity_acc_shuffled = nan(num_shuffles,6); % define empty variable that store decoders predictions for shuffled internal structure
    for shuffle = 1:num_shuffles % loop over realizations
       
       % determine the number of cells to be incuded in the analysis
        if cells_included <= length(cell_cutoff_list) % if below 700 cells than use 'cell_cutoff_list'
            breaker = cell_cutoff_list(cells_included);
        else %used the entire dataset
            breaker = 1;
        end
        
        % create realizations of pseudo-mice until all pseudo-areas contain 
        % more cells relative to 'cell_cutoff_list'
        pseudo_area_cell_num = zeros(2,6); % define a variable for of zeros that will store the number of units in each pseudo-area
        while   breaker > min(pseudo_area_cell_num(:)) % check if all pseudo-areas contain more cells relative to 'cell_cutoff_list'
            
            % split the dataset into two independent group of mice
            pseudo_mouseA_ind = sort(randperm(size(subset_population_vectors,1),size(subset_population_vectors,1)./2)); % indices for group A (pseudo-mouse A)
            pseudo_mouseB_ind = find(~ismember([1:size(subset_population_vectors,1)],pseudo_mouseA_ind)); % indices for group B (pseudo-mouse B)
            
            % for each visual area in each pseudo mouse, pool all units across mice
            % to create 12 pseudo areas (6 areas x 2 pseudo-mice)
            pseudo_mouseA = {}; % define an empty variable that will store the pooled units across mice for each pseudo-area of pseudo-mouse A
            pseudo_mouseB = {}; % define an empty variable that will store the pooled units across mice for each pseudo-area of pseudo-mouse B
            for area = 1:6 % loop over areas
                pseudo_mouseA{area} = cell2mat(subset_population_vectors(pseudo_mouseA_ind,area,nat_movie)); % pool units across mice for a single visual area (pseudo-mouse A)  
                pseudo_mouseB{area} = cell2mat(subset_population_vectors(pseudo_mouseB_ind,area,nat_movie)); % pool units across mice for a single visual area (pseudo-mouse B) 
                pseudo_area_cell_num(1,area) = size(pseudo_mouseA{area},1); % store number of cells in pseudo-area of pseudo-mouse A
                pseudo_area_cell_num(2,area) = size(pseudo_mouseB{area},1); % store number of cells in pseudo-area of pseudo-mouse B
            end
            % if the minimum value in "pseudo_area_cell_num" is higher than
            % the value in "breaker" than the loop breaks and the algorithm continues
        end
        % min_cell_num - the minimal number of cells across pseudo-areas and pseudo-mice.
        % this number will be used to randomly subsample the same number of cells for all pseudo-ares 
        if cells_included <= length(cell_cutoff_list) % when below 700
            min_cell_num = cell_cutoff_list(cells_included);
        else % when using the entire dataset
            min_cell_num = min(pseudo_area_cell_num(:));
            min_cell_num_all_shuffles(shuffle) = min_cell_num;
        end
        
         clc;
        disp(['Performing between Neuropixels pseudo-mice decoding. Cells included: ',num2str(min_cell_num),' | Realization: ',num2str(shuffle),'/',num2str(num_shuffles)])
       % subsampling randomly the same number of units for all pseudo-areas
        pseudo_mouseA_subset = {}; % define an empty variable that will store the subsampled units for each pseudo-area of pseudo-mouse A
        pseudo_mouseB_subset = {}; % define an empty variable that will store the subsampled units for each pseudo-area of pseudo-mouse B
        for area = 1:6 % loop over areas
            % random sampling #min_cell_num of units from each pseudo-area of pseudo-mouse A
            subset_cell_ids_mouseA = sort(randperm(pseudo_area_cell_num(1,area),min_cell_num));
            pseudo_mouseA_subset{area} = pseudo_mouseA{area}(subset_cell_ids_mouseA,:,:);
            
            % random sampling #min_cell_num of units from each pseudo-area of pseudo-mouse B
            subset_cell_ids_mouseB = sort(randperm(pseudo_area_cell_num(2,area),min_cell_num));
            pseudo_mouseB_subset{area} = pseudo_mouseB{area}(subset_cell_ids_mouseB,:,:);
        end
        
        internal_structures_pseudoA_shuffle = []; % define an empty variable the will store the shuffled internal structures of pseudo-mouse A
        internal_structures_pseudoB_shuffle = []; % define an empty variable the will store the shuffled internal structures of pseudo-mouse B
        internal_structures_labels = []; % define an empty variable that will store the label of each internal structure
        triu_ind = boolean(triu(ones(30),1)); % define a boolian matrix with true values in the upper half of the matrix
        for area = 1:6 % loop over areas
            rand_ind = circshift([1:30],randperm(30,1)); % temporally shuffled time bins id for pseudo-mouse A
            current_structure_mouseA_shuffle = corr(mean(pseudo_mouseA_subset{area}(:,rand_ind,:),3,'omitnan')); % shuffle neuronal activity, average across movie repeats and calculate internal structure
            internal_structures_pseudoA_shuffle(:,area) = current_structure_mouseA_shuffle(triu_ind); % vectorize the internal structure of pseudo-mouse A
            
            rand_ind = circshift([1:30],randperm(30,1)); % temporally shuffled time bins id for pseudo-mouse B
            current_structure_mouseB_shuffle = corr(mean(pseudo_mouseB_subset{area}(:,rand_ind,:),3,'omitnan')); % shuffle neuronal activity, average across movie repeats and calculate internal structure
            internal_structures_pseudoB_shuffle(:,area) = current_structure_mouseB_shuffle(triu_ind); % vectorize the internal structure of pseudo-mouse A
        end
        
         % calculating the similarity between the internal structures of a
         % reference mouse (pseudo-mouse A) and all 720 permutations of the test
         % mouse (pseudo-mouse B)
          similarity_sum_shuffled = []; % define an empty variable that will store the total similarity (correlation sum) between shuffled pseudo-mice
         permutations = flipud(perms([1:6])); % define a matrix with all possible permutations
         for perm = 1:size(permutations,1) % loop over permutations
             current_perm = permutations(perm,:); % set the current permutations
             similarity_between_pseudomice_shuffled = corr(internal_structures_pseudoA_shuffle,internal_structures_pseudoB_shuffle(:,current_perm)); % calculate  the correlation between the internal structures of reference and permutated shuffled pseudo-mice
             similarity_sum_shuffled(perm) = sum(diag(similarity_between_pseudomice_shuffled)); % calculate the sum of correlations between corresponding pseudo-areas for shuffled pseudo-mice
         end
      
         [B,I] = max(similarity_sum_shuffled); % find the permutation with the highest similarity between shuffled pseudo-mice
         between_similarity_acc_shuffled(shuffle,:) = permutations(I,:) == [1:6]; % assess decoder prediction for shuffled pseudo-mice
         
    end
    between_similarity_acc_shuffled_all_cutoffs(cells_included,:) = sum(between_similarity_acc_shuffled)./num_shuffles; % store calculate decoder overall performance across all realizations of shuffled pseudo-mice
    
end

plt = [];% define empty variable to store plot information for decoder performance
xvalues = [cell_cutoff_list,round(nanmean(min_cell_num_all_shuffles))]; % define x axis values
figure('units','normalized','position',[0.35 0.4 0.25 0.375]) % visualize decoder performance as function of the number of cells included
for area = 1:6
    hold on
    plt(area)= plot(xvalues,between_similarity_acc_shuffled_all_cutoffs(:,area)*100,'color',colors(area,:),'linewidth',3);
end
set(gca,'xtick',[0:100:800],'xticklabels',[0:100:800])
ylim([0 100])
xlim([-10 820])
ylabel('Successful classifications (%)')
xlabel('Number of cells included')
legend(plt,brain_areas(1:6),'Location','best')
legend('boxoff')

%% Figure S7E - Internal structure in the reduced space for both NM1 and SNM1
clearvars cell
subset_population_vectors = cell(58,6,2); % define an empty variable that will variable that will store the neuronal activity of mice from the 'Function connectivity' group
repeat_num = 30; % number of movie repeats for natural movie 1 in the functional connectivity group
for nat_movie = 1:2 % loop over natural movies
    for area = 1:6 % loop over areas
        for mouse = 1:size(neuropixels_population_vectors_tsne,1) % loop over mice
            if ~isempty(neuropixels_population_vectors_tsne{mouse,area,nat_movie}) && movie_repeats(mouse,1) == repeat_num % check if from functional connectivity group
                subset_population_vectors{mouse,area,nat_movie} = neuropixels_population_vectors_tsne{mouse,area,nat_movie}(:,:,1:10); % subset the neuronal activity during the first 10 movie repeats
            end
        end
    end
end

load('figureS7E.mat','state') % load seed for tsne visualization
rng(state) % set seed
dim_reduce_all = {}; % define empty variable that will store the tsne comp
cells_ind = [];
min_cell_num = 919; % subsample the same number of cells for each visual area
for area = 1:6 % loop over areas
    % dim reduction on neuronal activity during natural movie 1
    current_area = cell2mat(subset_population_vectors(:,area,1)); % subset neuronal activity of a single area
    % random subsampling the same number of cells 
    cells_ind(area,:) = randperm(size(current_area,1),min_cell_num); 
    current_area = current_area(cells_ind(area,:),:,:); 
    current_area_struct = reshape(current_area,[size(current_area,1),size(current_area,2)*size(current_area,3)]); % reshape neuronal activity from 3D into 2D
    clc;[area,1]
    dim_reduce = tsne(current_area_struct','Algorithm','exact','Distance','cosine',...
        'NumDimensions',3,'NumPCAComponents',20,'Perplexity',100);
    dim_reduce_all{1,area} = dim_reduce; % perform tsne dim reduction on the cells dimention
    
     % dim reduction on neuronal activity during shuffled natural movie 1
    current_area = cell2mat(subset_population_vectors(:,area,2)); % subset neuronal activity of a single area
    % random subsampling the same number of cells 
    current_area = current_area(cells_ind(area,:),:,:);
    current_area_struct = reshape(current_area,[size(current_area,1),size(current_area,2)*size(current_area,3)]); % reshape neuronal activity from 3D into 2D
    clc;[area,2]
    dim_reduce = tsne(current_area_struct','Algorithm','exact','Distance','cosine',...
        'NumDimensions',3,'NumPCAComponents',20,'Perplexity',100); % perform tsne dim reduction on the cells dimention
    dim_reduce_all{2,area} = dim_reduce;
end

% visualize the internal structures of natural movie 1 and shuffled natural
% movie 1 across areas
view_list_real = [-127.5 23.6; -49.1 10; -115.5 10;-37.5 30;5.3 12.4;-46.3 13.2];
view_list_shuffled = [-61.9 21.2;-119.1 13.2;-70.7 26.8;-24.7 12.4;-137.9 6.8;-51.9 11.6];
figure('units','normalized','position',[0.2 0.3 0.5 0.3])
for area = 1:6
    subplot(2,6,area)
    scatter3(dim_reduce_all{1,area}(:,1),dim_reduce_all{1,area}(:,2),dim_reduce_all{1,area}(:,3),20,repmat([1:90],[1 10]),'filled')
    view(view_list_real(area,:))
    set(gca,'xtick',[],'ytick',[])
    grid off
    
    if area == 3
        title('Natural movie 1 (original frame sequence)')
    end
    subplot(2,6,area+6)
    scatter3(dim_reduce_all{2,area}(:,1),dim_reduce_all{2,area}(:,2),dim_reduce_all{2,area}(:,3),20,repmat([1:90],[1 10]),'filled')
    set(gca,'xtick',[],'ytick',[])
    view(view_list_shuffled(area,:))
    
    grid off
    if area == 3
        title('Shuffled natural movie 1 (fixed temporally shuffled frame sequence)')
    end
    
end
colormap(new_jet_colormap)

%% Figure S7F - Between pseudo-mice decoder using NM1 compared to SNM1
clearvars cell
subset_population_vectors = cell(58,6,2); % define an empty variable that will variable that will store the neuronal activity of mice from the 'Function connectivity' group
repeat_num = 30; % number of movie repeats for natural movie 1 in the functional connectivity group
for nat_movie = 1:2 % loop over natural movies
    for area = 1:6 % loop over areas
        for mouse = 1:size(neuropixels_population_vectors,1) % loop over mice
            if ~isempty(neuropixels_population_vectors{mouse,area,nat_movie}) && movie_repeats(mouse,1) == repeat_num  % check if from functional connectivity group
                % in order to control for the passage of time the 20 subsequent movie repeats between blocks are subsampled
                if nat_movie == 1 % check if natural movie 1
                    subset_population_vectors{mouse,area,nat_movie} = neuropixels_population_vectors{mouse,area,nat_movie}(:,:,21:40); % subset the neuronal activity during the last and first 10 movie repeats in block A and block B respectivly
                elseif nat_movie == 2 % check if shuffled natural movie 1
                    subset_population_vectors{mouse,area,nat_movie} = neuropixels_population_vectors{mouse,area,nat_movie}(:,:,1:20); % subset the neuronal activity during the first 20 movie repeats
                end
            end
        end
    end
end

% decode the identity of each visual area based on the
% similarity between the internal structures of the two psedo-mice
accurate_classification_movA = []; % define an empty variable that will store the decoder results for the natural movie 1
accurate_classification_movB = []; % define an empty variable that will store the decoder results for the shuffled natural movie 1

num_shuffles = 1000; % number of pseudo-mice realizations
for shuffle = 1:num_shuffles % loop over realizations
     clc;
    disp(['Performing between pseudo-mice decoding. Realization: ',num2str(shuffle),'/',num2str(num_shuffles)])
    
    all_mice_ind = find(movie_repeats(:,1) == repeat_num); % indices of 'functional connectivity' mice
    pseudo_mouseA_ind = all_mice_ind(sort(randperm(length(all_mice_ind),round(length(all_mice_ind)./2)))); % indices for group A (pseudo-mouse A)
    pseudo_mouseB_ind = all_mice_ind(~ismember(all_mice_ind,pseudo_mouseA_ind)); % indices for group B (pseudo-mouse B)
    
    % for each visual area in each pseudo mouse, pool all units across mice
    % to create 12 pseudo areas (6 areas x 2 pseudo-mice)
    pseudo_area_cell_num = []; % define an empty variable that will store the number of units in each pseudo-area
    pseudo_mouseA = {}; % define an empty variable that will store the pooled units across mice for each pseudo-area of pseudo-mouse A
    pseudo_mouseB = {}; % define an empty variable that will store the pooled units across mice for each pseudo-area of pseudo-mouse B
    for nat_movie = 1:2
        for area = 1:6
            pseudo_mouseA{nat_movie,area} = cell2mat(subset_population_vectors(pseudo_mouseA_ind,area,nat_movie));
            pseudo_mouseB{nat_movie,area} = cell2mat(subset_population_vectors(pseudo_mouseB_ind,area,nat_movie));
            cell_num_all_areas(1,area) = size(pseudo_mouseA{nat_movie,area},1);
            cell_num_all_areas(2,area) = size(pseudo_mouseB{nat_movie,area},1);
        end
    end
    
    % min_cell_num - the minimal number of cells across pseudo-areas and pseudo-mice.
    % this number will be used to randomly subsample the same number of cells for all pseudo-ares 
    min_num_cell = min(cell_num_all_areas(:));
    subset_pseudo_mouseA = {}; % define an empty variable that will store the subsampled units for each pseudo-area of pseudo-mouse A
    subset_pseudo_mouseB = {}; % define an empty variable that will store the subsampled units for each pseudo-area of pseudo-mouse B
   
    % random sampling #min_cell_num of units from each pseudo-area
    for area = 1:6
        % define cells id only once for both movies
        rand_cellsA= randperm(size(pseudo_mouseA{1,area},1),min_num_cell);
        rand_cellsB = randperm(size(pseudo_mouseB{1,area},1),min_num_cell);
        
        for nat_movie = 1:2
            % loop over areas
            subset_pseudo_mouseA{nat_movie,area} = pseudo_mouseA{nat_movie,area}(rand_cellsA,:,:);
            subset_pseudo_mouseB{nat_movie,area} = pseudo_mouseB{nat_movie,area}(rand_cellsB,:,:);
        end
    end
    

    % calculating the internal structures for each movie repeat for all
    % pseudo-areas of both example pseudo-mice
    triu_ind = boolean(triu(ones(30),1)); % define a boolian matrix with true values in the upper half of the matrix
    movA_all_structures_pseudoA = []; % define an empty variable the will store the internal structures of natural movie for pseudo-mouse A
    movB_all_structures_pseudoA = []; % define an empty variable the will store the internal structures of shuffled natural movie for pseudo-mouse A
    movA_all_structures_pseudoB = []; % define an empty variable the will store the internal structures of natural movie for pseudo-mouse B
    movB_all_structures_pseudoB = []; % define an empty variable the will store the internal structures of shuffled natural movie for pseudo-mouse B
    for nat_movie = 1:2 % loop over natural movies
        for area = 1:6 % loop over areas
            pseudo_mouseA_struct = corr(mean(subset_pseudo_mouseA{nat_movie,area},3,'omitnan')); % average activity across movie repeats and calculate the internal structure for pseudo-mouse A
            pseudo_mouseB_struct = corr(mean(subset_pseudo_mouseB{nat_movie,area},3,'omitnan')); % average activity across movie repeats and calculate the internal structure for pseudo-mouse B
          
             % vectorize and store the internal structure
            if nat_movie ==1 % check if natural movie 1
                movA_all_structures_pseudoA(:,area) = pseudo_mouseA_struct(triu_ind);
                movA_all_structures_pseudoB(:,area) = pseudo_mouseB_struct(triu_ind);
            elseif nat_movie ==2  % check if shuffled natural movie 1
                movB_all_structures_pseudoA(:,area) = pseudo_mouseA_struct(triu_ind);
                movB_all_structures_pseudoB(:,area) = pseudo_mouseB_struct(triu_ind);
            end
        end
    end
    
     % calculating the similarity between the internal structures of a
    % reference mouse (pseudo-mouse A) and all 720 permutations of the test
    % mouse (pseudo-mouse B)
    permutations = flipud(perms([1:6])); % define a matrix with all possible permutations
    perms_similarity_movA = []; % define an empty variable that will store the total similarity (correlation sum) between non-shuffled pseudo-mice for natural movie 1
    perms_similarity_movB = [];  % define an empty variable that will store the total similarity (correlation sum) between non-shuffled pseudo-mice for shuffled natural movie 1
    
    for perm = 1:size(permutations,1) % loop over permutations
        current_perm = permutations(perm,:);  % set the current permutations
        
        perms_similarity_movA(perm) = sum(diag(corr(movA_all_structures_pseudoA,movA_all_structures_pseudoB(:,current_perm)))); % calculate the correlation sum between the natural movie 1 internal structures of reference and permutated pseudo-mice
        perms_similarity_movB(perm) = sum(diag(corr(movB_all_structures_pseudoA,movB_all_structures_pseudoB(:,current_perm)))); % calculate the correlation sum between the shuffled natural movie 1 internal structures of reference and permutated pseudo-mice
    end
    
    [~,best_perm_movA] = max(perms_similarity_movA); % find the permutation with the highest similarity for natural movie 1
    [~,best_perm_movB] = max(perms_similarity_movB); % find the permutation with the highest similarity for shuffled natural movie 1
    
    accurate_classification_movA(shuffle,:) = permutations(best_perm_movA,:); % assess decoder prediction for natural movie 1
    accurate_classification_movB(shuffle,:) = permutations(best_perm_movB,:); % assess decoder prediction for shuffled natural movie 1
    
end

% calculate decoder overall performance across all realizations
acc_NM1 = 100*sum(accurate_classification_movA==[1:6])./shuffle; % natural movie 1
acc_SNM1 = 100*sum(accurate_classification_movB==[1:6])./shuffle; % shuffled natural movie 1

plt = []; % define empty variable to store plot information for decoder performance
figure % visualize decoder performance
hold on
xticks = [];
for area = 1:6
    xticks(area) = nanmean([0.75+1.5*(area-1),1.25+1.5*(area-1)]);
    plt(1) = bar(0.75+1.5*(area-1),acc_NM1(area),'facecolor',[0.7 0.7 0.7]-0.2,'edgecolor','none','barwidth',0.5);
    plt(2) = bar(1.25+1.5*(area-1),acc_SNM1(area),'facecolor',[0.7 0.7 0.7],'edgecolor','none','barwidth',0.5);
end
plot([minmax(xticks)+[-1.5 0.75]],[100 100]./6,'--','color',[0.2 0.2 0.2],'linewidth',2)
text(0.85,0.16,'Chance','Units','normalized','FontSize',11)
xlim(minmax(xticks)+[-1.5 2.75])
ylim([0 100])
set(gca,'xtick',xticks,'xticklabel',brain_areas(1:6))
ylabel({'Successful classifications ';'between pseudo-mice (%)'})
lgd = legend(plt,{'NM1','SNM1'});
legend('boxoff')
lgd.Position = [0.75 0.85 0.05 0.05];

%% Figure S7G - Between pseudo-mice decoder using NM3 compared to DG

num_shuffle = 1000;
pseudo_mice_decoding_acc = nan(num_shuffle,6);
pseudo_mice_decoding_acc_movie3 = nan(num_shuffle,6);
for shuffle = 1:num_shuffle
    clc;[shuffle]
    
    pseudo_mouseA =  {};
    pseudo_mouseB = {};
    pseudo_mouseA_movie3 =  {};
    pseudo_mouseB_movie3 = {};
    pseudo_cell_num = [];
    for area = 1:6
        current_area_gratings = calcium_excitatory_drifting_gratings{area}';
        current_area_movie3 = calcium_excitatory_population_vectors{area}(:,4);
        pseudo_mouseA_ind = sort(randperm(length(current_area_gratings),round(length(current_area_gratings)/2)));
        pseudo_mouseB_ind = find(~ismember([1:length(current_area_gratings)],pseudo_mouseA_ind));
        
        pseudo_mouseA{area} = cell2mat(current_area_gratings(pseudo_mouseA_ind));
        pseudo_mouseB{area} = cell2mat(current_area_gratings(pseudo_mouseB_ind));
        
        pseudo_mouseA_movie3{area} = cell2mat(current_area_movie3(pseudo_mouseA_ind));
        pseudo_mouseB_movie3{area} = cell2mat(current_area_movie3(pseudo_mouseB_ind));
        
        pseudo_cell_num(:,area) = [size(pseudo_mouseA{area},1);size(pseudo_mouseB{area},1)];
    end
    
    min_num_cell = min(pseudo_cell_num(:));
    %min_num_cell = 150;
    subsampled_pseudo_mouseA = {};
    subsampled_pseudo_mouseB = {};
    subsampled_pseudo_mouseA_movie3 = {};
    subsampled_pseudo_mouseB_movie3 = {};
    for area = 1:6
        rand_cells_pseudoA = sort(randperm(pseudo_cell_num(1,area),min_num_cell));
        rand_cells_pseudoB = sort(randperm(pseudo_cell_num(2,area),min_num_cell));
        
        subsampled_pseudo_mouseA{area} = pseudo_mouseA{area}(rand_cells_pseudoA,:,:,:);
        subsampled_pseudo_mouseB{area} =  pseudo_mouseB{area}(rand_cells_pseudoB,:,:,:);
        
        subsampled_pseudo_mouseA_movie3{area} = pseudo_mouseA_movie3{area}(rand_cells_pseudoA,:,:,:);
        subsampled_pseudo_mouseB_movie3{area} =  pseudo_mouseB_movie3{area}(rand_cells_pseudoB,:,:,:);
    end
    
    all_pseudo_mouseA_structures = [];
    all_pseudo_mouseB_structures = [];
    
    all_pseudo_mouseA_structures_movie3 = [];
    all_pseudo_mouseB_structures_movie3 = [];
    
    for area = 1:6
        current_pseudoA_sorted = [];
        current_pseudoB_sorted = [];
        for ori = 1:8
            for freq = 1:size(subsampled_pseudo_mouseA{area},4)
                current_pseudoA_sorted = [current_pseudoA_sorted,subsampled_pseudo_mouseA{area}(:,ori,:,freq)];
                current_pseudoB_sorted = [current_pseudoB_sorted,subsampled_pseudo_mouseB{area}(:,ori,:,freq)];
            end
        end
        
        pseudo_mouseA_structure = corr(nanmean(current_pseudoA_sorted(:,:,1:2),3),'rows','pairwise');
        pseudo_mouseB_structure = corr(nanmean(current_pseudoB_sorted(:,:,1:2),3),'rows','pairwise');
        
        pseudo_mouseA_structure_movie3 = corr(nanmean(subsampled_pseudo_mouseA_movie3{area},3),'rows','pairwise');
        pseudo_mouseB_structure_movie3 = corr(nanmean(subsampled_pseudo_mouseB_movie3{area},3),'rows','pairwise');
        
        triu_ind = boolean(triu(ones(size(pseudo_mouseA_structure,1)),1));
        all_pseudo_mouseA_structures(:,area) = pseudo_mouseA_structure(triu_ind);
        all_pseudo_mouseB_structures(:,area) = pseudo_mouseB_structure(triu_ind);
        
        all_pseudo_mouseA_structures_movie3(:,area) = pseudo_mouseA_structure_movie3(triu_ind);
        all_pseudo_mouseB_structures_movie3(:,area) = pseudo_mouseB_structure_movie3(triu_ind);
        
        
        
    end
    
    
    permutations = flipud(perms([1:6]));
    similarity_vals = [];
    similarity_vals_movie3 = [];
    for perm = 1:size(permutations,1)
        current_perm = permutations(perm,:);
        similarity_vals(perm) =  sum(diag(corr(all_pseudo_mouseA_structures,...
            all_pseudo_mouseB_structures(:,current_perm))));
        
        similarity_vals_movie3(perm) =  sum(diag(corr(all_pseudo_mouseA_structures_movie3,...
            all_pseudo_mouseB_structures_movie3(:,current_perm))));
    end
    
    [~,I] = max(similarity_vals);
    pseudo_mice_decoding_acc(shuffle,:) = permutations(I,:) == [1:6];
    
    [~,I] = max(similarity_vals_movie3);
    pseudo_mice_decoding_acc_movie3(shuffle,:) = permutations(I,:) == [1:6];
    
end

acc_DG= 100*(sum(pseudo_mice_decoding_acc)./num_shuffle);
acc_NM3 = 100*(sum(pseudo_mice_decoding_acc_movie3)./num_shuffle);

figure
hold on
plt = [];
xticks = [];
for area = 1:6
    xticks(area) = nanmean([0.75+1.5*(area-1),1.25+1.5*(area-1)]);
    plt(1) = bar(0.75+1.5*(area-1),acc_NM3(area),'facecolor',[0.7 0.7 0.7]-0.2,'edgecolor','none','barwidth',0.5);
    plt(2) = bar(1.25+1.5*(area-1),acc_DG(area),'facecolor',[0.7 0.7 0.7],'edgecolor','none','barwidth',0.5);
end
plot([minmax(xticks)+[-1.5 0.75]],[100 100]./6,'--','color',[0.2 0.2 0.2],'linewidth',2)
text(0.85,0.16,'Chance','Units','normalized','FontSize',11)
xlim(minmax(xticks)+[-1.5 2.75])
ylim([0 100])
set(gca,'xtick',xticks,'xticklabel',brain_areas(1:6))
ylabel({'Successful classifications ';'between pseudo-mice (%)'})
lgd = legend(plt,{'NM3','DG'});
legend('boxoff')
lgd.Position = [0.75 0.85 0.05 0.05];


%% Figure S7H - Internal structure stability VS PV stability without normalization
subset_population_vectors = {};
for area = 1:6
    subset_population_vectors{area} = calcium_excitatory_population_vectors{area}(:,1);
end

triu_id = boolean(triu(ones(30),1));
num_shuffles = 1000;

elapsed_session_all_measurments = {};
for area = 1:6
    elapsed_session_pv = nan(num_shuffles,3);
    elapsed_session_struc = nan(num_shuffles,3);
    
    current_area = cell2mat(subset_population_vectors{area});
    
    valid_cells = [nanmean(nanmean(current_area(:,:,1:10),3),2) > 0 &...
        nanmean(nanmean(current_area(:,:,11:20),3),2) > 0 &...
        nanmean(nanmean(current_area(:,:,21:30),3),2) > 0];
    current_area = current_area(valid_cells,:,:);
    
    for shuffle = 1:num_shuffles
        clc;
        disp(['Generating calcium imaging pseudo-mice and calculating'])
        disp(['internal structure and population vectors stability:'])
        disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Realization: ',num2str(shuffle),'/',num2str(num_shuffles)])
        
        
        valid_cells = sort(randperm(size(current_area,1),round(size(current_area,1)*0.7)));
        subset_valid_current_area = current_area(valid_cells,:,:);
        
        all_structures = [];
        mean_pv = [];
        
        for half1 = 1:6
            current_half1 = nanmean(subset_valid_current_area(:,:,[1:5]+5*(half1-1)),3);
            for half2 = 1:6
                current_half2 = nanmean(subset_valid_current_area(:,:,[1:5]+5*(half2-1)),3);
                
                current_structure = corr(current_half1,current_half2);
                mean_pv(half1,half2) = nanmean(diag(current_structure));
                
                if half1==half2
                    all_structures(:,half1) =  current_structure(triu_id);
                end
            end
        end
        
        
        structure_stability = corr(all_structures);
        
        elapsed_session_pv(shuffle,1) = nanmean([mean_pv(1,2),mean_pv(3,4),mean_pv(5,6)]);
        elapsed_session_pv(shuffle,2) = nanmean([nanmean(nanmean(mean_pv(1:2,3:4))),nanmean(nanmean(mean_pv(3:4,5:6)))]);
        elapsed_session_pv(shuffle,3) = nanmean([nanmean(nanmean(mean_pv(1:2,5:6)))]);
        
        elapsed_session_struc(shuffle,1) = nanmean([structure_stability(1,2),structure_stability(3,4),structure_stability(5,6)]);
        elapsed_session_struc(shuffle,2) = nanmean([nanmean(nanmean(structure_stability(1:2,3:4))),nanmean(nanmean(structure_stability(3:4,5:6)))]);
        elapsed_session_struc(shuffle,3) = nanmean([nanmean(nanmean(structure_stability(1:2,5:6)))]);
        
    end
    elapsed_session_all_measurments(area,:) = {elapsed_session_pv,elapsed_session_struc};
    
end

plt = [];
figure('units','normalized','position',[0.35 0.4 0.25 0.35])
for area = 1:6
    current_area_pv = elapsed_session_all_measurments{area,1};
    current_area_structure = elapsed_session_all_measurments{area,2};
    
    
    mean_struct = nanmean(current_area_structure);
    std_struct = nanstd(current_area_structure);
    
    mean_pv_corr = nanmean(current_area_pv);
    std_pv_corr  = nanstd(current_area_pv);
    
    
    hold on
    plt(area) = errorbar(mean_pv_corr ,std_pv_corr ,'o','color',[0.2 0.2 0.2]+0.1*(area-1),...
        'markerfacecolor',[0.2 0.2 0.2]+0.1*(area-1),'capsize',0,'linestyle','-','linewidth',3,'markersize',3);
    
    plt(area+6) = errorbar(mean_struct,std_struct,'o','color',colors(area,:),...
        'markerfacecolor',colors(area,:),'capsize',0,'linestyle','-','linewidth',3,'markersize',3);
    
end
xlim([0.5 3.5])
ylim([0.3 1])
set(gca,'xtick',1:3,'xticklabel',{'Within session','Proximal sessions','Distal sessions'})
ylabel('Correlation (normalized)')
legend(plt,brain_areas(1:6),'Location','best')
legend('boxoff')

%% Figure S7I - Internal structure stability VS PV stability when shuffling cells� identities
subset_population_vectors = {};
for area = 1:6
    subset_population_vectors{area} = calcium_excitatory_population_vectors{area}(:,1);
end

triu_id = boolean(triu(ones(30),1));
num_shuffles = 1000;

elapsed_session_all_measurments = {};
for area = 1:6
    elapsed_session_pv = nan(num_shuffles,3);
    elapsed_session_struc = nan(num_shuffles,3);
    
    current_area = cell2mat(subset_population_vectors{area});
    
    valid_cells = [nanmean(nanmean(current_area(:,:,1:10),3),2) > 0 &...
        nanmean(nanmean(current_area(:,:,11:20),3),2) > 0 &...
        nanmean(nanmean(current_area(:,:,21:30),3),2) > 0];
    current_area = current_area(valid_cells,:,:);
    
    for shuffle = 1:num_shuffles
        clc;
        disp(['Generating calcium imaging pseudo-mice and calculating'])
        disp(['internal structure and population vectors stability:'])
        disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Realization: ',num2str(shuffle),'/',num2str(num_shuffles)])
        
        
        valid_cells = sort(randperm(size(current_area,1),round(size(current_area,1)*0.7)));
        subset_valid_current_area = current_area(valid_cells,:,:);
        
        
        rand_ind = [];
        for half = 1:6
            rand_ind(half,:) = randperm(size(subset_valid_current_area,1));
        end
        
        all_structures = [];
        mean_pv  = [];
        for half1 = 1:6
            current_half1 = nanmean(subset_valid_current_area(rand_ind(half1,:),:,[1:5]+5*(half1-1)),3);
            for half2 = 1:6
                current_half2 = nanmean(subset_valid_current_area(rand_ind(half2,:),:,[1:5]+5*(half2-1)),3);
                
                current_structure = corr(current_half1,current_half2);
                mean_pv(half1,half2) = nanmean(diag(current_structure));
                
                if half1==half2
                    all_structures(:,half1) =  current_structure(triu_id);
                end
            end
        end
        
        
        structure_stability = corr(all_structures);
        
        elapsed_session_pv(shuffle,1) = nanmean([mean_pv(1,2),mean_pv(3,4),mean_pv(5,6)]);
        elapsed_session_pv(shuffle,2) = nanmean([nanmean(nanmean(mean_pv(1:2,3:4))),nanmean(nanmean(mean_pv(3:4,5:6)))]);
        elapsed_session_pv(shuffle,3) = nanmean([nanmean(nanmean(mean_pv(1:2,5:6)))]);
        
        elapsed_session_struc(shuffle,1) = nanmean([structure_stability(1,2),structure_stability(3,4),structure_stability(5,6)]);
        elapsed_session_struc(shuffle,2) = nanmean([nanmean(nanmean(structure_stability(1:2,3:4))),nanmean(nanmean(structure_stability(3:4,5:6)))]);
        elapsed_session_struc(shuffle,3) = nanmean([nanmean(nanmean(structure_stability(1:2,5:6)))]);
        
    end
    elapsed_session_all_measurments(area,:) = {elapsed_session_pv,elapsed_session_struc};
    
end

plt = [];
figure('units','normalized','position',[0.35 0.4 0.25 0.35])
for area = 1:6
    current_area_pv = elapsed_session_all_measurments{area,1};
    current_area_structure = elapsed_session_all_measurments{area,2};
    
    
    mean_struct = nanmean(current_area_structure);
    std_struct = nanstd(current_area_structure);
    
    mean_pv_corr = nanmean(current_area_pv);
    std_pv_corr  = nanstd(current_area_pv);
    
    
    hold on
    plt(area) = errorbar(mean_pv_corr ,std_pv_corr ,'o','color',[0.2 0.2 0.2]+0.1*(area-1),...
        'markerfacecolor',[0.2 0.2 0.2]+0.1*(area-1),'capsize',0,'linestyle','-','linewidth',3,'markersize',3);
    
    plt(area+6) = errorbar(mean_struct,std_struct,'o','color',colors(area,:),...
        'markerfacecolor',colors(area,:),'capsize',0,'linestyle','-','linewidth',3,'markersize',3);
    
end
xlim([0.5 3.5])
ylim([-0.05 1])
set(gca,'xtick',1:3,'xticklabel',{'Within session','Proximal sessions','Distal sessions'})
ylabel('Correlation (normalized)')
legend(plt,brain_areas(1:6),'Location','best')
legend('boxoff')

%% Figure S7K - Internal structure stability VS Signal correlation stability
subset_population_vectors = {};
for area = 1:6
    subset_population_vectors{area} = calcium_excitatory_population_vectors{area}(:,1);
end

triu_id = boolean(triu(ones(30),1));
num_shuffles = 1000;

elapsed_session_all_measurments = {};
for area = 1:6
    elapsed_session_signal = nan(num_shuffles,3);
    elapsed_session_struc = nan(num_shuffles,3);
    
    current_area = cell2mat(subset_population_vectors{area});
    
    valid_cells = [nanmean(nanmean(current_area(:,:,1:10),3),2) > 0 &...
        nanmean(nanmean(current_area(:,:,11:20),3),2) > 0 &...
        nanmean(nanmean(current_area(:,:,21:30),3),2) > 0];
    current_area = current_area(valid_cells,:,:);
    
    for shuffle = 1:num_shuffles
        clc;
        disp(['Generating calcium imaging pseudo-mice and calculating'])
        disp(['internal structure and signal correlation stability:'])
        disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Realization: ',num2str(shuffle),'/',num2str(num_shuffles)])
        
        
        valid_cells = sort(randperm(size(current_area,1),round(size(current_area,1)*0.7)));
        subset_valid_current_area = current_area(valid_cells,:,:);
        signal_triu = boolean(triu(ones(length(valid_cells))));
        
        all_structures = [];
        signal_corr = [];
        
        for half1 = 1:6
            current_half1 = nanmean(subset_valid_current_area(:,:,[1:5]+5*(half1-1)),3);
            for half2 = 1:6
                current_half2 = nanmean(subset_valid_current_area(:,:,[1:5]+5*(half2-1)),3);
                
                if half1==half2
                    current_structure = corr(current_half1,current_half2);
                    current_signal = corr(current_half1',current_half2');
                    
                    all_structures(:,half1) =  current_structure(triu_id);
                    signal_corr(:,half1) = current_signal(signal_triu);
                end
            end
        end
        
        
        structure_stability = corr(all_structures);
        signal_stability = corr(signal_corr,'rows','pairwise');
        
        elapsed_session_signal(shuffle,1) = nanmean([signal_stability(1,2),signal_stability(3,4),signal_stability(5,6)]);
        elapsed_session_signal(shuffle,2) = nanmean([nanmean(nanmean(signal_stability(1:2,3:4))),nanmean(nanmean(signal_stability(3:4,5:6)))]);
        elapsed_session_signal(shuffle,3) = nanmean([nanmean(nanmean(signal_stability(1:2,5:6)))]);
        
        elapsed_session_struc(shuffle,1) = nanmean([structure_stability(1,2),structure_stability(3,4),structure_stability(5,6)]);
        elapsed_session_struc(shuffle,2) = nanmean([nanmean(nanmean(structure_stability(1:2,3:4))),nanmean(nanmean(structure_stability(3:4,5:6)))]);
        elapsed_session_struc(shuffle,3) = nanmean([nanmean(nanmean(structure_stability(1:2,5:6)))]);
        
    end
    elapsed_session_all_measurments(area,:) = {elapsed_session_signal,elapsed_session_struc};
    
end

plt = [];
figure('units','normalized','position',[0.35 0.4 0.25 0.35])
for area = 1:6
    current_area_signal = elapsed_session_all_measurments{area,1};
    current_area_structure = elapsed_session_all_measurments{area,2};
    
    
    elapsed_session_struc_norm = current_area_structure./current_area_structure(:,1);
    elapsed_session_signal_norm = current_area_signal./current_area_signal(:,1);
    
    mean_struct = nanmean(elapsed_session_struc_norm);
    std_struct = nanstd(elapsed_session_struc_norm);
    
    mean_signal_corr = nanmean(elapsed_session_signal_norm);
    std_signal_corr  = nanstd(elapsed_session_signal_norm);
    
    
    
    hold on
    plt(area) = errorbar(mean_signal_corr ,std_signal_corr ,'o','color',[0.2 0.2 0.2]+0.1*(area-1),...
        'markerfacecolor',[0.2 0.2 0.2]+0.1*(area-1),'capsize',0,'linestyle','-','linewidth',3,'markersize',3);
    
    plt(area+6) = errorbar(mean_struct,std_struct,'o','color',colors(area,:),...
        'markerfacecolor',colors(area,:),'capsize',0,'linestyle','-','linewidth',3,'markersize',3);
    
end
xlim([0.5 3.5])
ylim([0.55 1])
set(gca,'xtick',1:3,'xticklabel',{'Within session','Proximal sessions','Distal sessions'})
ylabel('Correlation (normalized)')
legend(plt,brain_areas(1:6),'Location','best')
legend('boxoff')


%% Figure S7L - Internal structure stability VS Signal correlation stability as function of cell count for area LM

subset_population_vectors = {};
for area = 1:6
    subset_population_vectors{area} = calcium_excitatory_population_vectors{area}(:,1);
end

triu_id = boolean(triu(ones(30),1));

cell_count_list_full = [25, 50, 100, 250,500,1000,2000,4000,6000, 7900];
elapsed_session_all_measurments = {};

area = 2;
current_area = cell2mat(subset_population_vectors{area});

valid_cells = [nanmean(nanmean(current_area(:,:,1:10),3),2) > 0 &...
    nanmean(nanmean(current_area(:,:,11:20),3),2) > 0 &...
    nanmean(nanmean(current_area(:,:,21:30),3),2) > 0];
current_area = current_area(valid_cells,:,:);
cell_count_list = cell_count_list_full(cell_count_list_full<=size(current_area,1));

shuffle_num = 1000;
elapsed_session_struc = nan(length(cell_count_list),3,shuffle_num);
elapsed_session_signal = nan(length(cell_count_list),3,shuffle_num);

for cell_count = 1:length(cell_count_list)
    
    for shuffle = 1:shuffle_num
        clc;[area,cell_count,shuffle]
        
        valid_cells = sort(randperm(size(current_area,1),cell_count_list_full(cell_count)));
        subset_valid_current_area = current_area(valid_cells,:,:);
        
        all_structures = [];
        all_signal_corr = [];
        
        for half1 = 1:6
            current_half1 = nanmean(subset_valid_current_area(:,:,[1:5]+5*(half1-1)),3);
            for half2 = 1:6
                current_half2 = nanmean(subset_valid_current_area(:,:,[1:5]+5*(half2-1)),3);
                
                if half1==half2
                    current_structure = corr(current_half1,current_half2);
                    all_structures(:,half1) =  current_structure(triu_id);
                    
                    signal_corr = corr(current_half1',current_half2');
                    
                    signal_triu_id = boolean(triu(ones(size(signal_corr,1)),1));
                    all_signal_corr(:,half1) = signal_corr(signal_triu_id);
                end
                
            end
        end
        
        structure_stability = corr(all_structures);
        signal_corr_stability = corr(all_signal_corr,'rows','pairwise');
        
        elapsed_session_struc(cell_count,1,shuffle) = nanmean([structure_stability(1,2),structure_stability(3,4),structure_stability(5,6)]);
        elapsed_session_struc(cell_count,2,shuffle) = nanmean([nanmean(nanmean(structure_stability(1:2,3:4))),nanmean(nanmean(structure_stability(3:4,5:6)))]);
        elapsed_session_struc(cell_count,3,shuffle) = nanmean([nanmean(nanmean(structure_stability(1:2,5:6)))]);
        
        elapsed_session_signal(cell_count,1,shuffle) = nanmean([signal_corr_stability(1,2),signal_corr_stability(3,4),signal_corr_stability(5,6)]);
        elapsed_session_signal(cell_count,2,shuffle) = nanmean([nanmean(nanmean(signal_corr_stability(1:2,3:4))),nanmean(nanmean(signal_corr_stability(3:4,5:6)))]);
        elapsed_session_signal(cell_count,3,shuffle) = nanmean([nanmean(nanmean(signal_corr_stability(1:2,5:6)))]);
    end
end

plt = [];
figure('units','normalized','position',[0.35 0.4 0.25 0.35])

curent_area_pv = nanmean(elapsed_session_signal,3);
elapsed_session_pv_norm = [curent_area_pv./curent_area_pv(:,1,:)];


curent_area_struct = nanmean(elapsed_session_struc,3);
elapsed_session_struct_norm = [curent_area_struct./curent_area_struct(:,1,:)];


for cell_count = 1:size(elapsed_session_struct_norm,1)
    hold on
    plot(elapsed_session_pv_norm(cell_count,:),'color',[0.2 0.2 0.2] + 0.075*(cell_count-1),'linewidth',3)
    plt(cell_count)=plot(elapsed_session_struct_norm(cell_count,:),'color',newmap3(end-6*(cell_count-1),:),'linewidth',3);
    
end
legend(plt,num2str(cell_count_list'),'Location','best')
legend('boxoff')
xlim([0.75 3.25])
set(gca,'xtick',1:3,'xticklabel',{'Within session','Proximal sessions','Distal sessions'})
ylabel('Correlation (normalized)')
