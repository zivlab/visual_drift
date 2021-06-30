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
brain_areas = {'VISp','VISl','VISal','VISpm','VISrl','VISam','LGd','LP','CA1','CA2','CA3','DG'};

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
        
        % Define binnization parameters for 'Natural movie 1' and 'Shuffled natural movie 1':
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
            
            
        % Define binnization parameters for 'Natural movie 3':
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
    results_path = ['E:\daniel_master\AllenBrainObservatory\calcium_imaging\results_files\excitatory3\',brain_areas{area},'\'];
    mat_list = dir([results_path,'*.mat']);
    mat_list = {mat_list.name};

    calcium_population_vectors_across_mice = {};
    calcium_drifting_gratings_across_mice = {};
    calcium_spont_population_vectors_across_mice = {};
    mean_cell_num_across_mice = [];
    binned_running_speed_across_mice = {};
    binned_pupil_size_across_mice = {};
    imaging_depth_all_mice = [];
    subject_id_all_mice = [];
    raw_calcium_population_vectors_across_mice = {};
    sorted_mouse_age = [];
    for file = 1:length(mat_list)
        clc;
        disp(['Loading calcium imaging excitatory data:'])
        disp(['Area: ',num2str(area),'\',num2str(6)])
        disp(['Mouse: ',num2str(file),'\',num2str(length(mat_list))])
        load([results_path,mat_list{file}])
        
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
        
        %tsne visualization
        frames_rate = 30;
        repeats_movie2 = 30;
        movie2_length = 30;
        movie2_frames = frames_rate * movie2_length;
        movie2_bin_size = 10;
        binned_movie2 = ones(movie2_frames,1);
        movie2_bin_edges = 1:movie2_bin_size:movie2_frames;
        for bin = 1:length(movie2_bin_edges)
            binned_movie2(movie2_bin_edges(bin):movie2_bin_edges(bin)+movie2_bin_size-1) = bin;
        end
        binned_movie_repeated2 = repmat(binned_movie2,[repeats_movie2,1]);
        
        frames_rate = 30;
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
        
        % NM3 VS DG decoding
        frames_rate = 30;
        repeats_movie3_decode = 10;
        movie3_bin_size_decode = 90;
        binned_movie3_decode = ones(movie3_frames,1);
        movie3_bin_edges_decode = 1:movie3_bin_size_decode:movie3_frames;
        for bin = 1:length(movie3_bin_edges_decode)
            binned_movie3_decode(movie3_bin_edges_decode(bin):movie3_bin_edges_decode(bin)+movie3_bin_size_decode-1) = bin;
        end
        binned_movie_repeated3_decode = repmat(binned_movie3_decode,[repeats_movie3_decode,1]);
        
        repeats_movie = [repeats_movie1,repeats_movie3,repeats_movie2,repeats_movie3];
        movie_frames = [movie1_frames,movie3_frames,movie2_frames,movie3_frames];
        movie_bin_edges = {movie1_bin_edges,movie3_bin_edges,movie2_bin_edges,movie3_bin_edges_decode};
        binned_movie_repeated = {binned_movie_repeated1,binned_movie_repeated3,binned_movie_repeated2,binned_movie_repeated3_decode};
        natural_movie1_traces = reshape(united_traces_days_events,size(united_traces_days_events,1),27000);
        natural_movie1_traces_sessA = filtered_traces_days_events{1,1};
        natural_movie3_traces = filtered_traces_days_events{1,2};
        natural_movie2_traces = filtered_traces_days_events{3,2};
        natural_movie_traces = {natural_movie1_traces,natural_movie3_traces,natural_movie1_traces,natural_movie3_traces};
        
        
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
        if ~isempty(natural_movie_running{1,2})
            natural_movie3_running = natural_movie_running{1,2}';
        else
            natural_movie3_running = nan(1,36000);
        end
        natural_movie_running = {natural_movie1_running,natural_movie3_running};
        
        
        natural_movie1_pupil= cell2mat(natural_movie_pupil_sorted(:,1))';
        if ~isempty(natural_movie_pupil{1,2})
            natural_movie3_pupil = natural_movie_pupil{1,2}';
        else
            natural_movie3_pupil = nan(1,36000);
        end
        natural_movie_pupil = {natural_movie1_pupil,natural_movie3_pupil};
        
        
        % binization for natural movies
        for nat_movie = 1:4
            pop_vector_info_trials = [];
            binned_running_speed = [];
            binned_pupil_size = [];
            sub = 1;
            for repeat = 1:repeats_movie(nat_movie)
                frames_temp = [1:movie_frames(nat_movie)] + (movie_frames(nat_movie)*(repeat-1));
                current_repeat = natural_movie_traces{nat_movie}(:,frames_temp);
                if nat_movie <3
                    current_repeat_running = natural_movie_running{nat_movie}(:,frames_temp);
                    current_repeat_pupil = natural_movie_pupil{nat_movie}(:,frames_temp);
                end
                for bin = 1:length(movie_bin_edges{nat_movie})
                    pop_vector_info_trials(:,bin,sub) = nanmean(current_repeat(:,binned_movie_repeated{nat_movie}(frames_temp) == bin),2);
                    if nat_movie <3
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
        raw_united_traces_days_sorted = reshape(raw_united_traces_days,[size(raw_united_traces_days,1),size(raw_united_traces_days,2)*size(raw_united_traces_days,3)]);
        nat_movie = 1;
        raw_pop_vector_info_trials = [];
        sub = 1;
        for repeat = 1:repeats_movie(nat_movie)
            frames_temp = [1:movie_frames(nat_movie)] + (movie_frames(nat_movie)*(repeat-1));
            current_repeat =  raw_united_traces_days(:,frames_temp);
            
            for bin = 1:length(movie_bin_edges{nat_movie})
                raw_pop_vector_info_trials(:,bin,sub) = nanmean(current_repeat(:,binned_movie_repeated{nat_movie}(frames_temp) == bin),2);
                
            end
            sub = sub + 1;
        end
        
        raw_calcium_population_vectors_across_mice(file) = {raw_pop_vector_info_trials};
        calcium_drifting_gratings_across_mice(file) = filtered_traces_days_events(1,4);
        
        mean_cell_num_across_mice(file,:) = cellfun(@length,cell_registration);
        imaging_depth_all_mice(file) = imaging_depth;
        subject_id_all_mice(file) = subject_id;
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

%% Load colormaps
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

%% Figure 1B - calcium imaging excitatory Cre line cell counts
average_cell_count = nan(94,6);
for area = 1:6
    current_area = nanmean(calcium_excitatory_cell_count{area},2);
    num_mice(area) = length(current_area);
    average_cell_count(1:num_mice(area),area) = current_area;
end

figure
figure_boxplot(average_cell_count(:,1:6))
set(gca,'xticklabels',brain_areas(1:6))
xtickangle(90)
ylabel('Cell count')

%% Figure 1D - Single cell tuning examples - seconds 3,9,14

nat_movie = 1;
area_list = [6,2,2];
mouse_list = [6,57,44];
cell_list = [40,21,78];
sub = 1;
for current_cell = 1:3
    current_unit = squeeze(calcium_excitatory_population_vectors{area_list(current_cell)}{mouse_list(current_cell),nat_movie}(cell_list(current_cell),:,:)*30)';
    smooth_current_unit = [];
    for repeat = 1:size(current_unit,1)
        smooth_current_unit(repeat,:) = imgaussfilt(current_unit(repeat,:),2);
    end
    norm_current_unit = smooth_current_unit ./ max(smooth_current_unit,[],2);
    [row,col] = find(norm_current_unit == 1);
    [B,I] = sort(row);
    subplot(2,3,sub)
    imagesc(norm_current_unit)
    hold on
    plot(xlim, [10 10]+0.5,'--',...
        'linewidth',2,'color','w')
    plot(xlim, [20 20]+0.5,'--',...
        'linewidth',2,'color','w')
    if sub == 1
        ylabel('Movie repeat')
    end
    colormap(newmap3)
    
    subplot(2,3,sub+3)
    mean_current_unit = [];
    mean_current_unit(1,:) = nanmean(current_unit(1:10,:));
    mean_current_unit(2,:) = nanmean(current_unit(11:20,:));
    mean_current_unit(3,:) = nanmean(current_unit(21:30,:));
    
    
    hold on
    
    plot(mean_current_unit(1,:),'-','markersize',20,'linewidth',1.5,'color',[0.4 0.4 0.4])
    plot(mean_current_unit(2,:)','-','markersize',20,'linewidth',1.5,'color',[0.6 0.6 0.6])
    plot(mean_current_unit(3,:)','-','markersize',20,'linewidth',1.5,'color',[0.8 0.8 0.8])
    [r,p] = corr(mean_current_unit(1,:)',mean_current_unit(2,:)');
    %     text(0.6, 0.95,['r = ',num2str(r)],'Units','normalized','color',[0.25 0.25 0.25])
    %     text(0.6, 0.85,['p = ',num2str(p)],'Units','normalized','color',[0.25 0.25 0.25])
    xlim([1 30])
    ylim([0 18])
    if sub == 1
        ylabel('Mean activity rate (events/sec)')
        legend({'Session 1','Session 2','Session 3'})
        legend('boxoff')
    elseif sub == 2
        xlabel('Time in movie (sec)')
    end
    
    sub = sub + 1;
end

%% Figure 1F - neuropixels cell counts
valid_cell_counts = neuropixels_cell_count(:,:,1);
valid_cell_counts(valid_cell_counts==0) = NaN;

figure
figure_boxplot(valid_cell_counts(:,1:6))
set(gca,'xticklabels',brain_areas(1:6))
xtickangle(90)
ylabel('Cell count')

%% Figure 1H -  Single cell tuning examples - seconds 3,9,14

nat_movie = 1;
area_list = [2,1,2];
mouse_list = [53,54,58];
unit_list = [22,39,7];
sub = 1;
for unit = 1:3
    current_unit = squeeze(neuropixels_population_vectors{mouse_list(unit),area_list(unit),nat_movie}(unit_list(unit),:,:)*30)';
    smooth_current_unit = [];
    for repeat = 1:size(current_unit,1)
        smooth_current_unit(repeat,:) = imgaussfilt(current_unit(repeat,:),2);
    end
    norm_current_unit = smooth_current_unit ./ max(smooth_current_unit,[],2);
    [row,col] = find(norm_current_unit == 1);
    [B,I] = sort(row);
    subplot(2,3,sub)
    imagesc(norm_current_unit)
    hold on
    title(['Unit #',num2str(unit)])
    plot(xlim, [30 30]+0.5,...
        'linewidth',2,'color','w')
    
    if sub == 1
        ylabel('Movie repeat')
    end
    colormap(newmap3)
    
    subplot(2,3,sub+3)
    mean_current_unit = [];
    mean_current_unit(1,:) = nanmean(current_unit(1:30,:));
    mean_current_unit(2,:) = nanmean(current_unit(31:60,:));
    
    mean_current_unit(3,:) = nanmean(current_unit);
    
    hold on
    plot(mean_current_unit(1,:),'-','markersize',20,'linewidth',1.5,'color',[0.4 0.4 0.4])
    plot(mean_current_unit(2,:)','-','markersize',20,'linewidth',1.5,'color',[0.6 0.6 0.6])
    
    xlim([1 30])
    if sub == 1
        ylabel('Mean activity rate (spike/sec)')
        
    elseif sub == 2
        xlabel('Time in movie (sec)')
    elseif sub == 3
        legend({'Block A','Block B'})
        legend('boxoff')
    end
    
    sub = sub + 1;
end


%% Figure 1D + Figure 1H - example frames from Natural movie 1

filename = ['D:\daniel-master\AllenBrainObservatory\Advanced\experimental_databoc\ophys_experiment_data\569792817.nwb'];
movie_name = 'natural_movie_one_image_stack';
nat_movie_frames = h5read(filename,['/stimulus/templates/',movie_name,'/data']);
relevent_frames = [120, 270, 420, 540, 740, 870];
time_in_movie = [3,9,14,19,24,29];

figure('Units','Normalized','Position',[0.2 0.4 0.7 0.15])
for frame = 1:length(relevent_frames)
    subplot(1,6,frame)
    imagesc(nat_movie_frames(:,:,relevent_frames(frame))')
    if frame == 1
        title([num2str(time_in_movie(frame)),'rd second'])
    else
        title([num2str(time_in_movie(frame)),'th second'])
    end
    set(gca,'xtick',[],'ytick',[])
    colormap(gray)
end


%% Figure 2A - PV correlation across repeats for a single mouse recorded from area PM
nat_movie = 1;
area = 4;
mouse = 53;
current_mouse = neuropixels_population_vectors{mouse,area,nat_movie}(:,:,1:30);
current_movie_pv = [];
for repeat = 1:10
    current_movie_pv = [current_movie_pv,current_mouse(:,:,repeat)];
end

current_movie_pv_corr = corr(current_movie_pv);
figure % main
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

current_mouse = neuropixels_population_vectors{mouse,area,nat_movie}(:,:,1:30);
all_repeat_pv = [];
sub = 1;
for repeat1 = 1:30
    for repeat2 = 1:30
        if ~(repeat2 == repeat1)
            all_repeat_pv(:,:,sub) = corr(current_mouse(:,:,repeat1),current_mouse(:,:,repeat2));
            sub = sub + 1;
        end
    end
end

% inset
subplot(1,2,2,'units','normalized','position',[0.6 0.7 0.15 0.2])

imagesc(nanmean(all_repeat_pv,3))
title('Time in movie')
set(gca,'xtick',[],'ytick',[])
colormap(newmap3)

%% Figure 2 B - Mean PV correlation for single PM example mouse
nat_movie = 1;
area = 4;
mouse = 53;
current_mouse = neuropixels_population_vectors{mouse,area,nat_movie}(:,:,1:30);

mean_current_movie_pv_corr = [];
for repeat1 = 1:30
    for repeat2 = 1:30
        current_pv = corr(current_mouse(:,:,repeat1),current_mouse(:,:,repeat2));
        mean_current_movie_pv_corr(repeat1,repeat2) = nanmean(diag(current_pv));
    end
end

figure
mean_current_movie_pv_corr(boolean(eye(size(mean_current_movie_pv_corr,1)))) = NaN;
mean_current_movie_pv_corr(isnan(mean_current_movie_pv_corr(:))) = max(mean_current_movie_pv_corr(:));
imagesc(mean_current_movie_pv_corr)
colormap(newmap3)
xlabel('Movie repeat')
ylabel('Movie repeat')
title('Single animal example')
cb = colorbar();
cb.Label.String = 'PV correlation';
cb.FontSize = 12;


%% Figure 2C - Mean PV correlation as fuction of interval for single PM example mouse
nat_movie = 1;
area = 4;
mouse = 53;
current_mouse = neuropixels_population_vectors{mouse,area,nat_movie}(:,:,1:30);

mean_current_movie_pv_corr = [];
for repeat1 = 1:30
    for repeat2 = 1:30
        current_pv = corr(current_mouse(:,:,repeat1),current_mouse(:,:,repeat2));
        mean_current_movie_pv_corr(repeat1,repeat2) = nanmean(diag(current_pv));
    end
end

mean_pv_elapse_repeat = [];
elapse_repeat_values_mat = [];
for diagonal = 1:size(mean_current_movie_pv_corr,1)-1
    diag_values = diag(mean_current_movie_pv_corr,diagonal);
    elapse_repeat_values_mat = [elapse_repeat_values_mat;[diag_values,ones(length(diag_values),1)*diagonal]];
    mean_pv_elapse_repeat(diagonal) = nanmean(diag_values);
end

figure
hold on
scatter(elapse_repeat_values_mat(:,2),elapse_repeat_values_mat(:,1),25,[0.7 0.7 0.7],'filled')
plot(mean_pv_elapse_repeat,'color',newmap3(1,:),'linewidth',4)
text(0.05,0.075,'1 repeat = 30 seconds','Units','normalized','FontSize',11)
ylabel('PV correlation')
xlabel('Elapsed time (# of movie repeats)')
title('Single animal example')
set(gca,'xtick',[1 10 20 29],'ytick',[0.72:0.04:0.92])

%% Figure 2 D - Mean PV correlation matrices across mice and areas
nat_movie = 1;
cell_cutoff = 15;
num_repeats = 30;
mean_pv_corr_all_areas = {};
clim = [];
for area = 1:6
    valid_mice = neuropixels_cell_count(:,area,nat_movie) >= cell_cutoff & movie_repeats(:,nat_movie) == num_repeats;
    current_area = neuropixels_population_vectors(valid_mice,area,nat_movie);
    mean_pv_corr_all_mice = [];
    
    for mouse = 1:length(current_area)
        clc;
        disp(['Calculating PV correlation between movie repeats for Neuropixels data:'])
        disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        current_mouse = current_area{mouse};
        current_mouse_blockA =  current_mouse(:,:,1:30);
        current_mouse_blockB =  current_mouse(:,:,31:60);
        
        current_mouse_mean_pv_corr = [];
        for repeat1 = 1:30
            for repeat2 = 1:30
                current_pv_blockA = corr(current_mouse_blockA(:,:,repeat1),current_mouse_blockA(:,:,repeat2));
                current_mouse_mean_pv_corr(repeat1,repeat2,1) = nanmean(diag(current_pv_blockA));
                
                current_pv_blockB = corr(current_mouse_blockB(:,:,repeat1),current_mouse_blockB(:,:,repeat2));
                current_mouse_mean_pv_corr(repeat1,repeat2,2) = nanmean(diag(current_pv_blockB));
            end
        end
        mean_pv_corr_all_mice(:,:,mouse) = nanmean(current_mouse_mean_pv_corr,3);
        
    end
    mean_pv_corr_all_areas{area} =  mean_pv_corr_all_mice;
    mean_pv_corr_all_mice = nanmean(mean_pv_corr_all_mice,3);
    mean_pv_corr_all_mice(boolean(eye(size(mean_pv_corr_all_mice,1)))) = NaN;
    clim(area,:) = minmax(mean_pv_corr_all_mice(:)');
end

figure('units','normalized','position',[0.3 0.3 0.26 0.3])
for area = 1:6
    subplot(2,3,area)
    current_area = nanmean(mean_pv_corr_all_areas{area},3);
    current_area(boolean(eye(size(current_area,1)))) = NaN;
    current_area(isnan(current_area(:))) = max(clim(:));
    imagesc(current_area,minmax(clim(:)'))
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
nat_movie = 1;
cell_cutoff = 15;
num_repeats = 30;
mean_pv_elapsed_repeat_per_area = {};
for area = 1:6
    valid_mice = neuropixels_cell_count(:,area,nat_movie) >= cell_cutoff & movie_repeats(:,nat_movie) == num_repeats;
    current_area = neuropixels_population_vectors(valid_mice,area,nat_movie);
    mean_pv_elapsed_repeat = [];
    
    for mouse = 1:length(current_area)
        clc;
        disp(['Calculating PV correlation between movie repeats for Neuropixels data:'])
        disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        current_mouse = current_area{mouse};
        current_mouse_blockA =  current_mouse(:,:,1:30);
        current_mouse_blockB =  current_mouse(:,:,31:60);
        
        current_mouse_mean_pv_corr = [];
        for repeat1 = 1:30
            for repeat2 = 1:30
                current_pv_blockA = corr(current_mouse_blockA(:,:,repeat1),current_mouse_blockA(:,:,repeat2));
                current_mouse_mean_pv_corr(repeat1,repeat2,1) = nanmean(diag(current_pv_blockA));
                
                current_pv_blockB = corr(current_mouse_blockB(:,:,repeat1),current_mouse_blockB(:,:,repeat2));
                current_mouse_mean_pv_corr(repeat1,repeat2,2) = nanmean(diag(current_pv_blockB));
            end
        end
        mean_pv_corr_across_blocks = nanmean(current_mouse_mean_pv_corr,3);
        
        for diagonal = 1:num_repeats-1
            mean_pv_elapsed_repeat(mouse,diagonal) = nanmean(diag(mean_pv_corr_across_blocks,diagonal));
        end
    end
    mean_pv_elapsed_repeat_per_area{area} = mean_pv_elapsed_repeat;
end

pval = [];
df = [];
chi = [];
figure
for area = 1:6
    current_area = mean_pv_elapsed_repeat_per_area{area};
    [pval(area),tbl,stats] = friedman(current_area ,[1],'off');
    df(area) = tbl{2,3};
    chi(area) = tbl{2,5};
    
    mean_across_mice = nanmean(current_area);
    std_across_mice = nanstd(current_area);
    ste_across_mice = std_across_mice./sqrt(size(current_area,1));
    
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

corrected_pval = bonf_holm(pval);

VarNames = {'area','df','chi','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas(1:6)'),df(:),chi(:),pval(:),corrected_pval(:),'VariableNames',VarNames);

clc;
disp(['PV correlation as a function of elapsed time'])
disp(['Friedman’s tests with Holm–Bonferroni correction:'])
disp(statistics)

%% Figure 2F - Activity rates of three PM units across movie repeats within a block
nat_movie = 1;
mouse = 53;
area = 4;
example_mouse = squeeze(nanmean(neuropixels_population_vectors{mouse,area,nat_movie}(:,:,1:30)*30,2));

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

nat_movie = 1;
mouse = 53;
area = 1;
example_mouse = neuropixels_population_vectors{mouse,area,nat_movie};
sub = 1;

for unit = [32,26,20]
    current_unit = squeeze(example_mouse(unit,:,:)*30)';
    smooth_current_unit = [];
    for repeat = 1:size(current_unit,1)
        smooth_current_unit(repeat,:) = imgaussfilt(current_unit(repeat,:),2);
    end
    norm_current_unit = smooth_current_unit ./ max(smooth_current_unit,[],2);
    [row,col] = find(norm_current_unit == 1);
    [B,I] = sort(row);
    subplot(2,3,sub)
    imagesc(norm_current_unit)
    hold on
    
    plot(xlim, [30 30]+0.5,'--','linewidth',2,'color','w')
    if sub ==1
        ylabel('Movie repeat')
        text(0.075, 0.9,['Block A'],'Units','normalized','color','w')
        text(0.075, 0.4,['Block B'],'Units','normalized','color','w')
    elseif sub == 3
        cb = colorbar;
        set(cb,'position',[0.925 0.575 0.04 0.35])
        
    end
    
    colormap(newmap3)
    title(['Unit #',num2str(unit)])
    
    subplot(2,3,sub+3)
    hold on
    mean_current_unit = [];
    mean_current_unit(1,:) = nanmean(current_unit(1:30,:));
    mean_current_unit(2,:) = nanmean(current_unit(31:60,:));
    plot(mean_current_unit(1,:),'-','markersize',20,'linewidth',1.5,'color',[0.3 0.3 0.3])
    plot(mean_current_unit(2,:)','-','markersize',20,'linewidth',1.5,'color',[0.6 0.6 0.6])
    [r,p] = corr(mean_current_unit(1,:)',mean_current_unit(2,:)');
    text(0.1, 0.95,['r = ',num2str(r)],'Units','normalized','color',[0.25 0.25 0.25])
    %text(0.1, 0.85,['p = ',num2str(p)],'Units','normalized','color',[0.25 0.25 0.25])
    xlim([1 30])
    if sub ==1
        ylabel('Mean activity rate (spike/sec)')
        legend({'Block A','Block B'})
        legend('boxoff')
        ylim([0 45])
    elseif sub ==2
        xlabel('Time in movie (sec)')
        ylim([0 60])
    else
        ylim([0 24])
    end
    sub = sub + 1;
end

%% Figure 2H - Ensemble rate correlation as function of elapsed time across mice
nat_movie = 1;
cell_cutoff = 15;
num_repeats = 30;
mean_rate_elapsed_repeat_per_area = {};
for area = 1:6
    valid_mice = neuropixels_cell_count(:,area,nat_movie) >= cell_cutoff & movie_repeats(:,nat_movie) == num_repeats;
    current_area = neuropixels_population_vectors(valid_mice,area,nat_movie);
    mean_rate_elapsed_repeat = [];
    
    for mouse = 1:length(current_area)
        clc;
        disp(['Calculating ensemble rate correlation between movie repeats for Neuropixels data:'])
        disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        current_mouse = current_area{mouse};
        current_mouse_blockA =  squeeze(nanmean(current_mouse(:,:,1:30),2));
        current_mouse_blockB =  squeeze(nanmean(current_mouse(:,:,31:60),2));
        rate_corr = [];
        rate_corr(:,:,1) = corr(current_mouse_blockA);
        rate_corr(:,:,2) = corr(current_mouse_blockB);
        
        mean_rate_corr_across_blocks = nanmean(rate_corr,3);
        
        for diagonal = 1:num_repeats-1
            mean_rate_elapsed_repeat(mouse,diagonal) = nanmean(diag(mean_rate_corr_across_blocks,diagonal));
        end
    end
    mean_rate_elapsed_repeat_per_area{area} = mean_rate_elapsed_repeat;
end

pval = [];
df = [];
chi = [];
figure
for area = 1:6
    current_area = mean_rate_elapsed_repeat_per_area{area};
    [pval(area),tbl,stats] = friedman(current_area ,[1],'off');
    df(area) = tbl{2,3};
    chi(area) = tbl{2,5};
    
    mean_across_mice = nanmean(current_area);
    std_across_mice = nanstd(current_area);
    ste_across_mice = std_across_mice./sqrt(size(current_area,1));
    
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

corrected_pval = bonf_holm(pval);

VarNames = {'area','df','chi','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas(1:6)'),df(:),chi(:),pval(:),corrected_pval(:),'VariableNames',VarNames);

clc;
disp(['Ensemble rate correlation as a function of elapsed time'])
disp(['Friedman’s tests with Holm–Bonferroni correction:'])
disp(statistics)

%% Figure 2I - Tuning curve correlation as function of elapsed time across mice
nat_movie = 1;
cell_cutoff = 15;
num_repeats = 30;
mean_tuning_elapsed_repeat_per_area = {};
for area = 1:6
    valid_mice = neuropixels_cell_count(:,area,nat_movie) >= cell_cutoff & movie_repeats(:,nat_movie) == num_repeats;
    current_area = neuropixels_population_vectors(valid_mice,area,nat_movie);
    mean_tuning_elapsed_repeat = [];
    
    for mouse = 1:length(current_area)
        clc;
        disp(['Calculating tuning curve correlation between movie repeats for Neuropixels data:'])
        disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        current_mouse = current_area{mouse};
        current_mouse_blockA =  current_mouse(:,:,1:30);
        current_mouse_blockB =  current_mouse(:,:,31:60);
        
        current_mouse_tuning_corr = [];
        for repeat1 = 1:30
            for repeat2 = 1:30
                current_tuning_blockA = corr(current_mouse_blockA(:,:,repeat1)',current_mouse_blockA(:,:,repeat2)');
                current_mouse_tuning_corr(repeat1,repeat2,1) = nanmedian(diag(current_tuning_blockA));
                
                current_tuning_blockB = corr(current_mouse_blockB(:,:,repeat1)',current_mouse_blockB(:,:,repeat2)');
                current_mouse_tuning_corr(repeat1,repeat2,2) = nanmedian(diag(current_tuning_blockB));
            end
        end
        mean_tuning_corr_across_blocks = nanmean(current_mouse_tuning_corr,3);
        
        for diagonal = 1:num_repeats-1
            mean_tuning_elapsed_repeat(mouse,diagonal) = nanmean(diag(mean_tuning_corr_across_blocks,diagonal));
        end
    end
    mean_tuning_elapsed_repeat_per_area{area} = mean_tuning_elapsed_repeat;
end

pval = [];
df = [];
chi = [];
figure
for area = 1:6
    current_area = mean_tuning_elapsed_repeat_per_area{area};
    [pval(area),tbl,stats] = friedman(current_area,1,'off');
    df(area) = tbl{2,3};
    chi(area) = tbl{2,5};
    
    mean_across_mice = nanmean(current_area);
    std_across_mice = nanstd(current_area);
    ste_across_mice = std_across_mice./sqrt(size(current_area,1));
    
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

corrected_pval = bonf_holm(pval);

VarNames = {'area','df','chi','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas(1:6)'),df(:),chi(:),pval(:),corrected_pval(:),'VariableNames',VarNames);

clc;
disp(['Tuning curve correlation as a function of elapsed time'])
disp(['Friedman’s tests with Holm–Bonferroni correction:'])
disp(statistics)

%% Figure 3A - Between blocks stability - PV correlation for area LM across mice
nat_movie = 2;
area = 2;
num_repeats = 5;
cell_cutoff = 15;
valid_mice = neuropixels_cell_count(:,area,nat_movie) >= cell_cutoff & movie_repeats(:,nat_movie) == num_repeats;
current_area = neuropixels_population_vectors(valid_mice,area,nat_movie);

pv_corr_across_mice = [];
for mouse = 1:length(current_area)
    current_mouse = current_area{mouse};
    
    
    current_mouse_blockA =  current_mouse(:,:,1:5);
    current_mouse_blockA_half1 = nanmean(current_mouse_blockA(:,:,1:2),3);
    current_mouse_blockA_half2 = nanmean(current_mouse_blockA(:,:,3:5),3);
    
    current_mouse_blockB =  current_mouse(:,:,6:10);
    current_mouse_blockB_half1 = nanmean(current_mouse_blockB(:,:,1:2),3);
    current_mouse_blockB_half2 = nanmean(current_mouse_blockB(:,:,3:5),3);
    
    pv_corr_across_mice(:,:,mouse) = corr([current_mouse_blockA_half1,current_mouse_blockA_half2,...
        current_mouse_blockB_half1,current_mouse_blockB_half2]);
    
end

mean_pv_corr_across_mice = nanmean(pv_corr_across_mice,3);

figure
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
nat_movie = 2;
num_repeats = 5;
cell_cutoff = 15;
pv_corr_between_blocks_across_areas = {};
for area = 1:6
    valid_mice = neuropixels_cell_count(:,area,nat_movie) >= cell_cutoff & movie_repeats(:,nat_movie) == num_repeats;
    current_area = neuropixels_population_vectors(valid_mice,area,nat_movie);
    
    within_between_stability = [];
    for mouse = 1:length(current_area)
        clc;
        disp(['Calculating PV correlation between blocks for Neuropixels data:'])
        disp(['Stimulus: Natural movie 3 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse};
        
        current_mouse_blockA =  current_mouse(:,:,1:5);
        current_mouse_blockB =  current_mouse(:,:,6:10);
        
        current_mouse_half_blocks = [];
        current_mouse_half_blocks(:,:,1) = nanmean(current_mouse_blockA(:,:,1:2),3);
        current_mouse_half_blocks(:,:,2) = nanmean(current_mouse_blockA(:,:,3:5),3);
        current_mouse_half_blocks(:,:,3) = nanmean(current_mouse_blockB(:,:,1:2),3);
        current_mouse_half_blocks(:,:,4) = nanmean(current_mouse_blockB(:,:,3:5),3);
        
        
        mean_pv_corr_across_mice = [];
        for half1 = 1:4
            for half2 = 1:4
                current_pv_corr = corr(current_mouse_half_blocks(:,:,half1),current_mouse_half_blocks(:,:,half2));
                mean_pv_corr_across_mice(half1,half2) = nanmean(diag(current_pv_corr));
            end
        end
        
        within_between_stability(mouse,1) = nanmean([mean_pv_corr_across_mice(1,2),mean_pv_corr_across_mice(3,4)]);
        within_between_stability(mouse,2) = nanmean(nanmean(mean_pv_corr_across_mice(1:2,3:4)));
        
        
    end
    pv_corr_between_blocks_across_areas{area} = within_between_stability;
end

pvalue = [];
zvalue = [];
figure('Units','Normalized','Position',[0.2 0.4 0.3 0.325])
for area = 1:6
    current_area = pv_corr_between_blocks_across_areas{area};
    mean_stability = nanmean(current_area);
    std_stability = nanstd(current_area);
    ste_stability = std_stability./sqrt(size(current_area,1));
    
    [pvalue(area),~,stats] = signrank(current_area(:,1),current_area(:,2));
    zvalue(area) = stats.zval;
    
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

corrected_pvalue = bonf_holm(pvalue);

VarNames = {'area','zvalue','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas(1:6)'),zvalue(:),pvalue(:),corrected_pvalue(:),'VariableNames',VarNames);

clc;
disp(['PV correlation within blocks compared to between blocks'])
disp(['Two-tailed Wilcoxon signed-rank tests with Holm–Bonferroni correction:'])
disp(statistics)

%% Figure 3C - Ensemble rate correlation within and between blocks of natural movie 3
nat_movie = 2;
num_repeats = 5;
cell_cutoff = 15;
pv_corr_between_blocks_across_areas = {};
for area = 1:6
    valid_mice = neuropixels_cell_count(:,area,nat_movie) >= cell_cutoff & movie_repeats(:,nat_movie) == num_repeats;
    current_area = neuropixels_population_vectors(valid_mice,area,nat_movie);
    
    within_between_stability = [];
    for mouse = 1:length(current_area)
        clc;
        disp(['Calculating ensemble rate correlation between blocks for Neuropixels data:'])
        disp(['Stimulus: Natural movie 3 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse};
        
        current_mouse_blockA =  current_mouse(:,:,1:5);
        current_mouse_blockB =  current_mouse(:,:,6:10);
        
        current_mouse_half_blocks = [];
        current_mouse_half_blocks(:,1) = nanmean(nanmean(current_mouse_blockA(:,:,1:2),3),2);
        current_mouse_half_blocks(:,2) = nanmean(nanmean(current_mouse_blockA(:,:,3:5),3),2);
        current_mouse_half_blocks(:,3) = nanmean(nanmean(current_mouse_blockB(:,:,1:2),3),2);
        current_mouse_half_blocks(:,4) = nanmean(nanmean(current_mouse_blockB(:,:,3:5),3),2);
        
        rate_corr = corr(current_mouse_half_blocks);
        
        within_between_stability(mouse,1) = nanmean([rate_corr(1,2),rate_corr(3,4)]);
        within_between_stability(mouse,2) = nanmean(nanmean(rate_corr(1:2,3:4)));
        
        
    end
    rate_corr_between_blocks_across_areas{area} = within_between_stability;
end

pvalue = [];
zvalue = [];
plt = [];
figure('Units','Normalized','Position',[0.2 0.4 0.2 0.325])
for area = 1:6
    current_area = rate_corr_between_blocks_across_areas{area};
    mean_stability = nanmean(current_area);
    std_stability = nanstd(current_area);
    ste_stability = std_stability./sqrt(size(current_area,1));
    
    [pvalue(area),~,stats] = signrank(current_area(:,1),current_area(:,2));
    zvalue(area) = stats.zval;
    
    
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

corrected_pvalue = bonf_holm(pvalue);

VarNames = {'area','zvalue','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas(1:6)'),zvalue(:),pvalue(:),corrected_pvalue(:),'VariableNames',VarNames);

clc;
disp(['Ensemble rate correlation within blocks compared to between blocks'])
disp(['Two-tailed Wilcoxon signed-rank tests with Holm–Bonferroni correction:'])
disp(statistics)

%% Figure 3D - Tuning curve correlation within and between blocks of natural movie 3
nat_movie = 2;
num_repeats = 5;
cell_cutoff = 15;
tuning_corr_between_blocks_across_areas = {};
for area = 1:6
    valid_mice = neuropixels_cell_count(:,area,nat_movie) >= cell_cutoff & movie_repeats(:,nat_movie) == num_repeats;
    current_area = neuropixels_population_vectors(valid_mice,area,nat_movie);
    
    within_between_stability = [];
    for mouse = 1:length(current_area)
        clc;
        disp(['Calculating tuning curve correlation between blocks for Neuropixels data:'])
        disp(['Stimulus: Natural movie 3 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse};
        
        current_mouse_blockA =  current_mouse(:,:,1:5);
        current_mouse_blockB =  current_mouse(:,:,6:10);
        
        current_mouse_half_blocks = [];
        current_mouse_half_blocks(:,:,1) = nanmean(current_mouse_blockA(:,:,1:2),3);
        current_mouse_half_blocks(:,:,2) = nanmean(current_mouse_blockA(:,:,3:5),3);
        current_mouse_half_blocks(:,:,3) = nanmean(current_mouse_blockB(:,:,1:2),3);
        current_mouse_half_blocks(:,:,4) = nanmean(current_mouse_blockB(:,:,3:5),3);
        
        
        tuning_corr_across_blocks = [];
        for half1 = 1:4
            for half2 = 1:4
                current_tuning_corr = corr(current_mouse_half_blocks(:,:,half1)',current_mouse_half_blocks(:,:,half2)');
                tuning_corr_across_blocks(half1,half2) = nanmedian(diag(current_tuning_corr));
            end
        end
        
        within_between_stability(mouse,1) = nanmean([tuning_corr_across_blocks(1,2),tuning_corr_across_blocks(3,4)]);
        within_between_stability(mouse,2) = nanmean(nanmean(tuning_corr_across_blocks(1:2,3:4)));
        
        
    end
    tuning_corr_between_blocks_across_areas{area} = within_between_stability;
end

pvalue = [];
zvalue = [];
plt = [];
figure('Units','Normalized','Position',[0.2 0.4 0.2 0.325])
for area = 1:6
    current_area = tuning_corr_between_blocks_across_areas{area};
    mean_stability = nanmean(current_area);
    std_stability = nanstd(current_area);
    ste_stability = std_stability./sqrt(size(current_area,1));
    
    [pvalue(area),~,stats] = signrank(current_area(:,1),current_area(:,2));
    zvalue(area) = stats.zval;
    
    
    hold on
    plt(area) = errorbar(mean_stability,ste_stability,'o','color',colors(area,:),...
        'markerfacecolor',colors(area,:),'capsize',0,'linestyle','-','linewidth',3,'markersize',5);
    
end
xlim([0.5 2.5])
ylim([0.35 0.7])
lgd = legend(plt,brain_areas(1:6));
legend('boxoff')
lgd.Position = [0.75 0.7 0.15 0.15];
set(gca,'xtick',[1,2],'xticklabel',{'Within block','Between blocks'},'ytick',0.4:0.1:0.1)
ylabel('Tuning curve correlation')
corrected_pvalue = bonf_holm(pvalue);

VarNames = {'area','zvalue','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas(1:6)'),zvalue(:),pvalue(:),corrected_pvalue(:),'VariableNames',VarNames);

clc;
disp(['Tuning curve correlation within blocks compared to between blocks'])
disp(['Two-tailed Wilcoxon signed-rank tests with Holm–Bonferroni correction:'])
disp(statistics)

%% Figure 3E - PV correlation across sessions for a single V1 mouse
area = 1;
nat_movie = 1;
mouse = 3;

current_mouse = calcium_excitatory_population_vectors{area}{mouse,nat_movie};
mean_current_mouse_sessions = [];
mean_current_mouse_sessions(:,:,1) =  nanmean(current_mouse(:,:,1:10),3);
mean_current_mouse_sessions(:,:,2)  =  nanmean(current_mouse(:,:,11:20),3);
mean_current_mouse_sessions(:,:,3)  =  nanmean(current_mouse(:,:,21:30),3);

pv_corr_pairwise_strict = [];
for sessionA = 1:3
    sessionA_activity = mean_current_mouse_sessions(:,:,sessionA);
    rows = [1:30] +30*(sessionA-1);
    for sessionB = 1:3
        cols = [1:30] +30*(sessionB-1);
        sessionB_activity = mean_current_mouse_sessions(:,:,sessionB);
        
        valid_cells = [nanmean(sessionA_activity,2) ~= 0] & [nanmean(sessionB_activity,2) ~= 0];
        valid_sessionA_activity = sessionA_activity(valid_cells,:);
        valid_sessionB_activity = sessionB_activity(valid_cells,:);
        
        pv_corr_pairwise_strict(rows,cols) = corr(valid_sessionA_activity,valid_sessionB_activity);
        
    end
end

figure('unit','normalized','position',[0.3 0.3 0.215 0.3])
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
nat_movie = 1;
cell_cutoff = 20;
within_between_session_stability_areas = {};
for area = 1:6
    valid_mice = min(calcium_excitatory_cell_count{area},[],2)>= cell_cutoff;
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie);
    within_between_session_stability = [];
    for mouse = 1:length(current_area)
        clc;
        disp(['Calculating PV correlation between sessions for Calcium imaging data:'])
        disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse};
        mean_activity_per_half = [];
        for half = 1:6
            current_half = [1:5] + 5*(half-1);
            mean_activity_per_half(:,:,half) = nanmean(current_mouse(:,:,current_half),3);
        end
        
        pv_corr_pairwise_strict = [];
        for halfA = 1:size(mean_activity_per_half,3)
            halfA_activity = mean_activity_per_half(:,:,halfA);
            for halfB = 1:size(mean_activity_per_half,3)
                halfB_activity = mean_activity_per_half(:,:,halfB);
                
                valid_cells = [nanmean(halfA_activity,2) ~= 0] & [nanmean(halfB_activity,2) ~= 0];
                valid_halfA_activity = halfA_activity(valid_cells,:);
                valid_halfB_activity = halfB_activity(valid_cells,:);
                
                pv_corr_pairwise_strict(halfA,halfB) = nanmean(diag(corr(valid_halfA_activity,valid_halfB_activity)));
            end
        end
        
        within_between_session_stability(mouse,1) = nanmean([pv_corr_pairwise_strict(1,2),pv_corr_pairwise_strict(3,4),pv_corr_pairwise_strict(5,6)]);
        within_between_session_stability(mouse,2) = nanmean([nanmean(nanmean(pv_corr_pairwise_strict(1:2,3:4))),...
            nanmean(nanmean(pv_corr_pairwise_strict(3:4,5:6)))]);
        within_between_session_stability(mouse,3) = nanmean(nanmean(pv_corr_pairwise_strict(1:2,5:6)));
    end
    within_between_session_stability_areas{area} = within_between_session_stability;
    
end


pvalue = [];
zvalue = [];
ylims = [0.375,0.6; 0.425,0.65; 0.4,0.6; 0.35 0.55;0.2,0.35;0.3,0.55];
figure('Units','Normalized','Position',[0.2 0.4 0.3 0.325])
for area = 1:6
    current_area = within_between_session_stability_areas{area};
    within_vs_between_sess = [current_area(:,1),nanmean(current_area(:,2:3),2)];
    
    mean_stability = nanmean(within_vs_between_sess);
    std_stability = nanstd(within_vs_between_sess);
    ste_stability = std_stability./sqrt(size(within_vs_between_sess,1));
    
    [pvalue(area),~,stats] = signrank(within_vs_between_sess(:,1),within_vs_between_sess(:,2));
    zvalue(area) = stats.zval;
    
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

corrected_pvalue = bonf_holm(pvalue);

VarNames = {'area','zvalue','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas(1:6)'),zvalue(:),pvalue(:),corrected_pvalue(:),'VariableNames',VarNames);

clc;
disp(['PV correlation within sessions compared to between sessions'])
disp(['Two-tailed Wilcoxon signed-rank tests with Holm–Bonferroni correction:'])
disp(statistics)


%% Figure 3G - Ensemble rate correlation within session and between sessions
nat_movie = 1;
cell_cutoff = 20;
within_between_session_stability_areas = {};
for area = 1:6
    valid_mice = min(calcium_excitatory_cell_count{area},[],2)>= cell_cutoff;
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie);
    within_between_session_stability = [];
    for mouse = 1:length(current_area)
        clc;
        disp(['Calculating ensemble rate correlation between sessions for Calcium imaging data:'])
        disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse};
        mean_activity_per_half = [];
        for half = 1:6
            current_half = [1:5] + 5*(half-1);
            mean_activity_per_half(:,:,half) = nanmean(current_mouse(:,:,current_half),3);
        end
        
        rate_corr_pairwise_strict = [];
        for halfA = 1:size(mean_activity_per_half,3)
            halfA_activity = mean_activity_per_half(:,:,halfA);
            for halfB = 1:size(mean_activity_per_half,3)
                halfB_activity = mean_activity_per_half(:,:,halfB);
                
                valid_cells = [nanmean(halfA_activity,2) ~= 0] & [nanmean(halfB_activity,2) ~= 0];
                valid_halfA_activity = nanmean(halfA_activity(valid_cells,:),2);
                valid_halfB_activity = nanmean(halfB_activity(valid_cells,:),2);
                
                rate_corr_pairwise_strict(halfA,halfB) = corr(valid_halfA_activity,valid_halfB_activity);
            end
        end
        
        within_between_session_stability(mouse,1) = nanmean([rate_corr_pairwise_strict(1,2),rate_corr_pairwise_strict(3,4),rate_corr_pairwise_strict(5,6)]);
        within_between_session_stability(mouse,2) = nanmean([nanmean(nanmean(rate_corr_pairwise_strict(1:2,3:4))),...
            nanmean(nanmean(rate_corr_pairwise_strict(3:4,5:6)))]);
        within_between_session_stability(mouse,3) = nanmean(nanmean(rate_corr_pairwise_strict(1:2,5:6)));
    end
    within_between_session_stability_areas{area} = within_between_session_stability;
    
end


pvalue = [];
zvalue = [];
plt = [];
figure('Units','Normalized','Position',[0.2 0.4 0.2 0.325])
for area = 1:6
    current_area = within_between_session_stability_areas{area};
    within_vs_between = [current_area(:,1),nanmean(current_area(:,2:3),2)];
    mean_stability = nanmean(within_vs_between);
    std_stability = nanstd(within_vs_between);
    ste_stability = std_stability./sqrt(size(within_vs_between,1));
    
    [pvalue(area),~,stats] = signrank(within_vs_between(:,1),within_vs_between(:,2));
    zvalue(area) = stats.zval;
    
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

corrected_pvalue = bonf_holm(pvalue);

VarNames = {'area','zvalue','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas(1:6)'),zvalue(:),pvalue(:),corrected_pvalue(:),'VariableNames',VarNames);

clc;
disp(['Ensemble rate correlation within sessions compared to between sessions'])
disp(['Two-tailed Wilcoxon signed-rank tests with Holm–Bonferroni correction:'])
disp(statistics)


%% Figure 3H - Tuning curve correlation within session and between sessions
nat_movie = 1;
cell_cutoff = 20;
within_between_session_stability_areas = {};
for area = 1:6
    valid_mice = min(calcium_excitatory_cell_count{area},[],2)>= cell_cutoff;
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie);
    within_between_session_stability = [];
    for mouse = 1:length(current_area)
        clc;
        disp(['Calculating tuning curve correlation between sessions for Calcium imaging data:'])
        disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse};
        mean_activity_per_half = [];
        for half = 1:6
            current_half = [1:5] + 5*(half-1);
            mean_activity_per_half(:,:,half) = nanmean(current_mouse(:,:,current_half),3);
        end
        
        tuning_corr_pairwise_strict = [];
        for halfA = 1:size(mean_activity_per_half,3)
            halfA_activity = mean_activity_per_half(:,:,halfA);
            for halfB = 1:size(mean_activity_per_half,3)
                halfB_activity = mean_activity_per_half(:,:,halfB);
                
                valid_cells = [nanmean(halfA_activity,2) ~= 0] & [nanmean(halfB_activity,2) ~= 0];
                valid_halfA_activity = halfA_activity(valid_cells,:)';
                valid_halfB_activity = halfB_activity(valid_cells,:)';
                
                tuning_corr_pairwise_strict(halfA,halfB) = nanmedian(diag(corr(valid_halfA_activity,valid_halfB_activity)));
            end
        end
        
        within_between_session_stability(mouse,1) = nanmean([tuning_corr_pairwise_strict(1,2),tuning_corr_pairwise_strict(3,4),tuning_corr_pairwise_strict(5,6)]);
        within_between_session_stability(mouse,2) = nanmean([nanmean(nanmean(tuning_corr_pairwise_strict(1:2,3:4))),...
            nanmean(nanmean(tuning_corr_pairwise_strict(3:4,5:6)))]);
        within_between_session_stability(mouse,3) = nanmean(nanmean(tuning_corr_pairwise_strict(1:2,5:6)));
    end
    within_between_session_stability_areas{area} = within_between_session_stability;
    
end


pvalue = [];
zvalue = [];
plt = [];
figure('Units','Normalized','Position',[0.2 0.4 0.2 0.325])
for area = 1:6
    current_area = within_between_session_stability_areas{area};
    within_vs_between = [current_area(:,1),nanmean(current_area(:,2:3),2)];
    mean_stability = nanmean(within_vs_between);
    std_stability = nanstd(within_vs_between);
    ste_stability = std_stability./sqrt(size(within_vs_between,1));
    
    [pvalue(area),~,stats] = signrank(within_vs_between(:,1),within_vs_between(:,2));
    zvalue(area) = stats.zval;
    
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

corrected_pvalue = bonf_holm(pvalue);

VarNames = {'area','zvalue','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas(1:6)'),zvalue(:),pvalue(:),corrected_pvalue(:),'VariableNames',VarNames);

clc;
disp(['Tuning curve correlation within sessions compared to between sessions'])
disp(['Two-tailed Wilcoxon signed-rank tests with Holm–Bonferroni correction:'])
disp(statistics)


%% Figure 3I - PV correlation between proximal and distal sessions
nat_movie = 1;
cell_cutoff = 20;
proximal_distal_session_stability_areas = {};
for area = 1:6
    valid_mice = min(calcium_excitatory_cell_count{area},[],2)>= cell_cutoff;
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie);
    proximal_distal_session_stability = [];
    for mouse = 1:length(current_area)
        clc;
        disp(['Calculating PV correlation between proximal and distal sessions for Calcium imaging data:'])
        disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse};
        mean_activity_per_half = [];
        for half = 1:3
            current_half = [1:10] + 10*(half-1);
            mean_activity_per_half(:,:,half) = nanmean(current_mouse(:,:,current_half),3);
        end
        
        pv_corr_pairwise_strict = [];
        for halfA = 1:size(mean_activity_per_half,3)
            halfA_activity = mean_activity_per_half(:,:,halfA);
            for halfB = 1:size(mean_activity_per_half,3)
                halfB_activity = mean_activity_per_half(:,:,halfB);
                
                valid_cells = [nanmean(halfA_activity,2) ~= 0] & [nanmean(halfB_activity,2) ~= 0];
                valid_halfA_activity = halfA_activity(valid_cells,:);
                valid_halfB_activity = halfB_activity(valid_cells,:);
                
                pv_corr_pairwise_strict(halfA,halfB) = nanmean(diag(corr(valid_halfA_activity,valid_halfB_activity)));
            end
        end
        
        proximal_distal_session_stability(mouse,:) = [nanmean(diag(pv_corr_pairwise_strict,1)),diag(pv_corr_pairwise_strict,2)];
    end
    proximal_distal_session_stability_areas{area} = proximal_distal_session_stability;
end


pvalue = [];
zvalue = [];
plt = [];
mean_stability = [];
ste_stability = [];
for area = 1:6
    current_area = proximal_distal_session_stability_areas{area};
    proximal_vs_distal = [current_area(:,1)-current_area(:,2)];
    mean_stability(area) = nanmean(proximal_vs_distal);
    std_stability = nanstd(proximal_vs_distal);
    ste_stability(area) = std_stability./sqrt(length(proximal_vs_distal));
    
    [pvalue(area),~,stats] = signrank(current_area(:,1),current_area(:,2),'tail','right');
    zvalue(area) = stats.zval;
end



figure('Units','Normalized','Position',[0.2 0.4 0.185 0.325])
bar(mean_stability,'facecolor',[0.8 0.8 0.8],'edgecolor','none')
hold on
errorbar(mean_stability,ste_stability,'capsize',0,'linewidth',2,'color',[0.6 0.6 0.6],'linestyle','none','markerfacecolor',[0.6 0.6 0.6])
xlim([0 7])
set(gca,'xtick',[1:6],'xticklabel',brain_areas(1:6))
ylabel('PV correlation difference')
set(gca, 'box','off')

corrected_pvalue = bonf_holm(pvalue);
VarNames = {'area','zvalue','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas(1:6)'),zvalue(:),pvalue(:),corrected_pvalue(:),'VariableNames',VarNames);

clc;
disp(['PV correlation between proximal sessions compared to distal sessions'])
disp(['One-tailed Wilcoxon signed-rank tests with Holm–Bonferroni correction:'])
disp(statistics)


%% Figure 3J - Ensemble rate correlation as a function of elapsed days
nat_movie = 1;
cell_cutoff = 20;
rate_corr_across_sessions_areas = {};
elapsed_days_areas = {};
for area = 1:6
    valid_mice = min(calcium_excitatory_cell_count{area},[],2)>= cell_cutoff;
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie);
    current_area_age = calcium_excitatory_sorted_mouse_age{area}(valid_mice,:);
    rate_corr_across_sessions = [];
    elapsed_days = [];
    for mouse = 1:length(current_area)
        clc;
        disp(['Calculating ensemble rate correlation between proximal and distal sessions for Calcium imaging data:'])
        disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse};
        current_mouse_area = current_area_age(mouse,:);
        mean_activity_per_half = [];
        for half = 1:3
            current_half = [1:10] + 10*(half-1);
            mean_activity_per_half(:,:,half) = nanmean(current_mouse(:,:,current_half),3);
        end
        
        rate_corr_pairwise_strict = [];
        for halfA = 1:size(mean_activity_per_half,3)
            halfA_activity = mean_activity_per_half(:,:,halfA);
            for halfB = 1:size(mean_activity_per_half,3)
                halfB_activity = mean_activity_per_half(:,:,halfB);
                
                valid_cells = [nanmean(halfA_activity,2) ~= 0] & [nanmean(halfB_activity,2) ~= 0];
                valid_halfA_activity = nanmean(halfA_activity(valid_cells,:),2);
                valid_halfB_activity = nanmean(halfB_activity(valid_cells,:),2);
                
                rate_corr_pairwise_strict(halfA,halfB) = corr(valid_halfA_activity,valid_halfB_activity);
            end
        end
        
        rate_corr_across_sessions(mouse,:) = [diag(rate_corr_pairwise_strict,1)',diag(rate_corr_pairwise_strict,2)];
        elapsed_days(mouse,:) = [current_mouse_area(2)-current_mouse_area(1),...
            current_mouse_area(3)-current_mouse_area(2),...
            current_mouse_area(3)-current_mouse_area(1)];
    end
    rate_corr_across_sessions_areas{area} = rate_corr_across_sessions;
    elapsed_days_areas{area} = elapsed_days;
end

r = [];
pvals=[];
figure
for area = 1:6
    current_area = rate_corr_across_sessions_areas{area};
    current_elapsed_day = elapsed_days_areas{area};
    
    within_day_intervals = sum(current_elapsed_day==0,2)>0;
    current_area = current_area(~within_day_intervals,:);
    current_elapsed_day = current_elapsed_day(~within_day_intervals,:);
    
    
    mean_correlations = [];
    mean_elapsed_day = [];
    for mouse = 1:size(current_area,1)
        current_mouse = current_area(mouse,:);
        current_days = current_elapsed_day(mouse,:);
        if current_days(1) == current_days(2)
            mean_correlations = [mean_correlations,[nanmean(current_mouse(1:2)),current_mouse(3)]];
            mean_elapsed_day  = [mean_elapsed_day,[nanmean(current_days(1:2)),current_days(3)]];
        else
            mean_correlations = [mean_correlations,current_mouse];
            mean_elapsed_day  = [mean_elapsed_day,current_days];
        end
    end
    
    
    max_interval = max(mean_elapsed_day);
    
    [~,i] = sort(mean_elapsed_day);
    x = mean_elapsed_day(i)';
    y = mean_correlations(i)';
    [p,s] = polyfit(x,y,1);
    [yfit,dy] = polyconf(p,x,s,'predopt','curve');
    
    [r(area),pvals(area)] = corr(mean_elapsed_day',mean_correlations','tail','left');
    
    subplot(2,3,area)
    hold on
    scatter(mean_elapsed_day',mean_correlations',20,[0.7 0.7 0.7],'filled','markerfacealpha',0.4)
    h = fill([x;flipud(x)],[yfit-dy;flipud(yfit+dy)],[0.6 0.6 0.6],'linestyle','none');
    line(x,yfit,'color',[0.6 0.6 0.6],'linewidth',2)
    set(h,'facealpha',.25)
    
    ylim([-0.2 1.2])
    xlim([0 max_interval+1])
    
    
    if area == 4
        ylabel('Ensemble rate correlation')
    elseif area == 5
        xlabel('Days between recording sessions')
    end
    
end

corrected_pvalue = bonf_holm(pvals);
for area = 1:6
    hold on
    subplot(2,3,area)
    text(0.1, 1,[brain_areas{area}],'Units','normalized','color',[0 0 0],'FontSize',12)
    text(0.6, 1,['r = ',num2str(round(r(area),2))],'Units','normalized','color',[0 0 0])
    text(0.6, 0.925,['p = ',num2str(corrected_pvalue(area))],'Units','normalized','color',[0 0 0])
end

%% Figure 3K - Tuning curve correlation as a function of elapsed days
nat_movie = 1;
cell_cutoff = 20;
rate_corr_across_sessions_areas = {};
elapsed_days_areas = {};
for area = 1:6
    valid_mice = min(calcium_excitatory_cell_count{area},[],2)>= cell_cutoff;
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie);
    current_area_age = calcium_excitatory_sorted_mouse_age{area}(valid_mice,:);
    tuning_corr_across_sessions = [];
    elapsed_days = [];
    for mouse = 1:length(current_area)
        clc;
        disp(['Calculating tuning curve correlation between proximal and distal sessions for Calcium imaging data:'])
        disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse};
        current_mouse_area = current_area_age(mouse,:);
        mean_activity_per_half = [];
        for half = 1:3
            current_half = [1:10] + 10*(half-1);
            mean_activity_per_half(:,:,half) = nanmean(current_mouse(:,:,current_half),3);
        end
        
        tuning_corr_pairwise_strict = [];
        for halfA = 1:size(mean_activity_per_half,3)
            halfA_activity = mean_activity_per_half(:,:,halfA);
            for halfB = 1:size(mean_activity_per_half,3)
                halfB_activity = mean_activity_per_half(:,:,halfB);
                
                valid_cells = [nanmean(halfA_activity,2) ~= 0] & [nanmean(halfB_activity,2) ~= 0];
                valid_halfA_activity = halfA_activity(valid_cells,:)';
                valid_halfB_activity = halfB_activity(valid_cells,:)';
                
                tuning_corr_pairwise_strict(halfA,halfB) = nanmedian(diag(corr(valid_halfA_activity,valid_halfB_activity)));
            end
        end
        
        tuning_corr_across_sessions(mouse,:) = [diag(tuning_corr_pairwise_strict,1)',diag(tuning_corr_pairwise_strict,2)];
        elapsed_days(mouse,:) = [current_mouse_area(2)-current_mouse_area(1),...
            current_mouse_area(3)-current_mouse_area(2),...
            current_mouse_area(3)-current_mouse_area(1)];
    end
    tuning_corr_across_sessions_areas{area} = tuning_corr_across_sessions;
    elapsed_days_areas{area} = elapsed_days;
end


r = [];
pvals=[];
figure
for area = 1:6
    current_area = tuning_corr_across_sessions_areas{area};
    current_elapsed_day = elapsed_days_areas{area};
    
    within_day_intervals = sum(current_elapsed_day==0,2)>0;
    current_area = current_area(~within_day_intervals,:);
    current_elapsed_day = current_elapsed_day(~within_day_intervals,:);
    
    
    mean_correlations = [];
    mean_elapsed_day = [];
    for mouse = 1:size(current_area,1)
        current_mouse = current_area(mouse,:);
        current_days = current_elapsed_day(mouse,:);
        if current_days(1) == current_days(2)
            mean_correlations = [mean_correlations,[nanmean(current_mouse(1:2)),current_mouse(3)]];
            mean_elapsed_day  = [mean_elapsed_day,[nanmean(current_days(1:2)),current_days(3)]];
        else
            mean_correlations = [mean_correlations,current_mouse];
            mean_elapsed_day  = [mean_elapsed_day,current_days];
        end
    end
    
    
    max_interval = max(mean_elapsed_day);
    
    [~,i] = sort(mean_elapsed_day);
    x = mean_elapsed_day(i)';
    y = mean_correlations(i)';
    [p,s] = polyfit(x,y,1);
    [yfit,dy] = polyconf(p,x,s,'predopt','curve');
    
    [r(area),pvals(area)] = corr(mean_elapsed_day',mean_correlations','tail','left');
    
    subplot(2,3,area)
    hold on
    scatter(mean_elapsed_day',mean_correlations',20,[0.7 0.7 0.7],'filled','markerfacealpha',0.4)
    h = fill([x;flipud(x)],[yfit-dy;flipud(yfit+dy)],[0.6 0.6 0.6],'linestyle','none');
    line(x,yfit,'color',[0.6 0.6 0.6],'linewidth',2)
    set(h,'facealpha',.25)
    
    ylim([-0.2 1.2])
    xlim([0 max_interval+1])
    
    
    if area == 4
        ylabel('Tuning curve correlation')
    elseif area == 5
        xlabel('Days between recording sessions')
    end
    
end

corrected_pvalue = bonf_holm(pvals);
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

depth_list = [150 250; 251 350; 351 500; 501 700];
layer_list = {'L2/3','L4','L5','L6'};
cell_cutoff = 20;
valid_mice_num = [];
for area = 1:6
    sub = 1;
    for depth = 1:size(depth_list)
        
        % valid mice
        valid_cell_num = min(calcium_excitatory_cell_count{area},[],2) >= cell_cutoff;
        valid_depth = calcium_excitatory_imaging_depth{area} >= depth_list(depth,1) & calcium_excitatory_imaging_depth{area} <= depth_list(depth,2);
        valid_mice = valid_cell_num & valid_depth';
        valid_mice_num(area,depth) = sum(valid_mice);
        
        elapsed_sess_pv = [];
        if valid_mice_num(area) >= 0
            current_area = calcium_excitatory_population_vectors{area}(valid_mice,1);
            
            for mouse = 1:size(current_area,1)
                clc;
                disp(['Calculating PV correlation across layers for Calcium imaging data:'])
                disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Depth: ',layer_list{depth},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
                
                current_mouse = current_area{mouse};
                
                pv_corr = [];
                for sess1 = 1:6
                    current_sess1 = nanmean(current_mouse(:,:,[1:5]+5*(sess1-1)),3);
                    for sess2 = 1:6
                        current_sess2 = nanmean(current_mouse(:,:,[1:5]+5*(sess2-1)),3);
                        valid_cells = nanmean(current_sess1,2) ~= 0 & nanmean(current_sess2,2) ~= 0;
                        
                        pv_corr(sess1,sess2) = nanmean(diag(corr(current_sess1(valid_cells,:),current_sess2(valid_cells,:))));
                    end
                end
                
                
                elapsed_sess_pv(mouse,1)= nanmean([pv_corr(1,2),pv_corr(3,4),pv_corr(5,6)]);
                elapsed_sess_pv(mouse,2)= nanmean(nanmean([pv_corr(1:2,3:4),pv_corr(3:4,5:6)]));
                elapsed_sess_pv(mouse,3)= nanmean(nanmean([pv_corr(1:2,5:6)]));
                
                
            end
            
        end
        elapsed_sess_pv_depth{area,depth} = elapsed_sess_pv;
    end
end


ylims = [0 0.7; 0.3 0.7; 0.3 0.7; 0.175 0.7; 0 0.5; 0.2 0.7];
figure('units','normalized','position',[0.3 0.3 0.4 0.4])
for area = 1:6
    for depth = 1:size(depth_list)
        curent_area_depth = elapsed_sess_pv_depth{area,depth};
        if ~isempty(curent_area_depth)
            
            mean_stability = nanmean(curent_area_depth);
            std_stability = nanstd(curent_area_depth);
            ste_stability=std_stability./sqrt(size(curent_area_depth,1));
            
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

depth_list = [150 250; 251 350; 351 500; 501 700];
layer_list = {'L2/3','L4','L5','L6'};
cell_cutoff = 20;
valid_mice_num = [];
for area = 1:6
    sub = 1;
    for depth = 1:size(depth_list)
        
        % valid mice
        valid_cell_num = min(calcium_excitatory_cell_count{area},[],2) >= cell_cutoff;
        valid_depth = calcium_excitatory_imaging_depth{area} >= depth_list(depth,1) & calcium_excitatory_imaging_depth{area} <= depth_list(depth,2);
        valid_mice = valid_cell_num & valid_depth';
        valid_mice_num(area,depth) = sum(valid_mice);
        
        elapsed_sess_pv = [];
        if valid_mice_num(area) >= 0
            current_area = calcium_excitatory_population_vectors{area}(valid_mice,1);
            
            for mouse = 1:size(current_area,1)
                clc;
                disp(['Calculating PV correlation across layers for Calcium imaging data:'])
                disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Depth: ',layer_list{depth},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
                
                current_mouse = current_area{mouse};
                
                pv_corr = [];
                for sess1 = 1:6
                    current_sess1 = nanmean(current_mouse(:,:,[1:5]+5*(sess1-1)),3);
                    for sess2 = 1:6
                        current_sess2 = nanmean(current_mouse(:,:,[1:5]+5*(sess2-1)),3);
                        valid_cells = nanmean(current_sess1,2) ~= 0 & nanmean(current_sess2,2) ~= 0;
                        
                        pv_corr(sess1,sess2) = nanmean(diag(corr(current_sess1(valid_cells,:),current_sess2(valid_cells,:))));
                    end
                end
                
                
                elapsed_sess_pv(mouse,1)= nanmean([pv_corr(1,2),pv_corr(3,4),pv_corr(5,6)]);
                elapsed_sess_pv(mouse,2)= nanmean(nanmean([pv_corr(1:2,3:4),pv_corr(3:4,5:6)]));
                elapsed_sess_pv(mouse,3)= nanmean(nanmean([pv_corr(1:2,5:6)]));
                
                
            end
            
        end
        elapsed_sess_pv_depth{area,depth} = elapsed_sess_pv;
    end
end


figure
for depth = 1:size(depth_list)
    curent_depth = cell2mat(elapsed_sess_pv_depth(:,depth));
    
    
    mean_stability = nanmean(curent_depth);
    std_stability = nanstd(curent_depth);
    ste_stability=std_stability./sqrt(size(curent_depth,1));
    
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

depth_list = [150 250; 251 350; 351 500; 501 700];
layer_list = {'L2/3','L4  ','L5  ','L6  '};
cell_cutoff = 20;
valid_mice_num = [];
for area = 1:6
    
    for depth = 1:size(depth_list)
        
        % valid mice
        valid_cell_num = min(calcium_excitatory_cell_count{area},[],2) >= cell_cutoff;
        valid_depth = calcium_excitatory_imaging_depth{area} >= depth_list(depth,1) & calcium_excitatory_imaging_depth{area} <= depth_list(depth,2);
        valid_mice = valid_cell_num & valid_depth';
        valid_mice_num(area,depth) = sum(valid_mice);
        
        elapsed_sess_pv = [];
        if valid_mice_num(area) >= 0
            current_area = calcium_excitatory_population_vectors{area}(valid_mice,1);
            
            for mouse = 1:size(current_area,1)
                clc;
                disp(['Calculating PV correlation across layers for Calcium imaging data:'])
                disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Depth: ',layer_list{depth},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
                
                current_mouse = current_area{mouse};
                
                pv_corr = [];
                for sess1 = 1:3
                    current_sess1 = nanmean(current_mouse(:,:,[1:10]+10*(sess1-1)),3);
                    for sess2 = 1:3
                        current_sess2 = nanmean(current_mouse(:,:,[1:10]+10*(sess2-1)),3);
                        valid_cells = nanmean(current_sess1,2) ~= 0 & nanmean(current_sess2,2) ~= 0;
                        
                        pv_corr(sess1,sess2) = nanmean(diag(corr(current_sess1(valid_cells,:),current_sess2(valid_cells,:))));
                    end
                end
                
                
                elapsed_sess_pv(mouse,:) = [nanmean(diag(pv_corr,1)),diag(pv_corr,2)];
                
            end
            
        end
        elapsed_sess_pv_depth{area,depth} = elapsed_sess_pv;
    end
end

pvalues = [];
mean_stability = [];
ste_stability = [];
pv_similarity_score_depth = {};
for depth = 1:4
    current_depth = cell2mat(elapsed_sess_pv_depth(:,depth));
    [pvalues(depth),~,stats] = signrank(current_depth(:,1),current_depth(:,2),'tail','right');
    
    current_depth(current_depth<0) = 0;
    pv_similarity_score = [current_depth(:,2)-current_depth(:,1)]./[current_depth(:,2)+current_depth(:,1)];
    mean_stability(depth) = nanmean(pv_similarity_score);
    std_stability = nanstd(pv_similarity_score);
    ste_stability(depth) = std_stability./sqrt(length(pv_similarity_score));
    
    pv_similarity_score_depth{depths} = pv_similarity_score;
end

figure
hold on
xlim([0.5 4.5])
for depth = 1:4
    errorbar(depth,mean_stability(depth),ste_stability(depth),'o','markerfacecolor',new_jet_colormap(depth+20*(depth-1),:),'color',new_jet_colormap(depth+20*(depth-1),:),'capsize',0,'linewidth',2,'linestyle','none')
end
plot(xlim,[0 0],'--','color',[0.4 0.4 0.4])
ylim([-0.2 0.05])
ylabel('Population vector similarity index')
set(gca,'xtick',[1:4],'xticklabel',{'L2/3','L4','L5','L6'})


corrected_pvalues = bonf_holm(pvalues);
VarNames = {'depth','pvalue','bonf_holm'};
statistics = table(cell2mat(layer_list'),pvalues(:),corrected_pvalues(:),'VariableNames',VarNames);

clc;
disp(['PV correlation between proximal sessions compared to distal sessions'])
disp(['One-tailed Wilcoxon signed-rank tests with Holm–Bonferroni correction:'])
disp(statistics)

%% Figure 4D - SST cells visual fields across movie repeats
nat_movie = 1;
mouse = 9;
area = 1;
example_mouse = calcium_inhibitory_population_vectors{area}{mouse,nat_movie}*30;
sub = 1;

figure('units','normalized','position',[0.3 0.3 0.3 0.5])
for unit = [5,4,7]
    current_unit = squeeze(example_mouse(unit,:,:))';
    smooth_current_unit = [];
    for repeat = 1:size(current_unit,1)
        smooth_current_unit(repeat,:) = imgaussfilt(current_unit(repeat,:),2);
    end
    norm_current_unit = smooth_current_unit ./ max(smooth_current_unit,[],2);
    [row,col] = find(norm_current_unit == 1);
    [B,I] = sort(row);
    subplot(2,3,sub)
    imagesc(norm_current_unit)
    hold on
    plot(xlim, [10 10]+0.5,'--','linewidth',2,'color','w')
    plot(xlim, [20 20]+0.5,'--','linewidth',2,'color','w')
    if sub ==1
        ylabel('Movie repeat')
        text(0.575, 0.925,['Sess 1'],'Units','normalized','color','w')
        text(0.575, 0.6,['Sess 2'],'Units','normalized','color','w')
        text(0.575, 0.275,['Sess 3'],'Units','normalized','color','w')
    elseif sub == 3
        cb = colorbar;
        set(cb,'position',[0.925 0.55 0.03 0.32])
        
    end
    colormap(newmap3)
    title(['Unit #',num2str(unit)])
    
    subplot(2,3,sub+3)
    hold on
    mean_current_unit = [];
    mean_current_unit(1,:) = nanmean(current_unit(1:10,:));
    mean_current_unit(2,:) = nanmean(current_unit(11:20,:));
    mean_current_unit(3,:) = nanmean(current_unit(21:30,:));
    plot(mean_current_unit(1,:),'-','markersize',20,'linewidth',1.5,'color',[0.3 0.3 0.3])
    plot(mean_current_unit(2,:)','-','markersize',20,'linewidth',1.5,'color',[0.5 0.5 0.5])
    plot(mean_current_unit(3,:)','-','markersize',20,'linewidth',1.5,'color',[0.7 0.7 0.7])
    
    
    xlim([1 30])
    if sub ==1
        ylabel('Mean activity rate (events/sec)')
        legend({'Sess 1','Sess 2','Sess 3'})
        legend('boxoff')
    elseif sub ==2
        xlabel('Time in movie (sec)')
    end
    sub = sub + 1;
end

suptitle('SST interneurons:')

%% Figure 4E - SST cells visual fields across movie repeats
nat_movie = 1;
mouse = 23;
area = 1;
example_mouse = calcium_inhibitory_population_vectors{area}{mouse,nat_movie}*30;
sub = 1;

figure('units','normalized','position',[0.3 0.3 0.3 0.5])
for unit = [4,10,16]
    current_unit = squeeze(example_mouse(unit,:,:))';
    smooth_current_unit = [];
    for repeat = 1:size(current_unit,1)
        smooth_current_unit(repeat,:) = imgaussfilt(current_unit(repeat,:),2);
    end
    norm_current_unit = smooth_current_unit ./ max(smooth_current_unit,[],2);
    [row,col] = find(norm_current_unit == 1);
    [B,I] = sort(row);
    subplot(2,3,sub)
    imagesc(norm_current_unit)
    hold on
    plot(xlim, [10 10]+0.5,'--','linewidth',2,'color','w')
    plot(xlim, [20 20]+0.5,'--','linewidth',2,'color','w')
    if sub ==1
        ylabel('Movie repeat')
        text(0.575, 0.925,['Sess 1'],'Units','normalized','color','w')
        text(0.575, 0.6,['Sess 2'],'Units','normalized','color','w')
        text(0.575, 0.275,['Sess 3'],'Units','normalized','color','w')
    elseif sub == 3
        cb = colorbar;
        set(cb,'position',[0.925 0.55 0.03 0.32])
        
    end
    colormap(newmap3)
    title(['Unit #',num2str(unit)])
    
    subplot(2,3,sub+3)
    hold on
    mean_current_unit = [];
    mean_current_unit(1,:) = nanmean(current_unit(1:10,:));
    mean_current_unit(2,:) = nanmean(current_unit(11:20,:));
    mean_current_unit(3,:) = nanmean(current_unit(21:30,:));
    plot(mean_current_unit(1,:),'-','markersize',20,'linewidth',1.5,'color',[0.3 0.3 0.3])
    plot(mean_current_unit(2,:)','-','markersize',20,'linewidth',1.5,'color',[0.5 0.5 0.5])
    plot(mean_current_unit(3,:)','-','markersize',20,'linewidth',1.5,'color',[0.7 0.7 0.7])
    
    
    xlim([1 30])
    if sub ==1
        ylabel('Mean activity rate (events/sec)')
        legend({'Sess 1','Sess 2','Sess 3'})
        legend('boxoff')
    elseif sub ==2
        xlabel('Time in movie (sec)')
    end
    sub = sub + 1;
end

suptitle('VIP interneurons:')

%% Figure 4F - Pvalb cells visual fields across movie repeats
nat_movie = 1;
mouse = 38;
area = 1;
example_mouse = calcium_inhibitory_population_vectors{area}{mouse,nat_movie}*30;
sub = 1;

figure('units','normalized','position',[0.3 0.3 0.3 0.5])
for unit = [1,15,16]
    current_unit = squeeze(example_mouse(unit,:,:))';
    smooth_current_unit = [];
    for repeat = 1:size(current_unit,1)
        smooth_current_unit(repeat,:) = imgaussfilt(current_unit(repeat,:),2);
    end
    norm_current_unit = smooth_current_unit ./ max(smooth_current_unit,[],2);
    [row,col] = find(norm_current_unit == 1);
    [B,I] = sort(row);
    subplot(2,3,sub)
    imagesc(norm_current_unit)
    hold on
    plot(xlim, [10 10]+0.5,'--','linewidth',2,'color','w')
    plot(xlim, [20 20]+0.5,'--','linewidth',2,'color','w')
    if sub ==1
        ylabel('Movie repeat')
        text(0.575, 0.925,['Sess 1'],'Units','normalized','color','w')
        text(0.575, 0.6,['Sess 2'],'Units','normalized','color','w')
        text(0.575, 0.275,['Sess 3'],'Units','normalized','color','w')
    elseif sub == 3
        cb = colorbar;
        set(cb,'position',[0.925 0.55 0.03 0.32])
        
    end
    colormap(newmap3)
    title(['Unit #',num2str(unit)])
    
    subplot(2,3,sub+3)
    hold on
    mean_current_unit = [];
    mean_current_unit(1,:) = nanmean(current_unit(1:10,:));
    mean_current_unit(2,:) = nanmean(current_unit(11:20,:));
    mean_current_unit(3,:) = nanmean(current_unit(21:30,:));
    plot(mean_current_unit(1,:),'-','markersize',20,'linewidth',1.5,'color',[0.3 0.3 0.3])
    plot(mean_current_unit(2,:)','-','markersize',20,'linewidth',1.5,'color',[0.5 0.5 0.5])
    plot(mean_current_unit(3,:)','-','markersize',20,'linewidth',1.5,'color',[0.7 0.7 0.7])
    
    
    xlim([1 30])
    if sub ==1
        ylabel('Mean activity rate (events/sec)')
        legend({'Sess 1','Sess 2','Sess 3'})
        legend('boxoff')
    elseif sub ==2
        xlabel('Time in movie (sec)')
    end
    sub = sub + 1;
end

suptitle('Pvalb interneurons:')

%% Figure 4G - Excitatory VS inhibitory - PV corr across movie repeats
nat_movie = 1;
num_repeats = 10;
cell_cutoff = 20;

cell_type = {'Excitatory','Inhibitory'};
area_list = [1,2,4];
elapse_repeat_pv_corr_areas = {};
for subtype = 1:2
    for area = 1:3
        if subtype == 1
            valid_mice = min(calcium_excitatory_cell_count{area_list(area)},[],2) >= cell_cutoff;
            current_area = calcium_excitatory_population_vectors{area_list(area)}(valid_mice,nat_movie);
        elseif subtype == 2
            current_area = calcium_inhibitory_population_vectors{area_list(area)}(:,nat_movie);
        end
        
        elapse_repeat_pv_corr = [];
        for mouse = 1:length(current_area)
            clc;
            disp(['Calculating PV correlation between movie repeats for calcium imaging data:'])
            disp(['Stimulus: Natural movie 1 | Cell type: ',cell_type{subtype},' | Area: ',brain_areas{area_list(area)},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
            
            current_mouse = current_area{mouse};
            
            current_mouse_sess1 =  current_mouse(:,:,1:10);
            valid_cells = nanmean(squeeze(nanmean(current_mouse_sess1,2)),2)>0;
            valid_current_mouse_sess1 = current_mouse_sess1(valid_cells,:,:);
            
            current_mouse_sess2 =  current_mouse(:,:,11:20);
            valid_cells = nanmean(squeeze(nanmean(current_mouse_sess2,2)),2)>0;
            valid_current_mouse_sess2 = current_mouse_sess2(valid_cells,:,:);
            
            current_mouse_sess3 =  current_mouse(:,:,21:30);
            valid_cells = nanmean(squeeze(nanmean(current_mouse_sess3,2)),2)>0;
            valid_current_mouse_sess3 = current_mouse_sess3(valid_cells,:,:);
            
            
            
            pv_corr_across_session = [];
            for repeat1 = 1:10
                for repeat2 = 1:10
                    pv_corr = corr(valid_current_mouse_sess1(:,:,repeat1),valid_current_mouse_sess1(:,:,repeat2));
                    pv_corr_across_session (repeat1,repeat2,1) = nanmean(diag(pv_corr));
                    
                    pv_corr = corr(valid_current_mouse_sess2(:,:,repeat1),valid_current_mouse_sess2(:,:,repeat2));
                    pv_corr_across_session(repeat1,repeat2,2) = nanmean(diag(pv_corr));
                    
                    pv_corr = corr(valid_current_mouse_sess3(:,:,repeat1),valid_current_mouse_sess3(:,:,repeat2));
                    pv_corr_across_session(repeat1,repeat2,3) = nanmean(diag(pv_corr));
                    
                end
            end
            mean_pv_corr_across_mice = nanmean(pv_corr_across_session,3);
            for diagonal = 1:9
                elapse_repeat_pv_corr(mouse,diagonal) = nanmean(diag(mean_pv_corr_across_mice,diagonal));
                
            end
            
        end
        elapse_repeat_pv_corr_areas{subtype,area} = elapse_repeat_pv_corr;
    end
end

pvalues = [];
df = [];
chi = [];
figure('units','normalized','position',[0.3 0.3 0.4 0.4])
for area = 1:3
    plt = [];
    for subtype = 1:2
        current_area = elapse_repeat_pv_corr_areas{subtype,area};
        
        [pvalues(subtype,area),tbl,stats] = friedman(current_area,[1],'off');
        df(subtype,area) = tbl{2,3};
        chi(subtype,area) = tbl{2,5};
        
        norm_elapse_repeat_pv_corr = current_area-current_area(:,1);
        
        subplot(1,3,area)
        hold on
        mean_norm_vals = nanmean(norm_elapse_repeat_pv_corr);
        std_norm_vals = nanstd(norm_elapse_repeat_pv_corr)./sqrt(mouse);
        x = [1:length(mean_norm_vals)]';
        y = mean_norm_vals';
        dy = std_norm_vals';
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

corrected_pval = bonf_holm(pvalues(2,:));

VarNames = {'area','df','chi','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas([1,2,4])'),df(2,:)',chi(2,:)',pvalues(2,:)',corrected_pval(:),'VariableNames',VarNames);

clc;
disp(['PV correlation as a function of elapsed time for inhibitory Cre lines'])
disp(['Friedman’s tests with Holm–Bonferroni correction:'])
disp(statistics)

%% Figure 4H - Excitatory VS inhibitory - PV corr across blocks
nat_movie = 2;
cell_cutoff = 20;
cell_type = {'Excitatory','Inhibitory'};
area_list = [1,2,4];
within_between_stability_area = {};
for subtype = 1:2
    for area = 1:3
        if subtype == 1
            valid_mice = calcium_excitatory_cell_count{area_list(area)}(:,1) >= cell_cutoff;
            current_area = calcium_excitatory_population_vectors{area_list(area)}(valid_mice,nat_movie);
        elseif subtype == 2
            current_area = calcium_inhibitory_population_vectors{area_list(area)}(:,nat_movie);
        end
        
        
        within_between_stability = [];
        for mouse = 1:length(current_area)
            clc;
            disp(['Calculating PV correlation between blocks for calcium imaging data:'])
            disp(['Stimulus: Natural movie 3 | Cell type: ',cell_type{subtype},' | Area: ',brain_areas{area_list(area)},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
            
            current_mouse = current_area{mouse};
            current_mouse_blocks = [];
            current_mouse_blocks(:,:,1) = nanmean(current_mouse(:,:,1:2),3);
            current_mouse_blocks(:,:,2) = nanmean(current_mouse(:,:,3:5),3);
            current_mouse_blocks(:,:,3) = nanmean(current_mouse(:,:,6:7),3);
            current_mouse_blocks(:,:,4) = nanmean(current_mouse(:,:,8:10),3);
            
            
            mean_pv_corr_across_mice = [];
            for half1 = 1:4
                for half2 = 1:4
                    current_pv = corr(current_mouse_blocks(:,:,half1),current_mouse_blocks(:,:,half2));
                    mean_pv_corr_across_mice(half1,half2) = nanmean(diag(current_pv));
                end
            end
            
            within_between_stability(mouse,1) = nanmean([mean_pv_corr_across_mice(1,2),mean_pv_corr_across_mice(3,4)]);
            within_between_stability(mouse,2) = nanmean(nanmean(mean_pv_corr_across_mice(1:2,3:4)));
            
        end
        within_between_stability_area{subtype,area} = within_between_stability;
    end
end

pvalues = [];
zvalues = [];
figure('units','normalized','position',[0.3 0.3 0.4 0.4])
for area = 1:3
    plt = [];
    for subtype = 1:2
        current_area = within_between_stability_area{subtype,area};
        
        [pvalues(subtype,area),~,stats] = signrank(current_area(:,1),current_area(:,2));
        zvalues(subtype,area) = stats.zval;
        
        mean_stability = nanmean(current_area);
        std_stability = nanstd(current_area);
        ste_stability = std_stability./sqrt(size(current_area,1));
        
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

corrected_pval = bonf_holm(pvalues(2,:));

VarNames = {'area','zvalue','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas([1,2,4])'),zvalues(2,:)',pvalues(2,:)',corrected_pval(:),'VariableNames',VarNames);

clc;
disp(['PV correlation within blocks compared to between blocks for inhibitory Cre lines'])
disp(['Two-tailed Wilcoxon signed-rank tests with Holm–Bonferroni correction:'])
disp(statistics)

%% Figure 4I - Excitatory VS inhibitory - PV corr across sessions
nat_movie = 1;
cell_cutoff = 20;
cell_type = {'Excitatory','Inhibitory'};
area_list = [1,2,4];
within_between_session_stability_areas = {};
for subtype = 1:2
    for area = 1:3
        if subtype == 1
            valid_mice = min(calcium_excitatory_cell_count{area_list(area)},[],2) >= cell_cutoff;
            current_area = calcium_excitatory_population_vectors{area_list(area)}(valid_mice,nat_movie);
        elseif subtype == 2
            current_area = calcium_inhibitory_population_vectors{area_list(area)}(:,nat_movie);
        end
        
        within_between_session_stability = [];
        for mouse = 1:length(current_area)
            clc;
            disp(['Calculating PV correlation between sessions for calcium imaging data:'])
            disp(['Stimulus: Natural movie 1 | Cell type: ',cell_type{subtype},' | Area: ',brain_areas{area_list(area)},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
            
            current_mouse = current_area{mouse};
            
            mean_activity_per_half = [];
            for half = 1:6
                current_half = [1:5] + 5*(half-1);
                mean_activity_per_half(:,:,half) = nanmean(current_mouse(:,:,current_half),3);
            end
            
            pv_corr = [];
            for halfA = 1:size(mean_activity_per_half,3)
                halfA_activity = mean_activity_per_half(:,:,halfA);
                for halfB = 1:size(mean_activity_per_half,3)
                    halfB_activity = mean_activity_per_half(:,:,halfB);
                    
                    valid_cells = [nanmean(halfA_activity,2) ~= 0] | [nanmean(halfB_activity,2) ~= 0];
                    valid_halfA_activity = halfA_activity(valid_cells,:);
                    valid_halfB_activity = halfB_activity(valid_cells,:);
                    
                    pv_corr(halfA,halfB) = nanmean(diag(corr(valid_halfA_activity,valid_halfB_activity)));
                end
            end
            
            within_between_session_stability(mouse,1) = nanmean([pv_corr(1,2),pv_corr(3,4),pv_corr(5,6)]);
            within_between_session_stability(mouse,2) = nanmean([nanmean(nanmean(pv_corr(1:2,3:4))),nanmean(nanmean(pv_corr(3:4,5:6)))]);
            within_between_session_stability(mouse,3) = nanmean(nanmean(pv_corr(1:2,5:6)));
        end
        within_between_session_stability_areas{subtype,area} = within_between_session_stability;
    end
end

pvalues = [];
zvalues = [];
figure('units','normalized','position',[0.3 0.3 0.4 0.4])
for area = 1:3
    plt = [];
    for subtype = 1:2
        current_area = within_between_session_stability_areas{subtype,area};
        
        [pvalues(subtype,area),~,stats] = signrank(current_area(:,2),current_area(:,3),'tail','right');
        zvalues(subtype,area) = stats.zval;
        
        mean_stability = nanmean(current_area);
        std_stability = nanstd(current_area);
        ste_stability = std_stability./sqrt(size(current_area,1));
        
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

corrected_pval = bonf_holm(pvalues(2,:));

VarNames = {'area','zvalue','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas([1,2,4])'),zvalues(2,:)',pvalues(2,:)',corrected_pval(:),'VariableNames',VarNames);

clc;
disp(['PV correlation  of proximal sessions compared to distal sessions for inhibitory Cre lines'])
disp(['One-tailed Wilcoxon signed-rank tests with Holm–Bonferroni correction:'])
disp(statistics)

%% Figure 5A - Visual hierarchy - ensemble rate correlation between movie repeats - Neuropixels
cell_cutoff = 15;
nat_movie = 1;
num_repeats = 30;
area_list = [1,2,7,8];

for area = 1:4
    valid_mice = (neuropixels_cell_count(:,area_list(area),nat_movie) >= cell_cutoff) & movie_repeats(:,nat_movie) == num_repeats;
    current_area = neuropixels_population_vectors(valid_mice,area_list(area),nat_movie);
    elapse_repeat_rate_corr = [];
    
    for mouse = 1:size(current_area,1)
        clc;
        disp(['Calculating ensemble rate correlation between repeats for Neuropixels data:'])
        disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area_list(area)},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse};
        current_mouse_blockA = squeeze(nanmean(current_mouse(:,:,1:30),2));
        current_mouse_blockB = squeeze(nanmean(current_mouse(:,:,31:60),2));
        
        mean_rate_corr = [];
        mean_rate_corr(:,:,1) = corr(current_mouse_blockA);
        mean_rate_corr(:,:,2) = corr(current_mouse_blockB);
        
        mean_rate_corr(mean_rate_corr<0) = 0;
        mean_rate_corr = nanmean(mean_rate_corr,3);
        
        for diagonal = 1:29
            elapse_repeat_rate_corr(mouse,diagonal) = nanmean(diag(mean_rate_corr,diagonal));
        end
    end
    elapse_repeat_rate_corr_area{area} = elapse_repeat_rate_corr;
end

new_colors = [0 0.7 0.8;0.6 0.6 0.6;0.6 0.2 0.6;0.4 0.4 0.4];
new_colors2 = [0 0.5 0.6;0.5 0.5 0.5;0.5 0.1 0.4;0.3 0.3 0.3];
plt = [];
rate_similarity_index = {};
figure('units','normalized','position',[0.3 0.3 0.25 0.4])
for area = 1:4
    current_area = elapse_repeat_rate_corr_area{area};
    rate_similarity_index{area} = [current_area-current_area(:,1)]./[current_area+current_area(:,1)];
    
    mean_elapse_repeat = nanmean(rate_similarity_index{area});
    std_elapse_repeat = nanstd(rate_similarity_index{area});
    ste_elapse_repeat = std_elapse_repeat./sqrt(size(current_area,1));
    
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


V1_rate_corr = rate_similarity_index{1};
V2_rate_corr = rate_similarity_index{2};

LGN_rate_corr = rate_similarity_index{3};
LP_rate_corr = rate_similarity_index{4};

for trial = 1:29
    p = ranksum(V1_rate_corr(:,trial),V2_rate_corr(:,trial));
    if p < 0.05
        hold on
        plot(trial,0,'*','color',new_colors(1,:),'markersize',8)
    end
    
    p = ranksum(LGN_rate_corr(:,trial),LP_rate_corr(:,trial));
    if p < 0.05
        hold on
        plot(trial,0.005,'*','color',new_colors(3,:),'markersize',8)
    end
end


%% Figure 5B - Visual hierarchy - tuning curve correlation between movie repeats - Neuropixels
cell_cutoff = 15;
nat_movie = 1;
num_repeats = 30;
area_list = [1,2,7,8];

for area = 1:4
    valid_mice = (neuropixels_cell_count(:,area_list(area),nat_movie) >= cell_cutoff) & movie_repeats(:,nat_movie) == num_repeats;
    current_area = neuropixels_population_vectors(valid_mice,area_list(area),nat_movie);
    elapse_repeat_tuning_corr = [];
    
    for mouse = 1:size(current_area,1)
        clc;
        disp(['Calculating tuning curve correlation between repeats for Neuropixels data:'])
        disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area_list(area)},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse};
        current_mouse_blockA = current_mouse(:,:,1:30);
        current_mouse_blockB = current_mouse(:,:,31:60);
        mean_tuning_corr = [];
        for repeat1 = 1:num_repeats
            for repeat2 = 1:num_repeats
                tuning_corr_blockA = corr(current_mouse_blockA(:,:,repeat1)',current_mouse_blockA(:,:,repeat2)');
                mean_tuning_corr(repeat1,repeat2,1) = nanmedian(diag(tuning_corr_blockA));
                tuning_corr_blockB = corr(current_mouse_blockB(:,:,repeat1)',current_mouse_blockB(:,:,repeat2)');
                mean_tuning_corr(repeat1,repeat2,2) = nanmedian(diag(tuning_corr_blockB));
            end
        end
        mean_tuning_corr(mean_tuning_corr<0) = 0;
        mean_tuning_corr = nanmean(mean_tuning_corr,3);
        
        for diagonal = 1:29
            elapse_repeat_tuning_corr(mouse,diagonal) = nanmean(diag(mean_tuning_corr,diagonal));
        end
    end
    elapse_repeat_tuning_corr_area{area} = elapse_repeat_tuning_corr;
end

new_colors = [0 0.7 0.8;0.6 0.6 0.6;0.6 0.2 0.6;0.4 0.4 0.4];
new_colors2 = [0 0.5 0.6;0.5 0.5 0.5;0.5 0.1 0.4;0.3 0.3 0.3];
plt = [];
tuning_similarity_index = {};
figure('units','normalized','position',[0.3 0.3 0.25 0.4])
for area = 1:4
    current_area = elapse_repeat_tuning_corr_area{area};
    tuning_similarity_index{area} = [current_area-current_area(:,1)]./[current_area+current_area(:,1)];
    
    mean_elapse_repeat = nanmean(tuning_similarity_index{area});
    std_elapse_repeat = nanstd(tuning_similarity_index{area});
    ste_elapse_repeat = std_elapse_repeat./sqrt(size(current_area,1));
    
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


V1_tuning_corr = tuning_similarity_index{1};
V2_tuning_corr = tuning_similarity_index{2};

LGN_tuning_corr = tuning_similarity_index{3};
LP_tuning_corr = tuning_similarity_index{4};

for trial = 1:29
    p = ranksum(V1_tuning_corr(:,trial),V2_tuning_corr(:,trial));
    if p < 0.05
        hold on
        plot(trial,0,'*','color',new_colors(1,:),'markersize',8)
    end
    
    p = ranksum(LGN_tuning_corr(:,trial),LP_tuning_corr(:,trial));
    if p < 0.05
        hold on
        plot(trial,-0.4,'*','color',new_colors(3,:),'markersize',8)
    end
end
lgd = legend(plt,brain_areas(area_list));
legend('boxoff')
lgd.Position = [0.2 0.225 0.1 0.1];

%% Figure 5C - Visual hierarchy - ensemble rate correlation between movie repeats - calcium imaging

nat_movie = 1;
cell_cutoff = 20;

for area = 1:2
    valid_mice = min(calcium_excitatory_cell_count{area},[],2) >= cell_cutoff;
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie);
    
    elapse_repeat_rate_corr = [];
    for mouse = 1:length(current_area)
        clc;
        disp(['Calculating ensemble rate correlation between repeats for calcium imaging data:'])
        disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area_list(area)},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse};
        
        current_mouse_sess1 =  current_mouse(:,:,1:10);
        valid_cells = nanmean(squeeze(nanmean(current_mouse_sess1,2)),2)>0;
        valid_current_mouse_sess1 = squeeze(nanmean(current_mouse_sess1(valid_cells,:,:),2));
        
        current_mouse_sess2 =  current_mouse(:,:,11:20);
        valid_cells = nanmean(squeeze(nanmean(current_mouse_sess2,2)),2)>0;
        valid_current_mouse_sess2 = squeeze(nanmean(current_mouse_sess2(valid_cells,:,:),2));
        
        current_mouse_sess3 =  current_mouse(:,:,21:30);
        valid_cells = nanmean(squeeze(nanmean(current_mouse_sess3,2)),2)>0;
        valid_current_mouse_sess3 = squeeze(nanmean(current_mouse_sess3(valid_cells,:,:),2));
        
        
        rate_corr_across_session = [];
        for repeat1 = 1:10
            for repeat2 = 1:10
                
                rate_corr_across_session (repeat1,repeat2,1) = corr(valid_current_mouse_sess1(:,repeat1),valid_current_mouse_sess1(:,repeat2));
                rate_corr_across_session(repeat1,repeat2,2) =  corr(valid_current_mouse_sess2(:,repeat1),valid_current_mouse_sess2(:,repeat2));
                rate_corr_across_session(repeat1,repeat2,3) = corr(valid_current_mouse_sess3(:,repeat1),valid_current_mouse_sess3(:,repeat2));
                
            end
        end
        rate_corr_across_session(rate_corr_across_session<0) = 0;
        mean_rate_corr_across_mice = nanmean(rate_corr_across_session,3);
        
        for diagonal = 1:9
            elapse_repeat_rate_corr(mouse,diagonal) = nanmean(diag(mean_rate_corr_across_mice,diagonal));
        end
        
    end
    elapse_repeat_rate_corr_area{area} = elapse_repeat_rate_corr;
end

new_colors = [0 0.7 0.8;0.7 0.7 0.7];
new_colors2 = [0 0.5 0.6;0.5 0.5 0.5];
plt = [];
rate_similarity_index = {};
figure('units','normalized','position',[0.3 0.3 0.25 0.4])
for area = 1:2
    current_area = elapse_repeat_rate_corr_area{area};
    rate_similarity_index{area} = [current_area-current_area(:,1)]./[current_area+current_area(:,1)];
    
    mean_elapse_repeat = nanmean(rate_similarity_index{area});
    std_elapse_repeat = nanstd(rate_similarity_index{area});
    ste_elapse_repeat = std_elapse_repeat./sqrt(size(current_area,1));
    
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


V1_rate_corr = rate_similarity_index{1};
V2_rate_corr = rate_similarity_index{2};

for trial = 1:9
    p = ranksum(V1_rate_corr(:,trial),V2_rate_corr(:,trial));
    if p < 0.05
        hold on
        plot(trial,-0.01,'*','color',[0.4 0.4 0.4],'markersize',8)
    end
    
end
lgd = legend(plt,brain_areas(1:2));
legend('boxoff')
lgd.Position = [0.175 0.19 0.1 0.1];

%% Figure 5D - Visual hierarchy - tuning curve correlation between movie repeats - calcium imaging

nat_movie = 1;
cell_cutoff = 20;

for area = 1:2
    valid_mice = min(calcium_excitatory_cell_count{area},[],2) >= cell_cutoff;
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie);
    
    elapse_repeat_tuning_corr = [];
    for mouse = 1:length(current_area)
        clc;
        disp(['Calculating tuning curve correlation between repeats for calcium imaging data:'])
        disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area_list(area)},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse};
        
        current_mouse_sess1 =  current_mouse(:,:,1:10);
        valid_cells = nanmean(squeeze(nanmean(current_mouse_sess1,2)),2)>0;
        valid_current_mouse_sess1 = current_mouse_sess1(valid_cells,:,:);
        
        current_mouse_sess2 =  current_mouse(:,:,11:20);
        valid_cells = nanmean(squeeze(nanmean(current_mouse_sess2,2)),2)>0;
        valid_current_mouse_sess2 = current_mouse_sess2(valid_cells,:,:);
        
        current_mouse_sess3 =  current_mouse(:,:,21:30);
        valid_cells = nanmean(squeeze(nanmean(current_mouse_sess3,2)),2)>0;
        valid_current_mouse_sess3 = current_mouse_sess3(valid_cells,:,:);
        
        tuning_corr_across_session = [];
        for repeat1 = 1:10
            for repeat2 = 1:10
                tuning_corr = corr(valid_current_mouse_sess1(:,:,repeat1)',valid_current_mouse_sess1(:,:,repeat2)');
                tuning_corr_across_session (repeat1,repeat2,1) = nanmean(diag(tuning_corr));
                
                tuning_corr = corr(valid_current_mouse_sess2(:,:,repeat1)',valid_current_mouse_sess2(:,:,repeat2)');
                tuning_corr_across_session(repeat1,repeat2,2) = nanmean(diag(tuning_corr));
                
                tuning_corr = corr(valid_current_mouse_sess3(:,:,repeat1)',valid_current_mouse_sess3(:,:,repeat2)');
                tuning_corr_across_session(repeat1,repeat2,3) = nanmean(diag(tuning_corr));
                
            end
        end
        tuning_corr_across_session(tuning_corr_across_session<0) = 0;
        mean_tuning_corr_across_mice = nanmean(tuning_corr_across_session,3);
        
        for diagonal = 1:9
            elapse_repeat_tuning_corr(mouse,diagonal) = nanmean(diag(mean_tuning_corr_across_mice,diagonal));
        end
        
        
    end
    elapse_repeat_tuning_corr_area{area} = elapse_repeat_tuning_corr;
end

new_colors = [0 0.7 0.8;0.7 0.7 0.7];
new_colors2 = [0 0.5 0.6;0.5 0.5 0.5];
plt = [];
tuning_similarity_index = {};
figure('units','normalized','position',[0.3 0.3 0.25 0.4])
for area = 1:2
    current_area = elapse_repeat_tuning_corr_area{area};
    tuning_similarity_index{area} = [current_area-current_area(:,1)]./[current_area+current_area(:,1)];
    
    mean_elapse_repeat = nanmean(tuning_similarity_index{area});
    std_elapse_repeat = nanstd(tuning_similarity_index{area});
    ste_elapse_repeat = std_elapse_repeat./sqrt(size(current_area,1));
    
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


V1_tuning_corr = tuning_similarity_index{1};
V2_tuning_corr = tuning_similarity_index{2};

for trial = 1:9
    p = ranksum(V1_tuning_corr(:,trial),V2_tuning_corr(:,trial));
    if p < 0.05
        hold on
        plot(trial,-0.005,'*','color',[0.4 0.4 0.4],'markersize',8)
    end
    
end
lgd = legend(plt,brain_areas(1:2));
legend('boxoff')
lgd.Position = [0.175 0.19 0.1 0.1];


%% Figure 5E - Visual hierarchy - Ensemble rate correlation between blocks - calcium imaging
nat_movie = 2;
cell_cutoff = 20;
within_between_stability_area = {};
for area = 1:2
    valid_mice = calcium_excitatory_cell_count{area}(:,1) >= cell_cutoff;
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie);
    
    within_between_stability = [];
    for mouse = 1:length(current_area)
        clc;
        disp(['Calculating ensemble rate correlation between blocks for calcium imaging data:'])
        disp(['Stimulus: Natural movie 3 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse};
        
        mean_activity = [];
        mean_activity(:,1) = nanmean(nanmean(current_mouse(:,:,1:2),3),2);
        mean_activity(:,2) = nanmean(nanmean(current_mouse(:,:,3:5),3),2);
        mean_activity(:,3) = nanmean(nanmean(current_mouse(:,:,6:7),3),2);
        mean_activity(:,4) = nanmean(nanmean(current_mouse(:,:,8:10),3),2);
        
        rate_corr = corr(mean_activity);
        rate_corr(rate_corr<0) = 0;
        
        within_between_stability(mouse,1) = nanmean([rate_corr(1,2),rate_corr(3,4)]);
        within_between_stability(mouse,2) = nanmean(nanmean(rate_corr(1:2,3:4)));
    end
    within_between_stability_area{area} = within_between_stability;
end

rate_similarity_index = [];
plt = [];
figure('units','normalized','position',[0.3 0.3 0.2 0.3])
for area = [2,1]
    current_area = within_between_stability_area{area};
    
    rate_similarity_index{area} = [current_area(:,2)-current_area(:,1)]./[current_area(:,2) + current_area(:,1)];
    
    mean_stability = nanmean(current_area);
    std_stability = nanstd(current_area);
    ste_stability = std_stability./sqrt(size(current_area,1));
    
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


figure('units','normalized','position',[0.5 0.4 0.1 0.2])
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

[pval,~,stat] = ranksum(rate_similarity_index{1},rate_similarity_index{2});
zval = stat.zval;

clc;
disp(['Ensemble rate similarity index of V1 compared to LM'])
disp(['Two-tailed Mann-Whitney rank-sum test: Z=',num2str(zval),' p=',num2str(pval)])


%% Figure 5F - Visual hierarchy - Tuning curve correlation between blocks - calcium imaging
nat_movie = 2;
cell_cutoff = 20;
within_between_stability_area = {};
for area = 1:2
    valid_mice = calcium_excitatory_cell_count{area}(:,1) >= cell_cutoff;
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie);
    
    within_between_stability = [];
    for mouse = 1:length(current_area)
        clc;
        disp(['Calculating tuning curve correlation between blocks for calcium imaging data:'])
        disp(['Stimulus: Natural movie 3 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse};
        
        tuning_corr_across_mice = [];
        tuning_corr_across_mice(:,:,1) = nanmean(current_mouse(:,:,1:2),3);
        tuning_corr_across_mice(:,:,2) = nanmean(current_mouse(:,:,3:5),3);
        tuning_corr_across_mice(:,:,3) = nanmean(current_mouse(:,:,6:7),3);
        tuning_corr_across_mice(:,:,4) = nanmean(current_mouse(:,:,8:10),3);
        
        mean_tuning_corr_across_mice = [];
        for half1 = 1:4
            for half2 = 1:4
                tuning_corr =  corr(tuning_corr_across_mice(:,:,half1)',tuning_corr_across_mice(:,:,half2)');
                mean_tuning_corr_across_mice(half1,half2) = nanmedian(diag(tuning_corr));
            end
        end
        mean_tuning_corr_across_mice(mean_tuning_corr_across_mice<0) = 0;
        
        within_between_stability(mouse,1) = nanmean([ mean_tuning_corr_across_mice(1,2), mean_tuning_corr_across_mice(3,4)]);
        within_between_stability(mouse,2) = nanmean(nanmean( mean_tuning_corr_across_mice(1:2,3:4)));
    end
    within_between_stability_area{area} = within_between_stability;
end

tuning_similarity_index = [];
plt = [];
figure('units','normalized','position',[0.3 0.3 0.2 0.3])
for area = [2,1]
    current_area = within_between_stability_area{area};
    
    tuning_similarity_index{area} = [current_area(:,2)-current_area(:,1)]./[current_area(:,2) + current_area(:,1)];
    
    mean_stability = nanmean(current_area);
    std_stability = nanstd(current_area);
    ste_stability = std_stability./sqrt(size(current_area,1));
    
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


figure('units','normalized','position',[0.5 0.4 0.1 0.2])
xlim([0 3])
hold on
plot(xlim,[0 0 ],'--','color',[0.2 0.2 0.2])
stability_diff_areas = nan(size(tuning_similarity_index{1},1),2);
stability_diff_areas(:,1) = tuning_similarity_index{1};
stability_diff_areas(1:size(tuning_similarity_index{2},1),2) = tuning_similarity_index{2};
figure_boxplot(stability_diff_areas)
ylabel('Similarity index')
set(gca,'xtick',1:2,'xticklabel',{'V1','LM'})
ylim([-0.15 0.125])

[pval,~,stat] = ranksum(tuning_similarity_index{1},tuning_similarity_index{2});
zval = stat.zval;

clc;
disp(['Tuning curve similarity index of V1 compared to LM'])
disp(['Two-tailed Mann-Whitney rank-sum test: Z=',num2str(zval),' p=',num2str(pval)])


%% Figure 5G - Visual hierarchy - Ensemble rate correlation between sessions - calcium imaging
nat_movie = 1;
within_between_session_stability_area = {};
for area = 1:2
    valid_mice = min(calcium_excitatory_cell_count{area},[],2)>= cell_cutoff;
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie);
    within_between_session_stability = [];
    
    for mouse = 1:length(current_area)
        clc;
        disp(['Calculating ensemble rate correlation between sessions for calcium imaging data:'])
        disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse};
        
        mean_activity_per_half = [];
        for half = 1:6
            current_half = [1:5] + 5*(half-1);
            mean_activity_per_half(:,:,half) = nanmean(current_mouse(:,:,current_half),3);
        end
        
        rate_corr = [];
        for halfA = 1:size(mean_activity_per_half,3)
            halfA_activity = mean_activity_per_half(:,:,halfA);
            for halfB = 1:size(mean_activity_per_half,3)
                halfB_activity = mean_activity_per_half(:,:,halfB);
                
                valid_cells = [nanmean(halfA_activity,2) ~= 0] & [nanmean(halfB_activity,2) ~= 0];
                valid_halfA_activity = nanmean(halfA_activity(valid_cells,:),2);
                valid_halfB_activity = nanmean(halfB_activity(valid_cells,:),2);
                
                rate_corr(halfA,halfB) = corr(valid_halfA_activity,valid_halfB_activity);
                
            end
        end
        rate_corr(rate_corr<0) = 0;
        
        
        within_between_session_stability(mouse,1) = nanmean([rate_corr(1,2),rate_corr(3,4),rate_corr(5,6)]);
        within_between_session_stability(mouse,2) = nanmean([nanmean(nanmean(rate_corr(1:2,3:4))),...
            nanmean(nanmean(rate_corr(3:4,5:6)))]);
        within_between_session_stability(mouse,3) = nanmean(nanmean(rate_corr(1:2,5:6)));
    end
    within_between_session_stability_area{area} = within_between_session_stability;
end

rate_similarity_index = {};
figure('units','normalized','position',[0.3 0.3 0.2 0.3])
for area = [2,1]
    current_area = within_between_session_stability_area{area};
    
    rate_similarity_index{area} = [current_area-current_area(:,1)]./[current_area+current_area(:,1)];
    
    mean_stability = nanmean(rate_similarity_index{area});
    std_stability = nanstd(rate_similarity_index{area});
    ste_stability = std_stability ./sqrt(size(rate_similarity_index{area},1));
    
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

p=[];
for sess = 1:3
    p(sess) = ranksum(rate_similarity_index{1}(:,sess),rate_similarity_index{2}(:,sess));
    if p(sess) < 0.05
        hold on
        plot(sess,0,'*','color',[0.4 0.4 0.4],'markersize',8)
    end
end

%% Figure 5H - Visual hierarchy - Tuning curve correlation between sessions - calcium imaging
nat_movie = 1;
within_between_session_stability_area = {};
for area = 1:2
    valid_mice = min(calcium_excitatory_cell_count{area},[],2)>= cell_cutoff;
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie);
    within_between_session_stability = [];
    
    for mouse = 1:length(current_area)
        clc;
        disp(['Calculating tuning curve correlation between sessions for calcium imaging data:'])
        disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse};
        
        mean_activity_per_half = [];
        for half = 1:6
            current_half = [1:5] + 5*(half-1);
            mean_activity_per_half(:,:,half) = nanmean(current_mouse(:,:,current_half),3);
        end
        
        tuning_corr = [];
        for halfA = 1:size(mean_activity_per_half,3)
            halfA_activity = mean_activity_per_half(:,:,halfA);
            for halfB = 1:size(mean_activity_per_half,3)
                halfB_activity = mean_activity_per_half(:,:,halfB);
                
                valid_cells = [nanmean(halfA_activity,2) ~= 0] & [nanmean(halfB_activity,2) ~= 0];
                valid_halfA_activity = halfA_activity(valid_cells,:);
                valid_halfB_activity = halfB_activity(valid_cells,:);
                
                tuning_corr(halfA,halfB) = nanmedian(diag(corr(valid_halfA_activity',valid_halfB_activity')));
                
            end
        end
        tuning_corr(tuning_corr<0) = 0;
        
        
        within_between_session_stability(mouse,1) = nanmean([tuning_corr(1,2),tuning_corr(3,4),tuning_corr(5,6)]);
        within_between_session_stability(mouse,2) = nanmean([nanmean(nanmean(tuning_corr(1:2,3:4))),...
            nanmean(nanmean(tuning_corr(3:4,5:6)))]);
        within_between_session_stability(mouse,3) = nanmean(nanmean(tuning_corr(1:2,5:6)));
    end
    within_between_session_stability_area{area} = within_between_session_stability;
end

rate_similarity_index = {};
figure('units','normalized','position',[0.3 0.3 0.2 0.3])
for area = [2,1]
    current_area = within_between_session_stability_area{area};
    
    tuning_similarity_index{area} = [current_area-current_area(:,1)]./[current_area+current_area(:,1)];
    
    mean_stability = nanmean(tuning_similarity_index{area});
    std_stability = nanstd(tuning_similarity_index{area});
    ste_stability = std_stability ./sqrt(size(tuning_similarity_index{area},1));
    
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

p=[];
for sess = 1:3
    p(sess) = ranksum(tuning_similarity_index{1}(:,sess),tuning_similarity_index{2}(:,sess));
    if p(sess) < 0.05
        hold on
        plot(sess,-0.1,'*','color',[0.4 0.4 0.4],'markersize',8)
    end
end

%% Figure 6A - Internal structure in reduced space - single animal v1 example
current_area = calcium_excitatory_population_vectors{1}{91,3}*30;
valid_cells = nanmean(nanmean(current_area(:,:,1:10),3),2)>0;
valid_cells = current_area(valid_cells,:,1:10);

current_trial_sorted = [];
for cell = 1:size(valid_cells,1)
    current_trial_sorted(cell,:) = imgaussfilt(valid_cells(cell,:,4),6);
end
[~,i] = max(current_trial_sorted,[],2);
[~,sorted_ind1] = sort(i);


current_trial_sorted = [];
for cell = 1:size(valid_cells,1)
    current_trial_sorted(cell,:) = imgaussfilt(valid_cells(cell,:,6),6);
end
[~,i] = max(current_trial_sorted,[],2);
[~,sorted_ind2] = sort(i);


current_trial_sorted = [];
for cell = 1:size(valid_cells,1)
    current_trial_sorted(cell,:) = imgaussfilt(valid_cells(cell,:,8),6);
end
[~,i] = max(current_trial_sorted,[],2);
[~,sorted_ind3] = sort(i);

temp =[4,6,8];

figure('units','normalized','position',[0.3 0.3 0.275 0.4])
for trial = 1:3
    subplot(2,3,trial)
    current_trial = valid_cells(:,:,temp(trial));
    
    current_trial_sorted = current_trial(sorted_ind1,:);
    imagesc(current_trial_sorted,[0 30])
    title(['Repeat #',num2str(temp(trial))])
    if trial == 1
        ylabel('Sorted by repeat 4')
    else
        set(gca,'ytick',[])
    end
    
    subplot(2,3,trial+3)
    current_trial = valid_cells(:,:,temp(trial));
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


current_area = calcium_excitatory_population_vectors{1}{91,3};
load('Figure5A_12.mat','state')
rng(state)
temp = [];
temp_colors = [];
for repeat = 1:10
    temp = [temp,current_area(:,1:90,repeat)];
    temp_colors = [temp_colors,[1:90]];
end

dim_reduce = tsne(temp','Algorithm','exact','Distance','cosine',...
    'NumDimensions',2,'NumPCAComponents',20,'Perplexity',100);
figure
scatter(-dim_reduce(:,1),dim_reduce(:,2),15,temp_colors,'filled')
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


%% Figure 5C - internal structures of a single neuropixels mouse
nat_movie = 1;
num_repeats = 30;
mouse = 35;
sub = 1;
structures_area_mats = [];
structures_area_labels = [];
triu_ind = boolean(triu(ones(30),1));
for area = 1:6
    current_area = neuropixels_population_vectors{mouse,area,nat_movie};
    for repeat = 1:size(current_area,3)
        current_repeat = current_area(:,:,repeat);
        current_structure = corr(current_repeat);
        structures_area_mats(:,sub) = current_structure(triu_ind);
        structures_area_labels(sub) = area;
        sub = sub + 1;
    end
end



%
% state = rng;
% save('Figure5C_1.mat','state')
load('Figure6C_1.mat','state')
rng(state)
dim_reduce = tsne(structures_area_mats','Algorithm','exact','Distance','Cosine',...
    'NumDimensions',3,'NumPCAComponents',10,'Perplexity',30);
figure('units','normalized','position',[0.4 0.55 0.225 0.35])
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

figure('units','normalized','position',[0.4 0.25 0.1 0.175])
imagesc(corr(neuropixels_population_vectors{mouse,1,nat_movie}(:,:,49)))
xlabel('Time in movie (sec)')
ylabel('Time in movie (sec)')
title('V1 - repeat 49')
colormap(newmap3)

figure('units','normalized','position',[0.525 0.25 0.1 0.175])
imagesc(corr(neuropixels_population_vectors{mouse,3,nat_movie}(:,:,37)))
xlabel('Time in movie (sec)')
ylabel('Time in movie (sec)')
title('AL - repeat 37')
colormap(newmap3)

%% Figure 6D - Internal structures of example neuropixels pseudo-mouse

load('Figure_6D_19.mat','state')
rng(state)

nat_movie = 1;
subset_population_vectors = neuropixels_population_vectors(movie_repeats(:,nat_movie) == 30,1:6,nat_movie);

pseudo_mouseA_ind = sort(randperm(size(subset_population_vectors,1),size(subset_population_vectors,1)./2));
pseudo_mouseB_ind = find(~ismember([1:size(subset_population_vectors,1)],pseudo_mouseA_ind));

pseudo_area_cell_num = [];
pseudo_mouseA = {};
pseudo_mouseB = {};
for area = 1:6
    pseudo_mouseA{area} = cell2mat(subset_population_vectors(pseudo_mouseA_ind,area,nat_movie));
    pseudo_mouseB{area} = cell2mat(subset_population_vectors(pseudo_mouseB_ind,area,nat_movie));
    pseudo_area_cell_num(1,area) = size(pseudo_mouseA{area},1);
    pseudo_area_cell_num(2,area) = size(pseudo_mouseB{area},1);
end

min_cell_num = min(pseudo_area_cell_num(:));
pseudo_mouseA_subset = {};
pseudo_mouseB_subset = {};
for area = 1:6
    subset_cell_ids_mouseA = sort(randperm(pseudo_area_cell_num(1,area),min_cell_num));
    pseudo_mouseA_subset{area} = pseudo_mouseA{area}(subset_cell_ids_mouseA,:,:);
    
    subset_cell_ids_mouseB = sort(randperm(pseudo_area_cell_num(2,area),min_cell_num));
    pseudo_mouseB_subset{area} = pseudo_mouseB{area}(subset_cell_ids_mouseB,:,:);
end

internal_structures_pseudoA = [];
internal_structures_pseudoB = [];
internal_structures_labels = [];
sub = 1;
triu_ind = boolean(triu(ones(30),1));
for area = 1:6
    temp = [];
    for repeat = 1:60
        current_structure_mouseA = corr(pseudo_mouseA_subset{area}(:,:,repeat));
        internal_structures_pseudoA(:,sub) = current_structure_mouseA(triu_ind);
        
        current_structure_mouseB = corr(pseudo_mouseB_subset{area}(:,:,repeat));
        internal_structures_pseudoB(:,sub) = current_structure_mouseB(triu_ind);
        
        internal_structures_labels(sub) = area;
        sub = sub + 1;
    end
    
end

all_real_internal_structures = [internal_structures_pseudoA,internal_structures_pseudoB];
all_internal_structures_labels = [internal_structures_labels,internal_structures_labels+6];

dim_reduce = tsne(all_real_internal_structures','Algorithm','exact','Distance','Cosine',...
    'NumDimensions',3,'NumPCAComponents',20,'Perplexity',30);

plt = [];
figure('units','normalized','position',[0.3 0.3 0.215 0.355])
for area = 1:6
    current_structuresA = [all_internal_structures_labels ==area];
    current_structuresB = [all_internal_structures_labels ==area+6];
    hold on
    plt(area) = scatter3(dim_reduce(current_structuresA,3),dim_reduce(current_structuresA,2),dim_reduce(current_structuresA,1),60,colors(area,:),'filled');
    scatter3(dim_reduce(current_structuresB,3),dim_reduce(current_structuresB,2),dim_reduce(current_structuresB,1),60,colors2(area,:),'filled')
end

grid on
ax = gca;
ax.GridAlpha = 0.1;
%set(gca, 'ydir', 'reverse')
xlabel('Component 1')
ylabel('Component 2')
zlabel('Component 3')
view([2 24.4])
legend(plt,brain_areas(1:6),'Location','best')
legend('boxoff')

%% Figure 6E - Internal structure similarity across technologies
pseudo_mice = {};
for area = 1:6
    pseudo_mice{1,area} = cell2mat(calcium_excitatory_population_vectors{area}(:,1));

    valid_mice = neuropixels_cell_count(:,area,1) >0;
    current_area = neuropixels_population_vectors(valid_mice,area);
    pseudo_mice{2,area} = cell2mat(cellfun(@(x) x(:,:,1:20),current_area,'UniformOutput',false));
  end

calcium_internal_structure = [];
calcium_area_labels = [];

neuropixels_internal_structure = [];
neuropixels_area_labels = [];
triu_ind = boolean(triu(ones(30),1));
for area = 1:6
    calcium_current_area = pseudo_mice{1,area};
    internal_structure_area = [];
    for repeat = 1:30
        current_structure = corr(calcium_current_area(:,:,repeat));
        internal_structure_area(:,repeat) = current_structure(triu_ind);
    end
    calcium_internal_structure = [calcium_internal_structure,internal_structure_area];
    calcium_area_labels = [calcium_area_labels,ones(1,30)*area];
    
    
    neuropixels_current_area = pseudo_mice{2,area};
    internal_structure_area = [];
    for repeat = 1:20
        current_structure = corr(neuropixels_current_area(:,:,repeat));
        internal_structure_area(:,repeat) = current_structure(triu_ind);
    end
    neuropixels_internal_structure = [neuropixels_internal_structure,internal_structure_area];
    neuropixels_area_labels = [neuropixels_area_labels,ones(1,20)*area];
    
end

neuropixels_internal_structure_zscore = zscore(neuropixels_internal_structure,[],2);
calcium_internal_structure_zscore = zscore(calcium_internal_structure,[],2);

calcium_median_internal_structure= [];
neuropixels_median_internal_structure = [];
for area = 1:6
    current_area_calcium = calcium_internal_structure_zscore(:,calcium_area_labels == area);
    calcium_median_internal_structure(:,area)=  nanmedian(current_area_calcium,2);
    
    current_area_neuropixels = neuropixels_internal_structure_zscore(:,neuropixels_area_labels == area);
    neuropixels_median_internal_structure(:,area) = nanmedian(current_area_neuropixels,2);
end


normalized_internal_structures_both_tech = [neuropixels_internal_structure_zscore,calcium_internal_structure_zscore];
area_labels_both_tech =[neuropixels_area_labels,calcium_area_labels+6];


load('figure6E_7.mat','state')
rng(state)
dim_reduce = tsne(normalized_internal_structures_both_tech','Algorithm','exact','Distance','correlation',...
    'NumDimensions',3,'NumPCAComponents',20,'Perplexity',40);
plt = [];
figure('units','normalized','position',[0.3 0.3 0.25 0.4])
for area =  1:6
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
%view([15 20 2.5])
view([15 10 3])

similarity_across_technologies = pdist2(calcium_median_internal_structure',neuropixels_median_internal_structure','correlation');
figure('units','normalized','position',[0.55 0.4 0.2 0.275])
imagesc(similarity_across_technologies)
colormap(flipud(newmap3))
set(gca,'xtick',1:6,'xticklabel',brain_areas(1:6),'ytick',1:6,'yticklabel',brain_areas(1:6))
xlabel('Calcium imaging')
ylabel('Neuropixels')
cb = colorbar;
cb.Label.String = 'Correlation distance';
cb.FontSize = 10;


%% Figure 6F - Internal structures of two pseudo-mice (reduced and matrices)
nat_movie = 1;
subset_population_vectors = {};
for area = 1:6
    for mouse = 1:size(neuropixels_population_vectors_tsne,1)
        if ~isempty(neuropixels_population_vectors_tsne{mouse,area,nat_movie}) && movie_repeats(mouse,nat_movie) == 30
            subset_population_vectors{mouse,area} = neuropixels_population_vectors_tsne{mouse,area,nat_movie}(:,:,31:60);
        end

    end
end
load('neuropixels_structure_ind.mat','pseudo_mouseA_ind','pseudo_mouseB_ind')
pseudo_mouseA = {};
pseudo_mouseB = {};
for area = 1:6
    pseudo_mouseA{area} = cell2mat(subset_population_vectors(pseudo_mouseA_ind,area));
    pseudo_mouseB{area} = cell2mat(subset_population_vectors(pseudo_mouseB_ind,area));
end

dim_reduce_areas = {};
for area = 1:6
 load(['neuropixels_structure_STATE',num2str(area),'.mat'],'state')
 rng(state)
    
    pseudo_mouseA_strcut = reshape(pseudo_mouseA{area},[size(pseudo_mouseA{area},1),90*30]);
    dim_reduceA = tsne(pseudo_mouseA_strcut','Algorithm','exact','Distance','cosine',...
        'NumDimensions',3,'NumPCAComponents',20,'Perplexity',200);
    dim_reduce_areas{1,area} = dim_reduceA;
   
    
    pseudo_mouseB_strcut = reshape(pseudo_mouseB{area},[size(pseudo_mouseB{area},1),90*30]);
    dim_reduceB = tsne(pseudo_mouseB_strcut','Algorithm','exact','Distance','cosine',...
        'NumDimensions',3,'NumPCAComponents',20,'Perplexity',200);
    dim_reduce_areas{2,area} = dim_reduceB;

end

view_list_pseudoA = [-127.1 24.4;-174.7 63.6;4.9 19.6;-96.3 57.2;-22.7 6;34.5 52.4];
view_list_pseudoB = [55.3 -32.4;79.7 70.8;-3.5 -25.2;133.2 26;31.6 11.2;61.2 13.2];
figure('units','normalized','position',[0.1 0.5 0.8 0.4])
for area = 1:6
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

figure('units','normalized','position',[0.1 0 0.8 0.4])
for area = 1:6
     pseudo_mouseA_strcut = corr(nanmean(cell2mat(pseudo_mouseA(area)),3));
    pseudo_mouseB_strcut = corr(nanmean(cell2mat(pseudo_mouseB(area)),3));
   
    subplot(2,6,area)
    imagesc(pseudo_mouseA_strcut)
    subplot(2,6,area+6)
    imagesc(pseudo_mouseB_strcut)
end
colormap(newmap3)

%% Figure 6G - between pseudo-mice permutation decoder

nat_movie = 1;
subset_population_vectors = {};
for area = 1:6
    for mouse = 1:size(neuropixels_population_vectors,1)
        if ~isempty(neuropixels_population_vectors{mouse,area,nat_movie})
            subset_population_vectors{mouse,area} = neuropixels_population_vectors{mouse,area,nat_movie}(:,:,1:20);
        end
    end
end

similarity_between_pseudomice = [];
similarity_between_shuffled_pseudomice = [];

num_shuffles = 1000;
for shuffle = 1:num_shuffles
    clc;shuffle
    pseudo_mouseA_ind = sort(randperm(size(subset_population_vectors,1),size(subset_population_vectors,1)./2));
    pseudo_mouseB_ind = find(~ismember([1:size(subset_population_vectors,1)],pseudo_mouseA_ind));
    
    pseudo_area_cell_num = [];
    pseudo_mouseA = {};
    pseudo_mouseB = {};
    for area = 1:6
        pseudo_mouseA{area} = cell2mat(subset_population_vectors(pseudo_mouseA_ind,area,nat_movie));
        pseudo_mouseB{area} = cell2mat(subset_population_vectors(pseudo_mouseB_ind,area,nat_movie));
        pseudo_area_cell_num(1,area) = size(pseudo_mouseA{area},1);
        pseudo_area_cell_num(2,area) = size(pseudo_mouseB{area},1);
    end
    
    min_cell_num = min(pseudo_area_cell_num(:));
    pseudo_mouseA_subset = {};
    pseudo_mouseB_subset = {};
    for area = 1:6
        subset_cell_ids_mouseA = sort(randperm(pseudo_area_cell_num(1,area),min_cell_num));
        pseudo_mouseA_subset{area} = pseudo_mouseA{area}(subset_cell_ids_mouseA,:,:);
        
        subset_cell_ids_mouseB = sort(randperm(pseudo_area_cell_num(2,area),min_cell_num));
        pseudo_mouseB_subset{area} = pseudo_mouseB{area}(subset_cell_ids_mouseB,:,:);
    end
    
    all_cells_pseudo_mouseA = cell2mat(pseudo_mouseA');
    all_cells_pseudo_mouseB = cell2mat(pseudo_mouseB');
    rand_cells_id_mouseA = randperm(size(all_cells_pseudo_mouseA,1));
    rand_cells_id_mouseB = randperm(size(all_cells_pseudo_mouseB,1));
    
    shuffled_pseudo_mouseA_subset = {};
    shuffled_pseudo_mouseB_subset = {};
    for area = 1:6
        current_pseudo_area = [1:min_cell_num] + min_cell_num*(area-1);
        
        shuffled_pseudo_mouseA_subset{area} = all_cells_pseudo_mouseA(rand_cells_id_mouseA(current_pseudo_area),:,:);
        shuffled_pseudo_mouseB_subset{area} = all_cells_pseudo_mouseB(rand_cells_id_mouseB(current_pseudo_area),:,:);
    end
    
    internal_structures_pseudoA = [];
    internal_structures_pseudoB = [];
    internal_structures_pseudoA_shuffle = [];
    internal_structures_pseudoB_shuffle = [];
    internal_structures_labels = [];
    triu_ind = boolean(triu(ones(30),1));
    for area = 1:6
        current_structure_mouseA = corr(nanmean(pseudo_mouseA_subset{area},3));
        internal_structures_pseudoA(:,area) = current_structure_mouseA(triu_ind);
        
        current_structure_mouseB = corr(nanmean(pseudo_mouseB_subset{area},3));
        internal_structures_pseudoB(:,area) = current_structure_mouseB(triu_ind);
        
        current_structure_mouseA_shuffle = corr(nanmean(shuffled_pseudo_mouseA_subset{area},3));
        internal_structures_pseudoA_shuffle(:,area) = current_structure_mouseA_shuffle(triu_ind);
        
        
        current_structure_mouseB_shuffle = corr(nanmean(shuffled_pseudo_mouseB_subset{area},3));
        internal_structures_pseudoB_shuffle(:,area) = current_structure_mouseB_shuffle(triu_ind);
    end
    
    
    similarity_sum = [];
    similarity_sum_shuffled = [];
    permutations = flipud(perms([1:6]));
    for perm = 1:size(permutations,1)
        current_perm = permutations(perm,:);
        similarity_between_pseudomice = corr(internal_structures_pseudoA,internal_structures_pseudoB(:,current_perm));
        similarity_between_pseudomice_shuffled = corr(internal_structures_pseudoA_shuffle,internal_structures_pseudoB_shuffle(:,current_perm));
        
        similarity_sum(perm) = sum(diag(similarity_between_pseudomice));
        similarity_sum_shuffled(perm) = sum(diag(similarity_between_pseudomice_shuffled));
    end
    
    [B,I] = max(similarity_sum);
    between_similarity_acc(shuffle,:) = permutations(I,:) == [1:6];
    
    [B,I] = max(similarity_sum_shuffled);
    between_similarity_acc_shuffled(shuffle,:) = permutations(I,:) == [1:6];
    
end

between_pseudo_mice_decoder_acc = (sum(between_similarity_acc)./size(between_similarity_acc,1))*100;
between_shuffled_pseudo_mice_decoder_acc = (sum(between_similarity_acc_shuffled)./size(between_similarity_acc_shuffled,1))*100;

plt = [];
figure('units','normalized','position',[0.3 0.3 0.2 0.275])
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
nat_movie = 1;
subset_population_vectors = {};
for area = 1:6
    for mouse = 1:size(neuropixels_population_vectors,1)
        if ~isempty(neuropixels_population_vectors{mouse,area,nat_movie})
            subset_population_vectors{mouse,area} = neuropixels_population_vectors{mouse,area,nat_movie}(:,:,1:20);
        end
    end
end

num_shuffles = 3000;
cell_cutoff_list = [5,10,25,50,75,100,150,200,250,300,350,400,500,600,700];
between_similarity_acc_all_cutoffs = [];
between_similarity_acc_shuffled_all_cutoffs = [];

min_cell_num_all_shuffles = [];
for cells_included = 1:length(cell_cutoff_list)+1
    between_similarity_acc = nan(num_shuffles,6);
    between_similarity_acc_shuffled = nan(num_shuffles,6);
    
    for shuffle = 1:num_shuffles
        clc;[cells_included,shuffle]
        pseudo_area_cell_num = 0;
        if cells_included <= length(cell_cutoff_list)
            breaker = cell_cutoff_list(cells_included);
        else
            breaker = 1;
        end
        
        while   breaker > min(pseudo_area_cell_num(:))
            pseudo_mouseA_ind = sort(randperm(size(subset_population_vectors,1),size(subset_population_vectors,1)./2));
            pseudo_mouseB_ind = find(~ismember([1:size(subset_population_vectors,1)],pseudo_mouseA_ind));
            
            
            pseudo_mouseA = {};
            pseudo_mouseB = {};
            for area = 1:6
                pseudo_mouseA{area} = cell2mat(subset_population_vectors(pseudo_mouseA_ind,area,nat_movie));
                pseudo_mouseB{area} = cell2mat(subset_population_vectors(pseudo_mouseB_ind,area,nat_movie));
                pseudo_area_cell_num(1,area) = size(pseudo_mouseA{area},1);
                pseudo_area_cell_num(2,area) = size(pseudo_mouseB{area},1);
            end
        end
        
        if cells_included <= length(cell_cutoff_list)
            min_cell_num = cell_cutoff_list(cells_included);
        else
            min_cell_num = min(pseudo_area_cell_num(:));
            min_cell_num_all_shuffles(shuffle) = min_cell_num;
        end
        
        pseudo_mouseA_subset = {};
        pseudo_mouseB_subset = {};
        for area = 1:6
            subset_cell_ids_mouseA = sort(randperm(pseudo_area_cell_num(1,area),min_cell_num));
            pseudo_mouseA_subset{area} = pseudo_mouseA{area}(subset_cell_ids_mouseA,:,:);
            
            subset_cell_ids_mouseB = sort(randperm(pseudo_area_cell_num(2,area),min_cell_num));
            pseudo_mouseB_subset{area} = pseudo_mouseB{area}(subset_cell_ids_mouseB,:,:);
        end
        
        all_cells_pseudo_mouseA = cell2mat(pseudo_mouseA');
        all_cells_pseudo_mouseB = cell2mat(pseudo_mouseB');
        rand_cells_id_mouseA = randperm(size(all_cells_pseudo_mouseA,1));
        rand_cells_id_mouseB = randperm(size(all_cells_pseudo_mouseB,1));
        
        shuffled_pseudo_mouseA_subset = {};
        shuffled_pseudo_mouseB_subset = {};
        for area = 1:6
            current_pseudo_area = [1:min_cell_num] + min_cell_num*(area-1);
            
            shuffled_pseudo_mouseA_subset{area} = all_cells_pseudo_mouseA(rand_cells_id_mouseA(current_pseudo_area),:,:);
            shuffled_pseudo_mouseB_subset{area} = all_cells_pseudo_mouseB(rand_cells_id_mouseB(current_pseudo_area),:,:);
        end
        
        internal_structures_pseudoA = [];
        internal_structures_pseudoB = [];
        internal_structures_pseudoA_shuffle = [];
        internal_structures_pseudoB_shuffle = [];
        
        internal_structures_labels = [];
        triu_ind = boolean(triu(ones(30),1));
        for area = 1:6
            current_structure_mouseA = corr(nanmean(pseudo_mouseA_subset{area},3));
            internal_structures_pseudoA(:,area) = current_structure_mouseA(triu_ind);
            
            current_structure_mouseB = corr(nanmean(pseudo_mouseB_subset{area},3));
            internal_structures_pseudoB(:,area) = current_structure_mouseB(triu_ind);
            
            current_structure_mouseA_shuffle = corr(nanmean(shuffled_pseudo_mouseA_subset{area},3));
            internal_structures_pseudoA_shuffle(:,area) = current_structure_mouseA_shuffle(triu_ind);
            
            current_structure_mouseB_shuffle = corr(nanmean(shuffled_pseudo_mouseB_subset{area},3));
            internal_structures_pseudoB_shuffle(:,area) = current_structure_mouseB_shuffle(triu_ind);
            
        end
        
        permutations = flipud(perms([1:6]));
        
        % real data
        between_mouse_similarity = corr(internal_structures_pseudoA,internal_structures_pseudoB);
        between_mouse_similarity_shuffled = corr(internal_structures_pseudoA_shuffle,internal_structures_pseudoB_shuffle);
        
        similarity_sum = [];
        similarity_sum_shuffled = [];
        
        for perm = 1:size(permutations,1)
            similarity_sum(perm) = sum(diag(between_mouse_similarity(:,permutations(perm,:))));
            similarity_sum_shuffled(perm) = sum(diag(between_mouse_similarity_shuffled(:,permutations(perm,:))));
            
        end
        
        % real data
        [B,I] = max(similarity_sum);
        between_similarity_acc(shuffle,:) = permutations(I,:) == [1:6];
        
        [B,I] = max(similarity_sum_shuffled);
        between_similarity_acc_shuffled(shuffle,:) = permutations(I,:) == [1:6];
        
    end
    between_similarity_acc_all_cutoffs(cells_included,:) = sum(between_similarity_acc)./num_shuffles;
    between_similarity_acc_shuffled_all_cutoffs(cells_included,:) = sum(between_similarity_acc_shuffled)./num_shuffles;
    
    
end

plt = [];
xvalues = [cell_cutoff_list,round(nanmean(min_cell_num_all_shuffles))];
figure('units','normalized','position',[0.35 0.4 0.25 0.375])
for area = 1:6
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
%% Internal structure in reduced space - single animal v1 example across sessions
% state = rng;
% save('Figure7A_2.mat','state')

load('Figure7A_1.mat','state')
rng(state)

current_area = calcium_excitatory_population_vectors{1}{89,3};
temp = [];
temp_colors = [];
for repeat = 1:30
    temp = [temp,current_area(:,:,repeat)];
    temp_colors = [temp_colors,[1:90]];
end

dim_reduce_sess1 = tsne(temp(:,1:900)','Algorithm','exact','Distance','Cosine',...
    'NumDimensions',2,'NumPCAComponents',10,'Perplexity',200);

dim_reduce_sess2 = tsne(temp(:,901:1800)','Algorithm','exact','Distance','Cosine',...
    'NumDimensions',2,'NumPCAComponents',10,'Perplexity',200);

dim_reduce_sess3 = tsne(temp(:,1801:2700)','Algorithm','exact','Distance','cosine',...
    'NumDimensions',2,'NumPCAComponents',10,'Perplexity',200);

figure('units','normalized','position',[0.25 0.4 0.45 0.225])
subplot(1,3,1)
scatter(dim_reduce_sess1(:,2),dim_reduce_sess1(:,1),5,temp_colors(:,1:900),'filled')
set(gca, 'xdir', 'reverse')
xlim([-10 10])
ylim([-10 10])
text(0.1, 0.1,['Natural movie 1'],'Units','normalized','color',[0.6 0.6 0.6],'fontsize',12)
xlabel('Component 1')
ylabel('Component 2')
title({'Sessio 1 - Day 93'})

subplot(1,3,2)
scatter(dim_reduce_sess2 (:,1),dim_reduce_sess2(:,2),5,temp_colors(:,1:900),'filled')
xlim([-10 10])
ylim([-10 10])
xlabel('Component 1')
ylabel('Component 2')
title({'Sessio 2 - Day 94'})

subplot(1,3,3)
scatter(dim_reduce_sess3(:,2),dim_reduce_sess3(:,1),5,temp_colors(:,1:900),'filled')
xlim([-10 10])
ylim([-8 8])
xlabel('Component 1')
ylabel('Component 2')
title({'Sessio 3 - Day 98'})
set(gca, 'xdir', 'reverse', 'ydir', 'reverse')
colormap(new_jet_colormap)

%% Internal structure versus PV for pseudo-AL

area = 3;
current_area = cell2mat(calcium_excitatory_population_vectors{area}(:,1));
valid_cells = find(sum([nanmean(sum(current_area(:,:,1:10),2),3),...
    nanmean(sum(current_area(:,:,11:20),2),3),...
    nanmean(sum(current_area(:,:,21:30),2),3)]>0,2) == 3);
current_area = current_area(valid_cells,:,:);

rep_list = [1:10;21:30];


figure('units','normalized','position',[0.4 0.4 0.225 0.4])
for sess = 1:2
    current_sess = nanmean(current_area(:,:,rep_list(sess,:)),3);
    
    subplot(2,2,sess,'units','normalized','position',[0.2+0.3*(sess-1) 0.65 0.275  0.26])
    
    internal_structure = corr(current_rep);
    internal_structure (boolean(eye(size(internal_structure ,1)))) = NaN;
    internal_structure(isnan(internal_structure)) = max(internal_structure(:),[],'omitnan');
    imagesc(internal_structure)
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
    
    smooth_current_rep = [];
    for cell = 1:size(current_sess,1)
        smooth_current_rep(cell,:) = imgaussfilt(current_sess(cell,:),3);
    end
    
    smooth_current_rep_norm = smooth_current_rep ./ max(smooth_current_rep,[],2);
    if sess == 1
        [row,col] = find(smooth_current_rep_norm == 1);
        [B,I] = sort(col);
    end
    subplot(2,2,sess+2,'units','normalized','position',[0.2+0.3*(sess-1) 0.15 0.275  0.45])
    imagesc(smooth_current_rep_norm(row(I),:))
    
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
    
    elapsed_session_struc_norm = current_area_structure./current_area_structure(:,1);
    elapsed_session_pv_norm = current_area_pv./current_area_pv(:,1);
    
    mean_struct = nanmean(elapsed_session_struc_norm);
    std_struct = nanstd(elapsed_session_struc_norm);
    
    mean_pv_corr = nanmean(elapsed_session_pv_norm);
    std_pv_corr  = nanstd(elapsed_session_pv_norm);
    
    
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

subset_population_vectors = {};
for area = 1:6
    subset_population_vectors{area} = calcium_excitatory_population_vectors{area}(:,1);
end

triu_id = boolean(triu(ones(30),1));

cell_count_list_full = [25, 50, 100, 250,500,1000,2000,4000,6000, 7900];
elapsed_session_all_measurments = {};

for area = 1:6
    current_area = cell2mat(subset_population_vectors{area});
    
    valid_cells = [nanmean(nanmean(current_area(:,:,1:10),3),2) > 0 &...
        nanmean(nanmean(current_area(:,:,11:20),3),2) > 0 &...
        nanmean(nanmean(current_area(:,:,21:30),3),2) > 0];
    current_area = current_area(valid_cells,:,:);
    cell_count_list = cell_count_list_full(cell_count_list_full<=size(current_area,1));
    
    shuffle_num = 1000;
    elapsed_session_pv = nan(length(cell_count_list),3,shuffle_num);
    elapsed_session_struc = nan(length(cell_count_list),3,shuffle_num);
    
    for cell_count = 1:length(cell_count_list)
        for shuffle = 1:shuffle_num
            clc;
            disp(['Generating calcium imaging pseudo-mice and calculating'])
            disp(['internal structure and population vectors stability:'])
            disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Cells included: ',num2str(cell_count_list_full(cell_count)),'/',num2str(cell_count_list(end)),...
                ' | Realization: ',num2str(shuffle),'/',num2str(num_shuffles)])
            
            
            valid_cells = sort(randperm(size(current_area,1),cell_count_list_full(cell_count)));
            subset_valid_current_area = current_area(valid_cells,:,:);
            
            all_structures = [];
            mean_pv = [];
            mean_tuning = [];
            all_signal_corr = [];
            mean_rate = [];
            
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
            
            elapsed_session_pv(cell_count,1,shuffle) = nanmean([mean_pv(1,2),mean_pv(3,4),mean_pv(5,6)]);
            elapsed_session_pv(cell_count,2,shuffle) = nanmean([nanmean(nanmean(mean_pv(1:2,3:4))),nanmean(nanmean(mean_pv(3:4,5:6)))]);
            elapsed_session_pv(cell_count,3,shuffle) = nanmean([nanmean(nanmean(mean_pv(1:2,5:6)))]);
            
            elapsed_session_struc(cell_count,1,shuffle) = nanmean([structure_stability(1,2),structure_stability(3,4),structure_stability(5,6)]);
            elapsed_session_struc(cell_count,2,shuffle) = nanmean([nanmean(nanmean(structure_stability(1:2,3:4))),nanmean(nanmean(structure_stability(3:4,5:6)))]);
            elapsed_session_struc(cell_count,3,shuffle) = nanmean([nanmean(nanmean(structure_stability(1:2,5:6)))]);
        end
    end
    elapsed_session_all_measurments(area,:) = {elapsed_session_pv,elapsed_session_struc};
    
end


figure('units','normalized','position',[0.35 0.4 0.45 0.45])
for area =  1:6
    curent_area_pv = elapsed_session_all_measurments{area,1};
    curent_area_pv = nanmean(curent_area_pv,3);
    elapsed_session_pv_norm = [curent_area_pv./curent_area_pv(:,1,:)];
    
    
    mean_elapsed_session_pv_norm = nanmean(elapsed_session_pv_norm,3);
    std_elapsed_session_pv_norm = nanstd(elapsed_session_pv_norm,[],3);
    
    curent_area_struct = elapsed_session_all_measurments{area,2};
    curent_area_struct = nanmean(curent_area_struct,3);
    
    elapsed_session_struct_norm = [curent_area_struct./curent_area_struct(:,1,:)];
    mean_elapsed_session_struct_norm = nanmean(elapsed_session_struct_norm,3);
    std_elapsed_session_struct_norm = nanstd(elapsed_session_struct_norm,[],3);
    
    
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
nat_movie = 1;
mouse = 23;
area = 1;
unit = 9;
example_mouse = neuropixels_population_vectors{mouse,area,nat_movie};
current_unit = squeeze(example_mouse(unit,:,:))'.*30;

figure('units','normalized','position',[0.575 0.4 0.1 0.2])
hold on
mean_current_unit = [];
mean_current_unit(1,:) = nanmean(current_unit(1:10,:));
mean_current_unit(2,:) = nanmean(current_unit(11:20,:));
plot(mean_current_unit(1,:),'-','markersize',20,'linewidth',2,'color',[0.2 0.2 0.2])
plot(mean_current_unit(2,:)','-','markersize',20,'linewidth',2,'color',colors(area,:))
[r,p] = corr(mean_current_unit(1,:)',mean_current_unit(2,:)');
text(0.325, 0.675,['r: ',num2str(r)],'Units','normalized','color',[0 0 0])
legend({'Block A','Block B'},'Location','best')
legend('boxoff')
xlim([1 30])
ylabel('Mean activity rate (HZ)')
xlabel('Time in movie (sec)')


current_unit = squeeze(example_mouse(unit,:,11:20))*30';
figure('units','normalized','position',[0.3 0.3 0.25 0.425])
for repeat = 1:10
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
cell_cutoff = 15;
nat_movie = 1;
repeats_num = 30;
area = 1;
current_area = neuropixels_population_vectors(:,area,nat_movie);
valid_mice = neuropixels_cell_count(:,area,nat_movie) >= cell_cutoff & movie_repeats(:,nat_movie) == repeats_num;
all_cells = cell2mat(current_area(valid_mice))*30;

activity_rate = [];
tuning_curve_across_blocks = [];
for cell = 1:size(all_cells,1)
    current_cellA = squeeze(all_cells(cell,:,1:30));
    current_cellB = squeeze(all_cells(cell,:,31:60));
    
    activity_rate(cell,:) = [nanmean(nanmean(current_cellA)),nanmean(nanmean(current_cellB))];
    tuning_curve_across_blocks(cell) = corr(nanmean(current_cellA,2),nanmean(current_cellB,2));
end

text_pos_list = [0.1, 0.95; 0.5, 0.95; 0.5, 0.95];
cell_list = [216,1650,72];
figure('units','normalized','position',[0.25 0.4 0.45 0.25])
for cell = 1:length(cell_list)
    
    current_cellA = nanmean(squeeze(all_cells(cell_list(cell),:,1:30)),2);
    current_cellB = nanmean(squeeze(all_cells(cell_list(cell),:,31:60)),2);
    subplot(1,3,cell)
    
    hold on
    plot(current_cellA,'color',[0.2 0.2 0.2],'linewidth',2)
    plot(current_cellB,'color',colors(1,:),'linewidth',2)
    
    text(text_pos_list(cell,1), text_pos_list(cell,2),['A: ',num2str(nanmean(current_cellA)),' Hz'],'Units','normalized','color',[0.2 0.2 0.2])
    text(text_pos_list(cell,1), text_pos_list(cell,2)-0.05,['B: ',num2str(nanmean(current_cellB)),' Hz'],'Units','normalized','color',colors(1,:))
    text(text_pos_list(cell,1), text_pos_list(cell,2)-0.1,['r: ',num2str(tuning_curve_across_blocks(cell_list(cell)))],'Units','normalized','color',[0 0 0])
    
    if cell ==1
        ylabel('Activity rate (spike/sec)')
        title({'Stable tuning';'stable firing rate'})
    elseif cell == 2
        xlabel('Time in movie (sec)')
        title({'Stable tuning';'unstable firing rate'})
    else
        title({'Unstable tuning';'stable firing rate'})
    end
    xlim([1 30])
end


%% Figure S1F - Relationship between single cell tuning curve stability and single cell activity rate - all V1 units
cell_cutoff = 15;
nat_movie = 1;
repeats_num = 30;
area = 1;
current_area = neuropixels_population_vectors(:,area,nat_movie);
valid_mice = neuropixels_cell_count(:,area,nat_movie) >= cell_cutoff & movie_repeats(:,nat_movie) == repeats_num;
all_cells = cell2mat(current_area(valid_mice))*30;

activity_rate = [];
tuning_curve_across_blocks = [];
for cell = 1:size(all_cells,1)
    current_cellA = squeeze(all_cells(cell,:,1:30));
    current_cellB = squeeze(all_cells(cell,:,31:60));
    activity_rate(cell,:) = [nanmean(nanmean(current_cellA)),nanmean(nanmean(current_cellB))];
    tuning_curve_across_blocks(cell) = corr(nanmean(current_cellA,2),nanmean(current_cellB,2));
end
r  = corr(tuning_curve_across_blocks.^2',nanmean(activity_rate,2));


figure('units','normalized','position',[0.35 0.4 0.25 0.35])
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
nat_movie = 1;
cell_cutoff = 15;
num_repeats = 30;
for area = 1:6
    valid_mice = [neuropixels_cell_count(:,area,nat_movie) >= cell_cutoff] & [movie_repeats(:,nat_movie) == num_repeats];
    current_area = neuropixels_population_vectors(valid_mice,area,nat_movie);
    rate_tuning_relationship = [];
    rate_pv_relationship = [];
    tuning_pv_relationship = [];
    pv_rate_tuning_relationship = [];
    
    for mouse = 1:size(current_area,1)
        
        clc;
        disp(['Calculating the relationship between population vectors,'])
        disp(['ensemble rate and tuning curve correlations'])
        disp(['across movie repeats within a block:'])
        disp(['Dataset: Neuropixels | Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        
        current_mouse = current_area{mouse};
        
        rate_corr_blockA= corr(squeeze(nanmean(current_mouse(:,:,1:30),2)));
        rate_corr_blockB= corr(squeeze(nanmean(current_mouse(:,:,31:60),2)));
        tuning_corr_blockA = [];
        tuning_corr_blockB = [];
        pv_corr_blockA = [];
        pv_corr_blockB = [];
        for repeat1 = 1:30
            for repeat2 = 1:30
                tuning_corr_blockA(repeat1,repeat2) = nanmedian(diag(corr(current_mouse(:,:,repeat1)',current_mouse(:,:,repeat2)')));
                tuning_corr_blockB(repeat1,repeat2) = nanmedian(diag(corr(current_mouse(:,:,repeat1+num_repeats)',current_mouse(:,:,repeat2+num_repeats)')));
                
                pv_corr_blockA(repeat1,repeat2) = nanmean(diag(corr(current_mouse(:,:,repeat1),current_mouse(:,:,repeat2))));
                pv_corr_blockB(repeat1,repeat2) = nanmean(diag(corr(current_mouse(:,:,repeat1+num_repeats),current_mouse(:,:,repeat2+num_repeats))));
            end
        end
        
        for diagonal = 1:20
            mdl = fitlm(diag(rate_corr_blockA,diagonal),diag(tuning_corr_blockA,diagonal));
            rate_tuning_relationship(mouse,diagonal,1) = mdl.Rsquared.Ordinary;
            
            mdl = fitlm(diag(rate_corr_blockB,diagonal),diag(tuning_corr_blockB,diagonal));
            rate_tuning_relationship(mouse,diagonal,2) = mdl.Rsquared.Ordinary;
            
            mdl = fitlm(diag(rate_corr_blockA,diagonal),diag(pv_corr_blockA,diagonal));
            rate_pv_relationship(mouse,diagonal,1) = mdl.Rsquared.Ordinary;
            
            mdl = fitlm(diag(rate_corr_blockB,diagonal),diag(pv_corr_blockB,diagonal));
            rate_pv_relationship(mouse,diagonal,2) = mdl.Rsquared.Ordinary;
            
            mdl = fitlm(diag(tuning_corr_blockA,diagonal),diag(pv_corr_blockA,diagonal));
            tuning_pv_relationship(mouse,diagonal,1) = mdl.Rsquared.Ordinary;
            
            mdl = fitlm(diag(tuning_corr_blockB,diagonal),diag(pv_corr_blockB,diagonal));
            tuning_pv_relationship(mouse,diagonal,2) = mdl.Rsquared.Ordinary;
            
            predictors_blockA = [diag(rate_corr_blockA,diagonal),diag(tuning_corr_blockA,diagonal)];
            mdl = fitlm(predictors_blockA,diag(pv_corr_blockA,diagonal));
            pv_rate_tuning_relationship(mouse,diagonal,1) = mdl.Rsquared.Ordinary;
            
            predictors_blockB = [diag(rate_corr_blockB,diagonal),diag(tuning_corr_blockB,diagonal)];
            mdl = fitlm(predictors_blockB,diag(pv_corr_blockB,diagonal));
            pv_rate_tuning_relationship(mouse,diagonal,2) = mdl.Rsquared.Ordinary;
            
        end
    end
    
    rate_pv_explained = nanmean(nanmean(rate_pv_relationship,2),3);
    tuning_pv_explained = nanmean(nanmean(tuning_pv_relationship,2),3);
    rate_tuning_explained = nanmean(nanmean(rate_tuning_relationship,2),3);
    pv_rate_tuning_explained = nanmean(nanmean(pv_rate_tuning_relationship,2),3);
    
    all_models{area} = [pv_rate_tuning_explained(:),rate_pv_explained(:),tuning_pv_explained(:),rate_tuning_explained(:)];
end

pos_list = [0.1 0.575 0.25 0.325;0.4 0.575 0.25 0.325;0.7 0.575 0.25 0.325;...
    0.1 0.2 0.25 0.325;0.4 0.2 0.25 0.325;0.7 0.2 0.25 0.325];
figure('units','normalized','position',[0.25 0.2 0.45 0.55])
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

cell_cutoff = 20;
for area = 1:6
    valid_mice = min(calcium_excitatory_cell_count{area},[],2) >= cell_cutoff;
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,1);
    pv_rate_tuning_relationships = [];
    for mouse = 1:length(current_area)
        clc;
        disp(['Calculating the relationship between population vectors,'])
        disp(['ensemble rate and tuning curve correlationsacross sessions:'])
        disp(['Dataset: Calcium imaging | Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse};
        
        rate_corr = [];
        pv_corr = [];
        tuning_corr = [];
        for sess1 = 1:6
            current_sess1 = nanmean(current_mouse(:,:,[1:5]+5*(sess1-1)),3);
            for sess2 = 1:6
                current_sess2 = nanmean(current_mouse(:,:,[1:5]+5*(sess2-1)),3);
                valid_cells = nanmean(current_sess1,2) > 0 & nanmean(current_sess2,2) > 0;
                
                rate_corr(sess1,sess2) = corr(nanmean(current_sess1(valid_cells,:),2),nanmean(current_sess2(valid_cells,:),2),'rows','pairwise');
                pv_corr(sess1,sess2) = nanmean(diag(corr(current_sess1(valid_cells,:),current_sess2(valid_cells,:),'rows','pairwise')));
                tuning_corr(sess1,sess2) = nanmedian(diag(corr(current_sess1(valid_cells,:)',current_sess2(valid_cells,:)','rows','pairwise')));
            end
        end
        
        rate_temp = [rate_corr(1:2,3:6),rate_corr(3:4,5:6)];
        pv_temp = [pv_corr(1:2,3:6),pv_corr(3:4,5:6)];
        tuning_temp = [tuning_corr(1:2,3:6),tuning_corr(3:4,5:6)];
        
        mdl = fitlm([rate_temp(:),tuning_temp(:)],pv_temp(:));
        
        pv_rate_tuning_relationships(mouse,1) = mdl.Rsquared.Ordinary;
        pv_rate_tuning_relationships(mouse,2) = corr(rate_temp(:),pv_temp(:),'rows','pairwise').^2;
        pv_rate_tuning_relationships(mouse,3) = corr(pv_temp(:),tuning_temp(:),'rows','pairwise').^2;
        pv_rate_tuning_relationships(mouse,4) = corr(rate_temp(:),tuning_temp(:),'rows','pairwise').^2;
        
        
    end
    pv_rate_tuning_relationships_areas{area} = pv_rate_tuning_relationships;
    
end


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
nat_movie = 1;
cell_cutoff = 15;
num_repeats = 30;
for area = 1:6
    valid_mice = neuropixels_cell_count(:,area,nat_movie) >= cell_cutoff & movie_repeats(:,nat_movie) == num_repeats;
    current_area = neuropixels_population_vectors(valid_mice,area,nat_movie);
    
    all_mice_relationships = [];
    for mouse = 1:size(current_area,1)
        
        clc;
        disp(['Calculating the relationship between activity rate,'])
        disp(['absolute activity rate difference or absolute activity rate difference score'])
        disp(['and tuning curve correlation across movie repeats within a block:'])
        disp(['Dataset: Neuropixels | Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        
        current_mouse = current_area{mouse};
        
        for repeat1 = 1:30
            for repeat2 = 1:30
                tuning_corr_blockA = diag(corr(current_mouse(:,:,repeat1)',current_mouse(:,:,repeat2)'));
                tuning_corr_blockB = diag(corr(current_mouse(:,:,repeat1+num_repeats)',current_mouse(:,:,repeat2+num_repeats)'));
                
                activity_rate_blockA = nanmean(nanmean(current_mouse(:,:,[repeat1,repeat2]),3),2);
                activity_rate_blockB = nanmean(nanmean(current_mouse(:,:,[repeat1+num_repeats,repeat2+num_repeats]),3),2);
                
                activity_diff_blockA = abs(nanmean(current_mouse(:,:,repeat1),2) - nanmean(current_mouse(:,:,repeat2),2));
                activity_diff_blockB =  abs(nanmean(current_mouse(:,:,repeat1+num_repeats),2) - nanmean(current_mouse(:,:,repeat2+num_repeats),2));
                
                
                
                activity_diff_score_blockA = abs([nanmean(current_mouse(:,:,repeat1),2) - nanmean(current_mouse(:,:,repeat2),2)]./...
                    [nanmean(current_mouse(:,:,repeat1),2) + nanmean(current_mouse(:,:,repeat2),2)]);
                activity_diff_score_blockB = abs([nanmean(current_mouse(:,:,repeat1+num_repeats),2) - nanmean(current_mouse(:,:,repeat2+num_repeats),2)]./...
                    [nanmean(current_mouse(:,:,repeat1+num_repeats),2) + nanmean(current_mouse(:,:,repeat2+num_repeats),2)]);
                
                
                
                activity_rate_VS_activity_diff_blockA(repeat1,repeat2) = corr(activity_rate_blockA,activity_diff_blockA,'rows','pairwise');
                activity_rate_VS_activity_diff_score_blockA(repeat1,repeat2) = corr(activity_rate_blockA,activity_diff_score_blockA,'rows','pairwise');
                activity_diff_VS_activity_diff_score_blockA(repeat1,repeat2) = corr(activity_diff_blockA,activity_diff_score_blockA,'rows','pairwise');
                
                tuning_VS_activity_rate_blockA(repeat1,repeat2) = corr(activity_rate_blockA,tuning_corr_blockA,'rows','pairwise');
                tuning_VS_activity_diff_blockA(repeat1,repeat2) = corr(activity_diff_blockA,tuning_corr_blockA,'rows','pairwise');
                tuning_VS_activity_diff_score_blockA(repeat1,repeat2) = corr(activity_diff_score_blockA,tuning_corr_blockA,'rows','pairwise');
                
                
                activity_rate_VS_activity_diff_blockB(repeat1,repeat2) = corr(activity_rate_blockB,activity_diff_blockB,'rows','pairwise');
                activity_rate_VS_activity_diff_score_blockB(repeat1,repeat2) = corr(activity_rate_blockB,activity_diff_score_blockB,'rows','pairwise');
                activity_diff_VS_activity_diff_score_blockB(repeat1,repeat2) = corr(activity_diff_blockB,activity_diff_score_blockB,'rows','pairwise');
                
                tuning_VS_activity_rate_blockB(repeat1,repeat2) = corr(activity_rate_blockB,tuning_corr_blockB,'rows','pairwise');
                tuning_VS_activity_diff_blockB(repeat1,repeat2) = corr(activity_diff_blockB,tuning_corr_blockB,'rows','pairwise');
                tuning_VS_activity_diff_score_blockB(repeat1,repeat2) = corr(activity_diff_score_blockB,tuning_corr_blockB,'rows','pairwise');
                
            end
        end
        
        triu_ind = boolean(triu(ones(30),1));
        
        current_mouse_relationships = [[tuning_VS_activity_rate_blockA(triu_ind);tuning_VS_activity_rate_blockB(triu_ind)],...
            [tuning_VS_activity_diff_blockA(triu_ind);tuning_VS_activity_diff_blockB(triu_ind)],...
            [tuning_VS_activity_diff_score_blockA(triu_ind);tuning_VS_activity_diff_score_blockB(triu_ind)]].^2;
        
        all_mice_relationships = [all_mice_relationships;nanmean(current_mouse_relationships)];
        
    end
    
    all_mice_relationships_areas{area} = all_mice_relationships;
    
end


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
cell_cutoff = 20;
for area = 1:6
    valid_mice = min(calcium_excitatory_cell_count{area},[],2) >= cell_cutoff;
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,1);
    
    
    all_relationships = [];
    for mouse = 1:length(current_area)
        clc;
        disp(['Calculating the relationship between activity rate,'])
        disp(['absolute activity rate difference or absolute activity rate '])
        disp(['difference score and tuning curve correlation across days:'])
        disp(['Dataset: Calcium imaging | Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse};
        
        tuning_vs_activity_rate = [];
        tuning_rate_vs_activity_diff = [];
        tuning_rate_vs_activity_diff = [];
        
        for sess1 = 1:6
            current_sess1 = nanmean(current_mouse(:,:,[1:5]+5*(sess1-1)),3);
            for sess2 = 1:6
                current_sess2 = nanmean(current_mouse(:,:,[1:5]+5*(sess2-1)),3);
                valid_cells = nanmean(current_sess1,2) > 0 & nanmean(current_sess2,2) > 0;
                
                activity_rate_both_sess = [nanmean(current_sess1(valid_cells,:),2),nanmean(current_sess2(valid_cells,:),2)];
                mean_activity_rate = nanmean(activity_rate_both_sess,2);
                abs_activity_diff = abs(activity_rate_both_sess(:,1) - activity_rate_both_sess(:,2));
                abs_activity_dff_score =  abs([activity_rate_both_sess(:,1) - activity_rate_both_sess(:,2)] ./ [activity_rate_both_sess(:,1) + activity_rate_both_sess(:,2)]);
                tuning_corr = diag(corr(current_sess1(valid_cells,:)',current_sess2(valid_cells,:)','rows','pairwise'));
                
                tuning_vs_activity_rate(sess1,sess2) = corr(mean_activity_rate,tuning_corr,'rows','pairwise');
                tuning_rate_vs_activity_diff(sess1,sess2) = corr(tuning_corr,abs_activity_diff,'rows','pairwise');
                tuning_rate_vs_activity_diff(sess1,sess2) = corr(tuning_corr,abs_activity_diff,'rows','pairwise');
                
            end
        end
        triu_id = boolean(triu(ones(6),1));
        relationships = [tuning_vs_activity_rate(triu_id),...
            tuning_rate_vs_activity_diff(triu_id),tuning_rate_vs_activity_diff(triu_id)].^2;
        all_relationships(mouse,:) = nanmean(relationships);
    end
    all_relationships_areas{area} = all_relationships;
end

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
neuropixels_elapse_repeat_decoder_acc_areas = {};
neuropixels_elapse_repeat_decoder_acc_areas_shuffle = {};
calcium_elapse_repeat_decoder_acc_areas = {};
calcium_elapse_repeat_decoder_acc_areas_shuffle = {};

cell_cutoff = 15;
nat_movie = 1;
num_repeats = 30;
for area = 1:6
    valid_mice = neuropixels_cell_count(:,area,nat_movie) >= cell_cutoff & movie_repeats(:,nat_movie) == num_repeats;
    current_area = neuropixels_population_vectors(valid_mice,area,nat_movie);
    
    elapse_repeat_decoder_acc = [];
    elapse_repeat_decoder_acc_shuffle  = [];
    for mouse = 1:size(current_area,1)
        
        current_mouse = current_area{mouse};
        current_mouse_blockA = current_mouse(:,:,1:30);
        current_mouse_blockB = current_mouse(:,:,31:60);
        
        mean_knn_decoding_accuracy = nan(30);
        mean_knn_decoding_accuracy_shuffle  = nan(30);
        for repeat1 = 1:num_repeats
            for repeat2 = 1:num_repeats
                if repeat1 < repeat2
                    clc;
                    disp(['Performing time-lapse decoding across repeats:'])
                    disp(['Dataset: Neuropixls | Stimulus: Natural movie 1 | Area: ',brain_areas{area},...
                        ' | Mouse: ',num2str(mouse),'/',num2str(length(current_area)),' | Movie repeat: ',num2str(repeat1),'/29'])
                    
                    % real data
                    train_data = current_mouse_blockA(:,:,repeat1);
                    test_data = current_mouse_blockA(:,:,repeat2);
                    
                    mdl = fitcknn(train_data',[1:30]);
                    prediction = [];
                    prediction(1,:) = predict(mdl,test_data');
                    
                    train_data = current_mouse_blockB(:,:,repeat1);
                    test_data = current_mouse_blockB(:,:,repeat2);
                    mdl = fitcknn(train_data',[1:30]);
                    prediction(2,:) = predict(mdl,test_data');
                    
                    accuracy = (sum(prediction == [1:30],2)./30)*100;
                    mean_knn_decoding_accuracy(repeat1,repeat2) = nanmean(accuracy);
                    
                    % shuffled data
                    train_data = current_mouse_blockA(:,:,repeat1);
                    test_data = current_mouse_blockA(:,:,repeat2);
                    
                    mdl = fitcknn(train_data',randperm(30));
                    prediction = [];
                    prediction(1,:) = predict(mdl,test_data');
                    
                    train_data = current_mouse_blockB(:,:,repeat1);
                    test_data = current_mouse_blockB(:,:,repeat2);
                    mdl = fitcknn(train_data',randperm(30));
                    prediction(2,:) = predict(mdl,test_data');
                    
                    accuracy = (sum(prediction == [1:30],2)./30)*100;
                    mean_knn_decoding_accuracy_shuffle(repeat1,repeat2) = nanmean(accuracy);
                end
                
            end
        end
        
        for diagonal = 1:29
            elapse_repeat_decoder_acc(mouse,diagonal) = nanmean(diag(mean_knn_decoding_accuracy,diagonal));
            elapse_repeat_decoder_acc_shuffle(mouse,diagonal) = nanmean(diag(mean_knn_decoding_accuracy_shuffle,diagonal));
        end
        
    end
    neuropixels_elapse_repeat_decoder_acc_areas{area} = elapse_repeat_decoder_acc;
    neuropixels_elapse_repeat_decoder_acc_areas_shuffle{area} = elapse_repeat_decoder_acc_shuffle;
end

nat_movie = 1;
cell_cutoff = 20;
for area = 1:6
    valid_mice = min(calcium_excitatory_cell_count{area},[],2) >= cell_cutoff;
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie);
    
    
    elapse_repeat_decoder_acc = [];
    elapse_repeat_decoder_acc_shuffle = [];
    for mouse = 1:length(current_area)
        
        current_mouse = current_area{mouse};
        
        current_mouse_sess1 =  current_mouse(:,:,1:10);
        valid_cells = nanmean(squeeze(nanmean(current_mouse_sess1,2)),2)>0;
        valid_current_mouse_sess1 = current_mouse_sess1(valid_cells,:,:);
        
        current_mouse_sess2 =  current_mouse(:,:,11:20);
        valid_cells = nanmean(squeeze(nanmean(current_mouse_sess2,2)),2)>0;
        valid_current_mouse_sess2 = current_mouse_sess2(valid_cells,:,:);
        
        current_mouse_sess3 =  current_mouse(:,:,21:30);
        valid_cells = nanmean(squeeze(nanmean(current_mouse_sess3,2)),2)>0;
        valid_current_mouse_sess3 = current_mouse_sess3(valid_cells,:,:);
        
        knn_decoder_across_sessions = nan(10);
        knn_decoder_across_sessions_shuffle = nan(10);
        for repeat1 = 1:10
            for repeat2 = 1:10
                if repeat1 < repeat2
                    clc;
                    disp(['Performing time-lapse decoding across repeats:'])
                    disp(['Dataset: Neuropixls | Stimulus: Natural movie 1 | Area: ',brain_areas{area},...
                        ' | Mouse: ',num2str(mouse),'/',num2str(length(current_area)),' |  Movie repeat: ',num2str(repeat1),'/9'])
                    
                    % real data
                    train_data = valid_current_mouse_sess1(:,:,repeat1);
                    test_data = valid_current_mouse_sess1(:,:,repeat2);
                    mdl = fitcknn(train_data',[1:30]);
                    pred = predict(mdl,test_data');
                    knn_decoder_across_sessions(repeat1,repeat2,1) = (sum(pred == [1:30]')./30)*100;
                    
                    train_data = valid_current_mouse_sess2(:,:,repeat1);
                    test_data = valid_current_mouse_sess2(:,:,repeat2);
                    mdl = fitcknn(train_data',[1:30]);
                    pred = predict(mdl,test_data');
                    knn_decoder_across_sessions(repeat1,repeat2,2) = (sum(pred == [1:30]')./30)*100;
                    
                    train_data = valid_current_mouse_sess3(:,:,repeat1);
                    test_data = valid_current_mouse_sess3(:,:,repeat2);
                    mdl = fitcknn(train_data',[1:30]);
                    pred = predict(mdl,test_data');
                    knn_decoder_across_sessions(repeat1,repeat2,3) = (sum(pred == [1:30]')./30)*100;
                    
                    % shuffle data
                    train_data = valid_current_mouse_sess1(:,:,repeat1);
                    test_data = valid_current_mouse_sess1(:,:,repeat2);
                    mdl = fitcknn(train_data',randperm(30));
                    pred = predict(mdl,test_data');
                    knn_decoder_across_sessions_shuffle(repeat1,repeat2,1) = (sum(pred == [1:30]')./30)*100;
                    
                    train_data = valid_current_mouse_sess2(:,:,repeat1);
                    test_data = valid_current_mouse_sess2(:,:,repeat2);
                    mdl = fitcknn(train_data',randperm(30));
                    pred = predict(mdl,test_data');
                    knn_decoder_across_sessions_shuffle(repeat1,repeat2,2) = (sum(pred == [1:30]')./30)*100;
                    
                    train_data = valid_current_mouse_sess3(:,:,repeat1);
                    test_data = valid_current_mouse_sess3(:,:,repeat2);
                    mdl = fitcknn(train_data',randperm(30));
                    pred = predict(mdl,test_data');
                    knn_decoder_across_sessions_shuffle(repeat1,repeat2,3) = (sum(pred == [1:30]')./30)*100;
                end
            end
        end
        
        mean_acc_across_mice = nanmean(knn_decoder_across_sessions,3);
        mean_acc_across_mice_shuffle = nanmean(knn_decoder_across_sessions_shuffle,3);
        
        for diagonal = 1:9
            elapse_repeat_decoder_acc(mouse,diagonal) = nanmean(diag(mean_acc_across_mice,diagonal));
            elapse_repeat_decoder_acc_shuffle(mouse,diagonal) = nanmean(diag(mean_acc_across_mice_shuffle,diagonal));
        end
        
    end
    calcium_elapse_repeat_decoder_acc_areas{area} = elapse_repeat_decoder_acc;
    calcium_elapse_repeat_decoder_acc_areas_shuffle{area} = elapse_repeat_decoder_acc_shuffle;
end

p =[];
z = [];
plt = [];
figure('units','normalized','position',[0.3 0.2 0.2 0.5])

for area = 1:6
    % neuropixels data
    neuropixels_current_area = neuropixels_elapse_repeat_decoder_acc_areas{area};
    neuropixels_mean_acc = nanmean(neuropixels_current_area);
    neuropixels_std_acc = nanstd(neuropixels_current_area);
    neuropixels_ste_acc = neuropixels_std_acc./sqrt(size(neuropixels_current_area,1));
    
    neuropixels_current_area_shuffle = neuropixels_elapse_repeat_decoder_acc_areas_shuffle{area};
    neuropixels_mean_acc_shuffle = nanmean(neuropixels_current_area_shuffle);
    neuropixels_std_acc_shuffle = nanstd(neuropixels_current_area_shuffle);
    neuropixels_ste_acc_shuffle = neuropixels_std_acc_shuffle./sqrt(size(neuropixels_current_area_shuffle,1));
    
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
    
    
    % calcium imaging data
    calcium_current_area = calcium_elapse_repeat_decoder_acc_areas{area};
    calcium_mean_acc = nanmean(calcium_current_area);
    calcium_std_acc = nanstd(calcium_current_area);
    calcium_ste_acc = calcium_std_acc./sqrt(size(calcium_current_area,1));
    
    calcium_current_area_shuffle = calcium_elapse_repeat_decoder_acc_areas_shuffle{area};
    calcium_mean_acc_shuffle = nanmean(calcium_current_area_shuffle);
    calcium_std_acc_shuffle = nanstd(calcium_current_area_shuffle);
    calcium_ste_acc_shuffle = calcium_std_acc_shuffle./sqrt(size(calcium_current_area_shuffle,1));
    
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
calcium_elapse_repeat_areas = {};
neuropixels_elapse_repeat_areas = {};

cell_cutoff = 15;
nat_movie = 1;
num_repeats = 30;
for area = 1:6
    valid_mice = neuropixels_cell_count(:,area,nat_movie) >= cell_cutoff & movie_repeats(:,nat_movie) == num_repeats;
    current_area = neuropixels_population_vectors(valid_mice,area,nat_movie);
    
    elapse_repeat_pv = [];
    elapse_repeat_rate  = [];
    elapse_repeat_tuning  = [];
    for mouse = 1:size(current_area,1)
        clc;
        disp(['Calculating PV, ensemble rate and tuning curve correlations across movie repeats:'])
        disp(['Dataset: Neuropixls | Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse};
        current_mouse_blockA = current_mouse(:,:,1:30);
        current_mouse_blockB = current_mouse(:,:,31:60);
        
        rate_corr = [];
        tuning_corr = [];
        for repeat1 = 1:num_repeats
            for repeat2 = 1:num_repeats
                
                pv_corr(repeat1,repeat2,1) = nanmean(diag(corr(current_mouse_blockA(:,:,repeat1),current_mouse_blockA(:,:,repeat2))));
                pv_corr(repeat1,repeat2,2) = nanmean(diag(corr(current_mouse_blockB(:,:,repeat1),current_mouse_blockB(:,:,repeat2))));
                
                rate_corr(repeat1,repeat2,1) = corr(nanmean(current_mouse_blockA(:,:,repeat1),2),nanmean(current_mouse_blockA(:,:,repeat2),2));
                rate_corr(repeat1,repeat2,2) = corr(nanmean(current_mouse_blockB(:,:,repeat1),2),nanmean(current_mouse_blockB(:,:,repeat2),2));
                
                tuning_corr(repeat1,repeat2,1) = nanmedian(diag(corr(current_mouse_blockA(:,:,repeat1)',current_mouse_blockA(:,:,repeat2)')));
                tuning_corr(repeat1,repeat2,2) = nanmedian(diag(corr(current_mouse_blockB(:,:,repeat1)',current_mouse_blockB(:,:,repeat2)')));
                
            end
        end
        pv_corr = nanmean(pv_corr,3);
        rate_corr = nanmean(rate_corr,3);
        tuning_corr = nanmean(tuning_corr,3);
        
        for diagonal = 1:29
            elapse_repeat_pv(mouse,diagonal) = nanmean(diag(pv_corr,diagonal));
            elapse_repeat_rate(mouse,diagonal) = nanmean(diag(rate_corr,diagonal));
            elapse_repeat_tuning(mouse,diagonal) = nanmean(diag(tuning_corr,diagonal));
        end
        
    end
    neuropixels_elapse_repeat_areas(area,:) = {elapse_repeat_pv,elapse_repeat_rate,elapse_repeat_tuning};
end


nat_movie = 1;
cell_cutoff = 20;
for area = 1:6
    valid_mice = min(calcium_excitatory_cell_count{area},[],2) >= cell_cutoff;
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie);
    
    
    elapse_repeat_pv = [];
    elapse_repeat_rate = [];
    elapse_repeat_tuning = [];
    for mouse = 1:length(current_area)
        clc;
        disp(['Calculating PV, ensemble rate and tuning curve correlations across movie repeats:'])
        disp(['Dataset: Calcium imaging | Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse};
        
        current_mouse_sess1 =  current_mouse(:,:,1:10);
        valid_cells = nanmean(squeeze(nanmean(current_mouse_sess1,2)),2)>0;
        valid_current_mouse_sess1 = current_mouse_sess1(valid_cells,:,:);
        
        current_mouse_sess2 =  current_mouse(:,:,11:20);
        valid_cells = nanmean(squeeze(nanmean(current_mouse_sess2,2)),2)>0;
        valid_current_mouse_sess2 = current_mouse_sess2(valid_cells,:,:);
        
        current_mouse_sess3 =  current_mouse(:,:,21:30);
        valid_cells = nanmean(squeeze(nanmean(current_mouse_sess3,2)),2)>0;
        valid_current_mouse_sess3 = current_mouse_sess3(valid_cells,:,:);
        
        pv_corr = [];
        rate_corr = [];
        tuning_corr = [];
        for repeat1 = 1:10
            for repeat2 = 1:10
                pv_corr(repeat1,repeat2,1) = nanmean(diag(corr(current_mouse_sess1(:,:,repeat1),current_mouse_sess1(:,:,repeat2))));
                pv_corr(repeat1,repeat2,2) = nanmean(diag(corr(current_mouse_sess2(:,:,repeat1),current_mouse_sess2(:,:,repeat2))));
                pv_corr(repeat1,repeat2,3) = nanmean(diag(corr(current_mouse_sess3(:,:,repeat1),current_mouse_sess3(:,:,repeat2))));
                
                rate_corr(repeat1,repeat2,1) = corr(nanmean(current_mouse_sess1(:,:,repeat1),2),nanmean(current_mouse_sess1(:,:,repeat2),2));
                rate_corr(repeat1,repeat2,2) = corr(nanmean(current_mouse_sess2(:,:,repeat1),2),nanmean(current_mouse_sess2(:,:,repeat2),2));
                rate_corr(repeat1,repeat2,3) = corr(nanmean(current_mouse_sess3(:,:,repeat1),2),nanmean(current_mouse_sess3(:,:,repeat2),2));
                
                tuning_corr(repeat1,repeat2,1) = nanmean(diag(corr(current_mouse_sess1(:,:,repeat1)',current_mouse_sess1(:,:,repeat2)')));
                tuning_corr(repeat1,repeat2,2) = nanmean(diag(corr(current_mouse_sess2(:,:,repeat1)',current_mouse_sess2(:,:,repeat2)')));
                tuning_corr(repeat1,repeat2,3) = nanmean(diag(corr(current_mouse_sess3(:,:,repeat1)',current_mouse_sess3(:,:,repeat2)')));
                
            end
        end
        pv_corr = nanmean(pv_corr,3);
        rate_corr = nanmean(rate_corr,3);
        tuning_corr = nanmean(tuning_corr,3);
        
        
        for diagonal = 1:9
            elapse_repeat_pv(mouse,diagonal) = nanmean(diag(pv_corr,diagonal));
            elapse_repeat_rate(mouse,diagonal) = nanmean(diag(rate_corr,diagonal));
            elapse_repeat_tuning(mouse,diagonal) = nanmean(diag(tuning_corr,diagonal));
        end
        
    end
    calcium_elapse_repeat_areas(area,:) = {elapse_repeat_pv,elapse_repeat_rate,elapse_repeat_tuning};
end

ylims = {[-0.095 0],[-0.1 0],[-0.12 0];...
    [-0.06 0],[-0.125 0],[-0.025 0.005]};
figure('units','normalized','position',[0.2 0.2 0.475 0.45])
for measurment = 1:3
    for area = 1:6
        % neuropixels data
        neuropixels_current_area = neuropixels_elapse_repeat_areas{area,measurment};
        neuropixels_current_area_diff = neuropixels_current_area-neuropixels_current_area(:,1);
        neuropixels_mean_stability = nanmean(neuropixels_current_area_diff);
        neuropixels_std_stability = nanstd(neuropixels_current_area_diff);
        neuropixels_ste_stability = neuropixels_std_stability./sqrt(size(neuropixels_current_area,1));
        
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
        
        
        
        % calcium imaging data
        
        calcium_current_area = calcium_elapse_repeat_areas{area,measurment};
        calcium_current_area_diff = calcium_current_area-calcium_current_area(:,1);
        calcium_mean_stability = nanmean(calcium_current_area_diff);
        calcium_std_stability = nanstd(calcium_current_area_diff);
        calcium_ste_stability = calcium_std_stability./sqrt(size(calcium_current_area_diff,1));
        
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
repeat_num = 30;
nat_movie = 1;
valid_mice = neuropixels_running_speed(movie_repeats(:,nat_movie) == repeat_num,nat_movie);

running_speed_across_mice = [];
for mouse = 1:size(valid_mice,1)
    running_speed_across_mice(:,:,mouse) = valid_mice{mouse,nat_movie};
end

mean_running_speed_across_blocks = squeeze(nanmean(running_speed_across_mice,2))';
mean_running_speed = nanmean(mean_running_speed_across_blocks);
std_running_speed = nanstd(mean_running_speed_across_blocks);
ste_running_speed = std_running_speed./sqrt(size(mean_running_speed_across_blocks,1));

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

sig_mat = [];
for trial1 = 1:30
    for trial2 = 1:30
        sig_mat(trial1,trial2) = ttest(mean_running_speed_across_blocks(:,trial1),mean_running_speed_across_blocks(:,trial2));
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
repeat_num = 30;
nat_movie = 1;
valid_mice = neuropixels_pupil_size(movie_repeats(:,nat_movie) == repeat_num,nat_movie);

pupil_area_across_mice = [];
sub = 1;
for mouse = 1:size(valid_mice,1)
    if ~isempty(valid_mice{mouse,nat_movie})
        pupil_area_across_mice(:,:,sub) = valid_mice{mouse,nat_movie};
        sub = sub + 1;
    end
end

mean_pupil_area_across_blocks = squeeze(nanmean(pupil_area_across_mice,2))';
mean_pupil_area = nanmean(mean_pupil_area_across_blocks);
std_pupil_area = nanstd(mean_pupil_area_across_blocks);
ste_pupil_area = std_pupil_area./sqrt(size(mean_pupil_area_across_blocks,1));

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

sig_mat = [];
for trial1 = 1:30
    for trial2 = 1:30
        sig_mat(trial1,trial2) = ttest(mean_pupil_area_across_blocks(:,trial1),mean_pupil_area_across_blocks(:,trial2));
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
cell_cutoff = 15;
nat_movie = 1;
num_repeats = 30;
mean_activity_all_areas = {};
for area = 1:6
    valid_mice = neuropixels_cell_count(:,area,nat_movie) >= cell_cutoff & movie_repeats(:,nat_movie) == num_repeats;
    current_area = neuropixels_population_vectors(valid_mice,area,nat_movie);
    
    mean_activity_all_mice = [];
    for mouse = 1:size(current_area,1)
        current_mouse = current_area{mouse}.*30;
        
        mean_activity_blockA = squeeze(nanmean(nanmean(current_mouse(:,:,1:30),2)));
        mean_activity_blockB = squeeze(nanmean(nanmean(current_mouse(:,:,31:60),2)));
        mean_activity_all_mice(mouse,:) = nanmean([mean_activity_blockA,mean_activity_blockB]');
    end
    
    mean_activity_all_areas{area} = mean_activity_all_mice;
end


figure('units','normalized','position',[0.3 0.3 0.35 0.375])
for area = 1:6
    current_area =  mean_activity_all_areas{area};
    mean_activity = nanmean(current_area);
    std_activity = std(current_area);
    ste_activity = std_activity./sqrt(size(current_area,1));
    
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
cell_cutoff = 15;
nat_movie = 1;
num_repeats = 30;
repeat_lim = {1:30,9:30};
subset_list = {'Full','Repeats 9-30'};
elapse_repeat_rate_corr_area = {};
for area = 1:6
    valid_mice = neuropixels_cell_count(:,area,nat_movie) >= cell_cutoff & movie_repeats(:,nat_movie) == num_repeats;
    current_area = neuropixels_population_vectors(valid_mice,area,nat_movie);
    for repeat_span = 1:2
        current_span = repeat_lim{repeat_span};
        
        elapse_repeat_rate_corr = [];
        for mouse = 1:size(current_area,1)
            
            clc;
            disp(['Calculating ensemble rate correlation between movie repeats:'])
            disp(['Dataset: Neuropixels | Stimulus: Natural movie 1 | Area: ',brain_areas{area},...
                ' | Subset: ',subset_list{repeat_span},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
            
            current_mouse = current_area{mouse};
            
            current_mouse_blockA = current_mouse(:,:,1:30);
            current_mouse_blockA = squeeze(nanmean(current_mouse_blockA(:,:,current_span),2));
            
            current_mouse_blockB = current_mouse(:,:,31:60);
            current_mouse_blockB = squeeze(nanmean(current_mouse_blockB(:,:,current_span),2));
            
            mean_rate_corr = [];
            mean_rate_corr(:,:,1) = corr(current_mouse_blockA);
            mean_rate_corr(:,:,2) = corr(current_mouse_blockB);
            mean_rate_corr = nanmean(mean_rate_corr,3);
            
            for diagonal = 1:length(current_span)-1
                elapse_repeat_rate_corr(mouse,diagonal) = nanmean(diag(mean_rate_corr,diagonal));
            end
            
        end
        elapse_repeat_rate_corr_area{repeat_span,area} = elapse_repeat_rate_corr;
    end
end

pvalues = [];
df = [];
chi = [];
ylims = [0.88 0.98;0.87 1;0.86 0.98;0.88 0.98;0.87 0.98;0.87 0.98];
figure('units','normalized','position',[0.3 0.3 0.35 0.385])
for area = 1:6
    for repeat_span = 1:2
        current_area = elapse_repeat_rate_corr_area{repeat_span,area};
        
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

corrected_pval = bonf_holm(pvalues(2,:));

VarNames = {'area','df','chi','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas(1:6)'),df(2,:)',chi(2,:)',pvalues(2,:)',corrected_pval(:),'VariableNames',VarNames);

clc;
disp(['Ensemble rate correlation as a function of elapsed time for subsampled data (repeats 9-30)'])
disp(['Friedman’s tests with Holm–Bonferroni correction:'])
disp(statistics)

%% Figure S2E - Tuning curve correlation on subset of movie repeats (repeats 9-30)
cell_cutoff = 15;
nat_movie = 1;
num_repeats = 30;
repeat_lim = {1:30,9:30};
subset_list = {'Full','Repeats 9-30'};
elapse_repeat_tuning_corr_area = {};
for area = 1:6
    valid_mice = neuropixels_cell_count(:,area,nat_movie) >= cell_cutoff & movie_repeats(:,nat_movie) == num_repeats;
    current_area = neuropixels_population_vectors(valid_mice,area,nat_movie);
    for repeat_span = 1:2
        current_span = repeat_lim{repeat_span};
        
        elapse_repeat_tuning_corr = [];
        for mouse = 1:size(current_area,1)
            
            clc;
            disp(['Calculating tuning curve correlation between movie repeats:'])
            disp(['Dataset: Neuropixels | Stimulus: Natural movie 1 | Area: ',brain_areas{area},...
                ' | Subset: ',subset_list{repeat_span},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
            
            current_mouse = current_area{mouse};
            
            current_mouse_blockA = current_mouse(:,:,1:30);
            current_mouse_blockB = current_mouse(:,:,31:60);
            
            tuning_corr = [];
            for repeat1 = 1:length(current_span)
                for repeat2 = 1:length(current_span)
                    tuning_corr(repeat1,repeat2,1) = nanmedian(diag(corr(current_mouse_blockA(:,:,repeat1)',current_mouse_blockA(:,:,repeat2)')));
                    tuning_corr(repeat1,repeat2,2) = nanmedian(diag(corr(current_mouse_blockB(:,:,repeat1)',current_mouse_blockB(:,:,repeat2)')));
                end
            end
            mean_tuning_corr = nanmean(tuning_corr,3);
            
            for diagonal = 1:length(current_span)-1
                elapse_repeat_tuning_corr(mouse,diagonal) = nanmean(diag(mean_tuning_corr,diagonal));
            end
            
        end
        elapse_repeat_tuning_corr_area{repeat_span,area} = elapse_repeat_tuning_corr;
    end
end

pvalues = [];
df = [];
chi = [];
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
disp(['Friedman’s tests with Holm–Bonferroni correction:'])
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
disp(['Friedman’s tests with Holm–Bonferroni correction:'])
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
disp(['Friedman’s tests with Holm–Bonferroni correction:'])
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

nat_movie = 2;
area = 6;
cell_cutoff = 20;
valid_mice = calcium_excitatory_cell_count{area}(:,1) >= cell_cutoff;
current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie);

pv_corr_across_mice = [];
for mouse = 1:length(current_area)
    current_mouse = current_area{mouse};
    
    current_mouse_blockA =  current_mouse(:,:,1:5);
    current_mouse_blockA_half1 = nanmean(current_mouse_blockA(:,:,1:2),3);
    current_mouse_blockA_half2 = nanmean(current_mouse_blockA(:,:,3:5),3);
    
    current_mouse_blockB =  current_mouse(:,:,6:10);
    current_mouse_blockB_half1 = nanmean(current_mouse_blockB(:,:,1:2),3);
    current_mouse_blockB_half2 = nanmean(current_mouse_blockB(:,:,3:5),3);
    
    pv_corr_across_mice(:,:,mouse) = corr([current_mouse_blockA_half1,current_mouse_blockA_half2,...
        current_mouse_blockB_half1,current_mouse_blockB_half2]);
    
end

mean_pv_corr_across_mice = nanmean(pv_corr_across_mice,3);
mean_pv_corr_across_mice(boolean(eye(size(mean_pv_corr_across_mice,1)))) = NaN;
mean_pv_corr_across_mice(isnan(mean_pv_corr_across_mice)) = max(mean_pv_corr_across_mice(:));

figure('units','normalized','position',[0.3 0.3 0.275 0.4])
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
nat_movie = 2;
cell_cutoff = 20;
rate_corr_between_blocks_across_areas = {};
for area = 1:6
    valid_mice = calcium_excitatory_cell_count{area}(:,1) >= cell_cutoff;
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie);
    
    rate_corr_across_blocks = [];
    for mouse = 1:length(current_area)
        
        clc;
        disp(['Calculating ensemble rate correlation between blocks:'])
        disp(['Dataset: Calcium imaging | Stimulus: Natural movie 3 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse};
        current_mouse_blockA =  current_mouse(:,:,1:5);
        current_mouse_blockA_half1 = nanmean(nanmean(current_mouse_blockA(:,:,1:2),3),2);
        current_mouse_blockA_half2 = nanmean(nanmean(current_mouse_blockA(:,:,3:5),3),2);
        
        current_mouse_blockB =  current_mouse(:,:,6:10);
        current_mouse_blockB_half1 = nanmean(nanmean(current_mouse_blockB(:,:,1:2),3),2);
        current_mouse_blockB_half2 = nanmean(nanmean(current_mouse_blockB(:,:,3:5),3),2);
        
        rate_corr = corr([current_mouse_blockA_half1,current_mouse_blockA_half2,...
            current_mouse_blockB_half1,current_mouse_blockB_half2]);
        
        rate_corr_across_blocks(mouse,1) = nanmean([rate_corr(1,2),rate_corr(3,4)]);
        rate_corr_across_blocks(mouse,2) = nanmean(nanmean(rate_corr(1:2,3:4)));
        
    end
    rate_corr_between_blocks_across_areas{area} = rate_corr_across_blocks;
    
end


pvalue = [];
zvalue = [];
plt = [];
figure('Units','Normalized','Position',[0.2 0.4 0.2 0.325])
for area = 1:6
    current_area = rate_corr_between_blocks_across_areas{area};
    mean_stability = nanmean(current_area);
    std_stability = nanstd(current_area);
    ste_stability = std_stability./sqrt(size(current_area,1));
    
    [pvalue(area),~,stats] = signrank(current_area(:,1),current_area(:,2));
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

corrected_pvalue = bonf_holm(pvalue);

VarNames = {'area','zvalue','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas(1:6)'),zvalue(:),pvalue(:),corrected_pvalue(:),'VariableNames',VarNames);

clc;
disp(['Ensemble rate correlation within blocks compared to between blocks'])
disp(['Two-tailed Wilcoxon signed-rank tests with Holm–Bonferroni correction:'])
disp(statistics)

%% Figure S3C - Tuning curve correlation within and between blocks of natural movie 3 (calcium imaging)
nat_movie = 2;
cell_cutoff = 20;
tuning_corr_between_blocks_across_areas = {};
for area = 1:6
    valid_mice = calcium_excitatory_cell_count{area}(:,1) >= cell_cutoff;
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie);
    
    tuning_corr_across_blocks = [];
    for mouse = 1:length(current_area)
        
        clc;
        disp(['Calculating tuning curve correlation between blocks:'])
        disp(['Dataset: Calcium imaging | Stimulus: Natural movie 3 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse};
        current_mouse_halves = [];
        current_mouse_halves(:,:,1) = nanmean(current_mouse(:,:,1:2),3);
        current_mouse_halves(:,:,2) = nanmean(current_mouse(:,:,3:5),3);
        current_mouse_halves(:,:,3) = nanmean(current_mouse(:,:,6:7),3);
        current_mouse_halves(:,:,4) = nanmean(current_mouse(:,:,8:10),3);
        
        tuning_corr = [];
        for half1 = 1:4
            for half2 = 1:4
                tuning_corr(half1,half2) = nanmedian(diag(corr(current_mouse_halves(:,:,half1)',current_mouse_halves(:,:,half2)')));
            end
        end
        
        tuning_corr_across_blocks(mouse,1) = nanmean([tuning_corr(1,2),tuning_corr(3,4)]);
        tuning_corr_across_blocks(mouse,2) = nanmean(nanmean(tuning_corr(1:2,3:4)));
        
    end
    tuning_corr_between_blocks_across_areas{area} = tuning_corr_across_blocks;
    
end


pvalue = [];
zvalue = [];
plt = [];
figure('Units','Normalized','Position',[0.2 0.4 0.2 0.325])
for area = 1:6
    current_area = tuning_corr_between_blocks_across_areas{area};
    mean_stability = nanmean(current_area);
    std_stability = nanstd(current_area);
    ste_stability = std_stability./sqrt(size(current_area,1));
    
    [pvalue(area),~,stats] = signrank(current_area(:,1),current_area(:,2));
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

corrected_pvalue = bonf_holm(pvalue);

VarNames = {'area','zvalue','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas(1:6)'),zvalue(:),pvalue(:),corrected_pvalue(:),'VariableNames',VarNames);

clc;
disp(['Tuning curve correlation within blocks compared to between blocks'])
disp(['Two-tailed Wilcoxon signed-rank tests with Holm–Bonferroni correction:'])
disp(statistics)

%% Figure S3D - Ensmble rate and tuning curve correlation difference between blocks (calcium imaging)
nat_movie = 2;
cell_cutoff = 20;
rate_tuning_corr_diff_across_areas = {};
for area = 1:6
    valid_mice = calcium_excitatory_cell_count{area}(:,1) >= cell_cutoff;
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie);
    
    rate_corr_diff = [];
    tuning_corr_diff = [];
    for mouse = 1:length(current_area)
        
        clc;
        disp(['Calculating ensemble rate and tuning curve correlation between blocks:'])
        disp(['Dataset: Calcium imaging | Stimulus: Natural movie 3 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse};
        current_mouse_halves = [];
        current_mouse_halves(:,:,1) = nanmean(current_mouse(:,:,1:2),3);
        current_mouse_halves(:,:,2) = nanmean(current_mouse(:,:,3:5),3);
        current_mouse_halves(:,:,3) = nanmean(current_mouse(:,:,6:7),3);
        current_mouse_halves(:,:,4) = nanmean(current_mouse(:,:,8:10),3);
        
        rate_corr = [];
        tuning_corr = [];
        for half1 = 1:4
            for half2 = 1:4
                rate_corr(half1,half2) = corr(nanmean(current_mouse_halves(:,:,half1),2),nanmean(current_mouse_halves(:,:,half2),2));
                tuning_corr(half1,half2) = nanmedian(diag(corr(current_mouse_halves(:,:,half1)',current_mouse_halves(:,:,half2)')));
            end
        end
        
        tuning_corr_across_blocks = [];
        tuning_corr_across_blocks(1) = nanmean([tuning_corr(1,2),tuning_corr(3,4)]);
        tuning_corr_across_blocks(2) = nanmean(nanmean(tuning_corr(1:2,3:4)));
        
        rate_corr_across_blocks = [];
        rate_corr_across_blocks(1) = nanmean([rate_corr(1,2),rate_corr(3,4)]);
        rate_corr_across_blocks(2) = nanmean(nanmean(rate_corr(1:2,3:4)));
        
        rate_corr_diff(mouse) = rate_corr_across_blocks(1) - rate_corr_across_blocks(2);
        tuning_corr_diff(mouse) = tuning_corr_across_blocks(1) - tuning_corr_across_blocks(2);
    end
    
    rate_tuning_corr_diff_across_areas{area} = [rate_corr_diff',tuning_corr_diff'];
    
end


pvalue = [];
zvalue = [];
mean_stability = [];
ste_stability = [];

for area = 1:6
    current_area = rate_tuning_corr_diff_across_areas{area};
    mean_stability(area,:) = nanmean(current_area);
    std_stability = nanstd(current_area);
    ste_stability(area,:) = std_stability./sqrt(size(current_area,1));
    
    [pvalue(area,1),~,stats] = signrank(current_area(:,1));
    zvalue(area,1) = stats.zval;
    
    [pvalue(area,2),~,stats] = signrank(current_area(:,2));
    zvalue(area,2) = stats.zval;
end

plt = [];
figure('Units','Normalized','Position',[0.3 0.4 0.225 0.325])
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

corrected_pvalue = [];
corrected_pvalue(:,1) = bonf_holm(pvalue(:,1));
corrected_pvalue(:,2) = bonf_holm(pvalue(:,2));

VarNames = {'area','Rate_zvalue','Rate_pvalue','Rate_bonf_holm','Tuning_zvalue','Tuning_pvalue','Tuning_bonf_holm'};
statistics = table(cell2mat(brain_areas(1:6)'),zvalue(:,1),pvalue(:,1),corrected_pvalue(:,1),...
    zvalue(:,2),pvalue(:,2),corrected_pvalue(:,2),'VariableNames',VarNames);

clc;
disp(['Ensemble rate and tuning curve correlation within blocks compared to between blocks'])
disp(['Two-tailed Wilcoxon signed-rank tests with Holm–Bonferroni correction:'])
disp(statistics)

%% Figure S3E - Average ensemble rate correlation across V1 animals in the Brain Observatory group
cell_cutoff = 15;
nat_movie = 1;

reliability_cutoff = 0.6;
area = 1;
valid_mice = neuropixels_cell_count(:,area,nat_movie) >= cell_cutoff & movie_repeats(:,1) == 10;
current_area = squeeze(neuropixels_population_vectors(valid_mice,area,:));
across_blocks_ensembles_corr = [];
sub = 1;
for mouse = 1:length(current_area)
    current_mouse = current_area(mouse,:);
    
    movie1_blockA = current_mouse{1}(:,:,1:10);
    movie1_blockB = current_mouse{1}(:,:,11:20);
    
    movie3_blockA = current_mouse{2}(:,:,1:5);
    movie3_blockB = current_mouse{2}(:,:,6:10);
    
    mean_movie1_blockA = nanmean(movie1_blockA,3);
    mean_movie1_blockB = nanmean(movie1_blockB,3);
    
    mean_movie3_blockA = nanmean(movie3_blockA,3);
    mean_movie3_blockB = nanmean(movie3_blockB,3);
    
    
    movie1_blocks_stability = diag(corr(mean_movie1_blockA',mean_movie1_blockB'));
    movie3_blocks_stability = diag(corr(mean_movie3_blockA',mean_movie3_blockB'));
    
    valid_units_both_blocks = (movie1_blocks_stability >= reliability_cutoff) & (movie3_blocks_stability >= reliability_cutoff);
    if sum(valid_units_both_blocks) >= cell_cutoff
        valid_movie1_blockA_stability = movie1_blockA(valid_units_both_blocks,:,:);
        valid_movie1_blockB_stability = movie1_blockB(valid_units_both_blocks,:,:);
        valid_movie3_blockA_stability = movie3_blockA(valid_units_both_blocks,:,:);
        valid_movie3_blockB_stability = movie3_blockB(valid_units_both_blocks,:,:);
        
        mean_valid_movie1_blockA_stability = squeeze(nanmean(valid_movie1_blockA_stability,2));
        mean_valid_movie1_blockB_stability = squeeze(nanmean(valid_movie1_blockB_stability,2));
        mean_valid_movie3_blockA_stability = squeeze(nanmean(valid_movie3_blockA_stability,2));
        mean_valid_movie3_blockB_stability = squeeze(nanmean(valid_movie3_blockB_stability,2));
        
        mean_valid_movie1_half_blockA_stability = [nanmean(mean_valid_movie1_blockA_stability(:,1:5),2),...
            nanmean(mean_valid_movie1_blockA_stability(:,6:10),2)];
        mean_valid_movie1_half_blockB_stability = [nanmean(mean_valid_movie1_blockB_stability(:,1:5),2),...
            nanmean(mean_valid_movie1_blockB_stability(:,6:10),2)];
        
        mean_valid_movie3_half_blockA_stability = [nanmean(mean_valid_movie3_blockA_stability(:,1:2),2),...
            nanmean(mean_valid_movie3_blockA_stability(:,3:5),2)];
        mean_valid_movie3_half_blockB_stability = [nanmean(mean_valid_movie3_blockB_stability(:,1:2),2),...
            nanmean(mean_valid_movie3_blockB_stability(:,3:5),2)];
        
        
        
        across_blocks_ensembles = [mean_valid_movie3_half_blockA_stability,mean_valid_movie1_half_blockA_stability,...
            mean_valid_movie3_half_blockB_stability,mean_valid_movie1_half_blockB_stability];
        
        
        across_blocks_ensembles_corr(:,:,sub) = corr(across_blocks_ensembles);
        
        sub = sub + 1;
    end
end

figure('units','normalized','position',[0.3 0.3 0.25 0.35])
imagesc(nanmean(across_blocks_ensembles_corr,3))
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
cell_cutoff = 15;
nat_movie = 1;

reliability_cutoff = 0.6;
area = 1;
valid_mice = neuropixels_cell_count(:,area,nat_movie) >= cell_cutoff & movie_repeats(:,1) == 10;
current_area = squeeze(neuropixels_population_vectors(valid_mice,area,:));
movie1_blocks_stability_all_units = [];
movie3_blocks_stability_all_units = [];

for mouse = 1:length(current_area)
    current_mouse = current_area(mouse,:);
    
    movie1_blockA = current_mouse{1}(:,:,1:10);
    movie1_blockB = current_mouse{1}(:,:,11:20);
    
    movie3_blockA = current_mouse{2}(:,:,1:5);
    movie3_blockB = current_mouse{2}(:,:,6:10);
    
    mean_movie1_blockA = nanmean(movie1_blockA,3);
    mean_movie1_blockB = nanmean(movie1_blockB,3);
    
    mean_movie3_blockA = nanmean(movie3_blockA,3);
    mean_movie3_blockB = nanmean(movie3_blockB,3);
    
    
    movie1_blocks_stability = diag(corr(mean_movie1_blockA',mean_movie1_blockB'));
    movie3_blocks_stability = diag(corr(mean_movie3_blockA',mean_movie3_blockB'));
    
    movie1_blocks_stability_all_units = [movie1_blocks_stability_all_units;movie1_blocks_stability];
    movie3_blocks_stability_all_units = [movie3_blocks_stability_all_units;movie3_blocks_stability];
end

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

cell_cutoff = 15;
reliability_cutoff = 0.6;
for area = 1:6
    valid_mice = (movie_repeats(:,1) == 10) &(neuropixels_cell_count(:,area,1) >= cell_cutoff);
    current_area = squeeze(neuropixels_population_vectors(valid_mice,area,:)) ;
    across_blocks_ensembles_corr = [];
    diff_maps = [];
    same_maps = [];
    sub = 1;
    for mouse = 1:length(current_area)
        clc;
        disp(['Calculating ensemble rate correlation between blocks:'])
        disp(['Dataset: Neuropixels | Stimulus: Natural movie 1 & Natural movie 3 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area(mouse,:);
        
        movie1_blockA = current_mouse{1}(:,:,1:10);
        movie1_blockB = current_mouse{1}(:,:,11:20);
        
        movie3_blockA = current_mouse{2}(:,:,1:5);
        movie3_blockB = current_mouse{2}(:,:,6:10);
        
        mean_movie1_blockA = nanmean(movie1_blockA,3);
        mean_movie1_blockB = nanmean(movie1_blockB,3);
        
        mean_movie3_blockA = nanmean(movie3_blockA,3);
        mean_movie3_blockB = nanmean(movie3_blockB,3);
        
        
        movie1_blocks_stability = diag(corr(mean_movie1_blockA',mean_movie1_blockB'));
        movie3_blocks_stability = diag(corr(mean_movie3_blockA',mean_movie3_blockB'));
        
        valid_units_both_blocks = (movie1_blocks_stability >= reliability_cutoff) & (movie3_blocks_stability >= reliability_cutoff);
        if sum(valid_units_both_blocks) >= cell_cutoff
            
            valid_movie1_blockA_stability = movie1_blockA(valid_units_both_blocks,:,:);
            valid_movie1_blockB_stability = movie1_blockB(valid_units_both_blocks,:,:);
            valid_movie3_blockA_stability = movie3_blockA(valid_units_both_blocks,:,:);
            valid_movie3_blockB_stability = movie3_blockB(valid_units_both_blocks,:,:);
            
            mean_valid_movie1_blockA_stability = squeeze(nanmean(valid_movie1_blockA_stability,2));
            mean_valid_movie1_blockB_stability = squeeze(nanmean(valid_movie1_blockB_stability,2));
            mean_valid_movie3_blockA_stability = squeeze(nanmean(valid_movie3_blockA_stability,2));
            mean_valid_movie3_blockB_stability = squeeze(nanmean(valid_movie3_blockB_stability,2));
            
            
            
            mean_valid_movie1_half_blockA_stability = [nanmean(mean_valid_movie1_blockA_stability(:,1:5),2),...
                nanmean(mean_valid_movie1_blockA_stability(:,6:10),2)];
            mean_valid_movie1_half_blockB_stability = [nanmean(mean_valid_movie1_blockB_stability(:,1:5),2),...
                nanmean(mean_valid_movie1_blockB_stability(:,6:10),2)];
            
            mean_valid_movie3_half_blockA_stability = [nanmean(mean_valid_movie3_blockA_stability(:,1:2),2),...
                nanmean(mean_valid_movie3_blockA_stability(:,3:5),2)];
            mean_valid_movie3_half_blockB_stability = [nanmean(mean_valid_movie3_blockB_stability(:,1:2),2),...
                nanmean(mean_valid_movie3_blockB_stability(:,3:5),2)];
            
            
            
            across_blocks_ensembles = [mean_valid_movie3_half_blockA_stability,mean_valid_movie1_half_blockA_stability,...
                mean_valid_movie3_half_blockB_stability,mean_valid_movie1_half_blockB_stability];
            
            
            across_blocks_ensembles_corr = corr(across_blocks_ensembles);
            mean_across_blocks_ensembles_corr = [];
            for half1 = 1:4
                rows = [1:2] +2*(half1-1);
                for half2 = 1:4
                    cols = [1:2] +2*(half2-1);
                    current_half = across_blocks_ensembles_corr(rows,cols);
                    mean_across_blocks_ensembles_corr(half1,half2) = nanmean(current_half(:));
                    
                end
            end
            diff_maps(sub,:) = [diag(mean_across_blocks_ensembles_corr,1)',mean_across_blocks_ensembles_corr(1,4)];
            same_maps(sub,:) = [mean_across_blocks_ensembles_corr(1,3),mean_across_blocks_ensembles_corr(2,4)];
            
            sub = sub +1;
            
        end
    end
    diff_maps_areas{area} = diff_maps;
    same_maps_areas{area} = same_maps;
end


same_map_ticks = [20 72];
diff_map_ticks = [0 15 47 77];
ylims=[0.85 1;0.685 1;0.85 1;0.75 1;0.825 1;0.825 1];
plt = [];
figure('units','normalized','position',[0.3 0.3 0.4 0.45])
for area = 1:6
    current_area_diff_maps = diff_maps_areas{area};
    current_area_same_maps = same_maps_areas{area};
    
    mean_diff_maps = nanmean(current_area_diff_maps);
    ste_diff_maps = nanstd(current_area_diff_maps)./sqrt(size(current_area_diff_maps,1));
    
    mean_same_maps = nanmean(current_area_same_maps);
    ste_same_maps = nanstd(current_area_same_maps)./sqrt(size(current_area_same_maps,1));
    
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

cell_cutoff = 15;
reliability_cutoff = 0.5;
area = 1;
valid_mice = neuropixels_cell_count(:,area,1) >= cell_cutoff & movie_repeats(:,1) == 30;
current_area = squeeze(neuropixels_population_vectors(valid_mice,area,:));

across_blocks_ensembles_corr = [];
sub = 1;
for mouse = 1:length(current_area)
    current_mouse = current_area(mouse,:);
    nm1_blockA = current_mouse{1}(:,:,1:30);
    nm1_blockB = current_mouse{1}(:,:,31:60);
    
    snm1_blockA = current_mouse{2}(:,:,1:10);
    snm1_blockB = current_mouse{2}(:,:,11:20);
    
    mean_nm1_blockA = nanmean(nm1_blockA,3);
    mean_nm1_blockB = nanmean(nm1_blockB,3);
    
    mean_snm1_blockA = nanmean(snm1_blockA,3);
    mean_snm1_blockB = nanmean(snm1_blockB,3);
    
    
    nm1_blocks_stability = diag(corr(mean_nm1_blockA',mean_nm1_blockB'));
    snm1_blocks_stability = diag(corr(mean_snm1_blockA',mean_snm1_blockB'));
    
    valid_units_both_blocks = (nm1_blocks_stability >= reliability_cutoff);
    if sum(valid_units_both_blocks) >= cell_cutoff
        valid_nm1_blockA_stability = nm1_blockA(valid_units_both_blocks,:,:);
        valid_nm1_blockB_stability = nm1_blockB(valid_units_both_blocks,:,:);
        valid_snm1_blockA_stability = snm1_blockA(valid_units_both_blocks,:,:);
        valid_snm1_blockB_stability = snm1_blockB(valid_units_both_blocks,:,:);
        
        mean_valid_nm1_blockA_stability = squeeze(nanmean(valid_nm1_blockA_stability,2));
        mean_valid_nm1_blockB_stability = squeeze(nanmean(valid_nm1_blockB_stability,2));
        mean_valid_snm1_blockA_stability = squeeze(nanmean(valid_snm1_blockA_stability,2));
        mean_valid_snm1_blockB_stability = squeeze(nanmean(valid_snm1_blockB_stability,2));
        
        mean_valid_nm1_half_blockA_stability = [nanmean(mean_valid_nm1_blockA_stability(:,1:5),2),...
            nanmean(mean_valid_nm1_blockA_stability(:,6:10),2)];
        mean_valid_nm1_half_blockB_stability = [nanmean(mean_valid_nm1_blockB_stability(:,1:5),2),...
            nanmean(mean_valid_nm1_blockB_stability(:,6:10),2)];
        
        mean_valid_snm1_half_blockA_stability = [nanmean(mean_valid_snm1_blockA_stability(:,1:2),2),...
            nanmean(mean_valid_snm1_blockA_stability(:,3:5),2)];
        mean_valid_snm1_half_blockB_stability = [nanmean(mean_valid_snm1_blockB_stability(:,1:2),2),...
            nanmean(mean_valid_snm1_blockB_stability(:,3:5),2)];
        
        
        
        across_blocks_ensembles = [mean_valid_nm1_half_blockA_stability,mean_valid_snm1_half_blockA_stability...
            mean_valid_snm1_half_blockB_stability,mean_valid_nm1_half_blockB_stability];
        
        
        across_blocks_ensembles_corr(:,:,sub) = corr(across_blocks_ensembles);
        sub = sub + 1;
    end
end



figure
imagesc(nanmean(across_blocks_ensembles_corr,3))
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
cell_cutoff = 15;
area = 1;
valid_mice = neuropixels_cell_count(:,area,1) >= cell_cutoff & movie_repeats(:,1) == 30;
current_area = squeeze(neuropixels_population_vectors(valid_mice,area,:));
reliability_cutoff = 0.5;

movie1_blocks_stability = [];
for mouse = 1:length(current_area)
    current_mouse = current_area(mouse,:);
    movie1_blockA = current_mouse{1}(:,:,1:30);
    movie1_blockB = current_mouse{1}(:,:,31:60);
    mean_movie1_blockA = nanmean(movie1_blockA,3);
    mean_movie1_blockB = nanmean(movie1_blockB,3);
    
    movie1_blocks_stability = [movie1_blocks_stability;diag(corr(mean_movie1_blockA',mean_movie1_blockB'))];
end

figure('units','normalized','position',[0.3 0.3 0.275 0.4])
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
cell_cutoff = 15;
reliability_cutoff = 0.5;

diff_maps_areas = {};
same_maps_areas = {};
for area = 1:6
    valid_mice = neuropixels_cell_count(:,area,1) >= cell_cutoff & movie_repeats(:,1) == 30;
    current_area = squeeze(neuropixels_population_vectors(valid_mice,area,:));
    
    sub = 1;
    diff_maps = [];
    same_maps = [];
    for mouse = 1:length(current_area)
        clc;
        disp(['Calculating ensemble rate correlation between blocks:'])
        disp(['Dataset: Neuropixels | Stimulus: Natural movie 1 & Shuffled natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area(mouse,:);
        nm1_blockA = current_mouse{1}(:,:,1:30);
        nm1_blockB = current_mouse{1}(:,:,31:60);
        
        snm1_blockA = current_mouse{2}(:,:,1:10);
        snm1_blockB = current_mouse{2}(:,:,11:20);
        
        mean_nm1_blockA = nanmean(nm1_blockA,3);
        mean_nm1_blockB = nanmean(nm1_blockB,3);
        
        mean_snm1_blockA = nanmean(snm1_blockA,3);
        mean_snm1_blockB = nanmean(snm1_blockB,3);
        
        
        nm1_blocks_stability = diag(corr(mean_nm1_blockA',mean_nm1_blockB'));
        snm1_blocks_stability = diag(corr(mean_snm1_blockA',mean_snm1_blockB'));
        
        valid_units_both_blocks = (nm1_blocks_stability >= reliability_cutoff);
        if sum(valid_units_both_blocks) >= cell_cutoff
            valid_nm1_blockA_stability = nm1_blockA(valid_units_both_blocks,:,:);
            valid_nm1_blockB_stability = nm1_blockB(valid_units_both_blocks,:,:);
            valid_snm1_blockA_stability = snm1_blockA(valid_units_both_blocks,:,:);
            valid_snm1_blockB_stability = snm1_blockB(valid_units_both_blocks,:,:);
            
            mean_valid_nm1_blockA_stability = squeeze(nanmean(valid_nm1_blockA_stability,2));
            mean_valid_nm1_blockB_stability = squeeze(nanmean(valid_nm1_blockB_stability,2));
            mean_valid_snm1_blockA_stability = squeeze(nanmean(valid_snm1_blockA_stability,2));
            mean_valid_snm1_blockB_stability = squeeze(nanmean(valid_snm1_blockB_stability,2));
            
            mean_valid_nm1_half_blockA_stability = [nanmean(mean_valid_nm1_blockA_stability(:,1:5),2),...
                nanmean(mean_valid_nm1_blockA_stability(:,6:10),2)];
            mean_valid_nm1_half_blockB_stability = [nanmean(mean_valid_nm1_blockB_stability(:,1:5),2),...
                nanmean(mean_valid_nm1_blockB_stability(:,6:10),2)];
            
            mean_valid_snm1_half_blockA_stability = [nanmean(mean_valid_snm1_blockA_stability(:,1:2),2),...
                nanmean(mean_valid_snm1_blockA_stability(:,3:5),2)];
            mean_valid_snm1_half_blockB_stability = [nanmean(mean_valid_snm1_blockB_stability(:,1:2),2),...
                nanmean(mean_valid_snm1_blockB_stability(:,3:5),2)];
            
            
            
            across_blocks_ensembles = [mean_valid_nm1_half_blockA_stability,mean_valid_snm1_half_blockA_stability...
                mean_valid_snm1_half_blockB_stability,mean_valid_nm1_half_blockB_stability];
            
            
            across_blocks_ensembles_corr = corr(across_blocks_ensembles);
            
            mean_across_blocks_ensembles_corr = [];
            for half1 = 1:4
                rows = [1:2] +2*(half1-1);
                for half2 = 1:4
                    cols = [1:2] +2*(half2-1);
                    current_half = across_blocks_ensembles_corr(rows,cols);
                    mean_across_blocks_ensembles_corr(half1,half2) = nanmean(current_half(:));
                    
                end
            end
            diff_maps(sub,:) = [nanmean([mean_across_blocks_ensembles_corr(1,2),mean_across_blocks_ensembles_corr(3,4)]),...
                nanmean([mean_across_blocks_ensembles_corr(1,3),mean_across_blocks_ensembles_corr(2,4)])];
            same_maps(sub,:) = [mean_across_blocks_ensembles_corr(2,3),mean_across_blocks_ensembles_corr(1,4)];
            
            sub = sub + 1;
        end
    end
    diff_maps_areas{area} = diff_maps;
    same_maps_areas{area} = same_maps;
end


same_map_ticks = [60 70];
diff_map_ticks = [0 65];
plt = [];
ylims=[0.69 1;0.7 1;0.65 1;0.69 1;0.7 1;0.775 1];
figure('units','normalized','position',[0.3 0.3 0.4 0.45])
for area = 1:6
    current_area_diff_maps = diff_maps_areas{area};
    current_area_same_maps = same_maps_areas{area};
    
    mean_diff_maps = nanmean(current_area_diff_maps);
    ste_diff_maps = nanstd(current_area_diff_maps)./sqrt(size(current_area_diff_maps,1));
    
    mean_same_maps = nanmean(current_area_same_maps);
    ste_same_maps = nanstd(current_area_same_maps)./sqrt(size(current_area_same_maps,1));
    
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

area = 1;
mouse = 4;
current_mouse = neuropixels_drifting_gratings{mouse,area};
ylims = [0 20; 0 20; 0 3.25];
sub = 1;
plt = [];
figure('units','normalized','position',[0.3 0.3 0.4 0.3])
for cell = [63,104,98]
    current_cell_ori = [];
    for sess =1:3
        current_cell_ori(sess,:) = nanmean(current_mouse(cell,:,sess,:)./2,4);
    end
    
    subplot(1,3,sub)
    hold on
    plt(1) = plot(current_cell_ori(1,:),'color',[0.3 0.3 0.3],'linewidth',2);
    plt(2) = plot(current_cell_ori(2,:),'color',[0.55 0.55 0.55],'linewidth',2);
    plt(3) = plot(current_cell_ori(3,:),'color',[0.8 0.8 0.8],'linewidth',2);
    xlim([1 8])
    ylim(ylims(sub,:))
    title(['Unit #',num2str(cell)])
    if sub == 1
        ylabel('Mean activity rate (Hz)')
    elseif sub == 2
        xlabel('Drifting gratings direction')
    elseif sub == 3
        legend(plt,{'Block A','Block B','Block C'},'Location','northeast')
        legend('boxoff')
    end
    set(gca,'xtick',1:8,'xticklabel',0:45:315)
    xtickangle(45)
    sub = sub + 1;
end

%% Figure S4B - Drifting gratings - PV correlation across V1 mice
area = 1;
cell_cutoff = 20;
valid_mice = calcium_excitatory_cell_count{area}(:,1)>=cell_cutoff;
current_area = calcium_excitatory_drifting_gratings{area}(valid_mice);
pv_corr_across_mice = [];
for mouse = 1:length(current_area)
    current_mouse = current_area{mouse};
    
    current_mouse_sorted = [];
    for sess =1:3
        for ori =1:8
            for freq = 1:5
                current_mouse_sorted = [current_mouse_sorted,current_mouse(:,ori,sess,freq)];
            end
        end
    end
    pv_corr_across_mice(:,:,mouse) = corr(current_mouse_sorted,'rows','pairwise');
end

mean_pv_corr_across_mice = nanmean(pv_corr_across_mice,3);
mean_pv_corr_across_mice(boolean(eye(size(mean_pv_corr_across_mice,1)))) = NaN;
mean_pv_corr_across_mice(isnan(mean_pv_corr_across_mice)) = max(mean_pv_corr_across_mice(:));

figure
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

pv_ori = [];
for row = 1:16
    rows_ind = [41:45]+5*(row-1);
    pv_ori(:,:,row) = mean_pv_corr_across_mice(rows_ind-40,rows_ind);
    
    if row<9
        pv_ori(:,:,row+16) = mean_pv_corr_across_mice(rows_ind-40,rows_ind+40);
    end
    
end

subplot(1,2,2,'units','normalized','position',[0.6 0.7 0.15 0.2])
imagesc(nanmean(pv_ori,3))
set(gca,'xtick',1:5,'xticklabel',[1,2,4,8,15],'ytick',[])
title({'Temporal';'frequancy (Hz)'})
colormap(newmap3)

%% Figure S4C - Drifting gratings - PV correlation across directions
cell_cutoff = 20;

for area = 1:6
    valid_mice = calcium_excitatory_cell_count{area}(:,1)>=cell_cutoff;
    current_area = calcium_excitatory_drifting_gratings{area}(valid_mice);
    mean_elapsed_pv_shuffled = [];
    mean_elapsed_pv = [];
    for mouse = 1:length(current_area)
        clc;
        disp(['Calculating population vector correlation between blocks:'])
        disp(['Dataset: Calcium imaging | Stimulus: Drifting gratings | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse};
        
        current_mouse_sorted = [];
        
        
        for ori =1:size(current_mouse,2)
            for freq = 1:size(current_mouse,4)
                current_mouse_sorted = [current_mouse_sorted,current_mouse(:,ori,:,freq)];
            end
        end
        
        current_mouse_sorted_shuffled = current_mouse_sorted;
        for cell = 1:size(current_mouse_sorted,1)
            for sess = 1:3
                current_mouse_sorted_shuffled(cell,:,sess) = current_mouse_sorted(cell,circshift([1:40],randperm(40,1)),sess);
            end
        end
        
        
        elapsed_pv = [];
        elapsed_pv_shuffled = [];
        ori_list = [5:5:35];
        ori_list2 = [-35:5:-5];
        block = 1;
        for block1 = 1:size(current_mouse_sorted,3)
            current_mouse_block1 = current_mouse_sorted(:,:,block1);
            current_mouse_block1_shuffled = current_mouse_sorted_shuffled(:,:,block1);
            
            for block2 = 1:size(current_mouse_sorted,3)
                current_mouse_block2 = current_mouse_sorted(:,:,block2);
                current_mouse_block2_shuffled = current_mouse_sorted_shuffled(:,:,block2);
                
                structure = corr(current_mouse_block1,current_mouse_block2);
                structure_shuffled = corr(current_mouse_block1_shuffled,current_mouse_block2_shuffled);
                
                if block1 < block2
                    for diagonal = 1:length(ori_list)+1
                        if diagonal == 1
                            elapsed_pv(block,diagonal) = nanmean(diag(structure));
                            elapsed_pv_shuffled(block,diagonal) = nanmean(diag(structure_shuffled));
                            
                        else
                            elapsed_pv(block,diagonal) = nanmean([diag(structure,ori_list(diagonal-1));diag(structure,ori_list2(diagonal-1))]);
                            elapsed_pv_shuffled(block,diagonal) = nanmean([diag(structure_shuffled,ori_list(diagonal-1));diag(structure_shuffled,ori_list2(diagonal-1))]);
                            
                        end
                    end
                    block = block + 1;
                end
                
            end
        end
        
        mean_elapsed_pv(:,:,mouse) = [nanmean(elapsed_pv([1,3],:,:));elapsed_pv(2,:,:)];
        mean_elapsed_pv_shuffled(:,:,mouse) = [nanmean(elapsed_pv_shuffled([1,3],:,:));elapsed_pv_shuffled(2,:,:)];
    end
    
    mean_elapsed_pv_area{area} = mean_elapsed_pv;
    mean_elapsed_pv_area_shuffled{area} = mean_elapsed_pv_shuffled;
end


figure('units','normalized','position',[0.3 0.3 0.3 0.3])
for area = 1:6
    subplot(2,3,area)
    
    current_area = nanmean([mean_elapsed_pv_area{area}(:,5:8,:),mean_elapsed_pv_area{area}(:,1:5,:)]);
    mean_current_area_shuffled = nanmean([mean_elapsed_pv_area_shuffled{area}(:,5:8,:),mean_elapsed_pv_area_shuffled{area}(:,1:5,:)]);
    
    mean_pv = nanmean(current_area,3);
    ste_pv = nanstd(current_area,[],3)./sqrt(size(current_area,3));
    
    mean_pv_shuffled = nanmean(mean_current_area_shuffled,3);
    ste_pv_shuffled = nanstd(mean_current_area_shuffled,[],3)./sqrt(size(mean_current_area_shuffled,3));
    
    hold on
    errorbar(mean_pv_shuffled(1,:),ste_pv_shuffled,'o','color',[0.7 0.7 0.7],...
        'markerfacecolor',[0.7 0.7 0.7],'capsize',0,'linestyle','-','linewidth',2,'markersize',2)
    
    errorbar(mean_pv(1,:),ste_pv,'o','color',colors(area,:),...
        'markerfacecolor',colors(area,:),'capsize',0,'linestyle','-','linewidth',2,'markersize',2)
    
    text(0.1,0.95,brain_areas{area},'Units','normalized','FontSize',12)
    text(0.05,0.85,['N=',num2str(sum(size(current_area,3),2))],'Units','normalized','FontSize',10)
    
    text(0.8,0.95,['Real'],'Units','normalized','FontSize',10,'Color',colors(area,:))
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
cell_cutoff = 20;
for area = 1:6
    valid_mice = calcium_excitatory_cell_count{area}(:,1)>=cell_cutoff;
    current_area = calcium_excitatory_drifting_gratings{area}(valid_mice);
    
    elapsed_block_rate_corr = [];
    for mouse = 1:length(current_area)
        clc;
        disp(['Calculating ensemble rate correlation between blocks:'])
        disp(['Dataset: Calcium imaging | Stimulus: Drifting gratings | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse};
        
        current_mouse_sorted = [];
        for freq = 1:size(current_mouse,4)
            current_mouse_sorted = [current_mouse_sorted,current_mouse(:,:,:,freq)];
        end
        
        rate_corr = corr(squeeze(nanmean(current_mouse_sorted,2)));
        
        elapsed_block_rate_corr(mouse,:) = [nanmean(diag(rate_corr,1)),diag(rate_corr,2)];
        
    end
    
    elapsed_block_rate_corr_area{area} = elapsed_block_rate_corr;
    
end

pvalues = [];
zvalues = [];
plt = [];
figure('units','normalized','position',[0.35 0.35 0.2 0.3])
for area = 1:6
    current_area = elapsed_block_rate_corr_area{area};
    
    [pvalues(area),~,stats] = signrank(current_area(:,1),current_area(:,2));
    zvalues(area) = stats.zval;
    
    mean_stability = nanmean(current_area);
    std_stability = nanstd(current_area)./sqrt(size(current_area,1));
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


corrected_pval = bonf_holm(pvalues);
VarNames = {'area','zvalue','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas(1:6)'),zvalues(:),pvalues(:),corrected_pval(:),'VariableNames',VarNames);

clc;
disp(['Ensemble rate correlation between proximal blocks compared to between distal blocks'])
disp(['Two-tailed Wilcoxon signed-rank tests with Holm–Bonferroni correction:'])
disp(statistics)

%% Figure S4E - Drifting gratings - Tuning curve correlation across blocks (calcium imaging)
cell_cutoff = 20;
for area = 1:6
    valid_mice = calcium_excitatory_cell_count{area}(:,1)>=cell_cutoff;
    current_area = calcium_excitatory_drifting_gratings{area}(valid_mice);
    
    elapsed_block_tuning_corr = [];
    for mouse = 1:length(current_area)
        clc;
        disp(['Calculating tuning curve correlation between blocks:'])
        disp(['Dataset: Calcium imaging | Stimulus: Drifting gratings | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse};
        
        current_mouse_sorted = [];
        for freq = 1:size(current_mouse,4)
            current_mouse_sorted = [current_mouse_sorted,current_mouse(:,:,:,freq)];
        end
        
        mean_tuning_corr = [];
        for block1 = 1:size(current_mouse_sorted,3)
            current_mouse_block1 = current_mouse_sorted(:,:,block1);
            for block2 = 1:size(current_mouse_sorted,3)
                current_mouse_block2 = current_mouse_sorted(:,:,block2);
                mean_tuning_corr(block1,block2) = nanmedian(diag(corr(current_mouse_block1',current_mouse_block2','rows','complete')));
            end
        end
        elapsed_block_tuning_corr(mouse,:) = [nanmean(diag(mean_tuning_corr,1)),diag(mean_tuning_corr,2)];
    end
    elapsed_block_tuning_corr_area{area} = elapsed_block_tuning_corr;
end

pvalues = [];
zvalues = [];
figure('units','normalized','position',[0.35 0.35 0.2 0.3])
for area = 1:6
    current_area = elapsed_block_tuning_corr_area{area};
    
    [pvalues(area),~,stats] = signrank(current_area(:,1),current_area(:,2));
    zvalues(area) = stats.zval;
    
    mean_stability = nanmean(current_area);
    ste_stability = nanstd(current_area)./sqrt(size(current_area,1));
    hold on
    errorbar(mean_stability,ste_stability,'o','color',colors(area,:),...
        'markerfacecolor',colors(area,:),'capsize',0,'linestyle','-','linewidth',3,'markersize',3);
    
end

ylabel('Tuning curve correlation')
set(gca,'xtick',1:2,'xticklabels',[15,30])
xlabel('Elapsed time (min)')
xlim([0.5 2.5])


corrected_pval = bonf_holm(pvalues);
VarNames = {'area','zvalue','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas(1:6)'),zvalues(:),pvalues(:),corrected_pval(:),'VariableNames',VarNames);

clc;
disp(['Tuning curve correlation between proximal blocks compared to between distal blocks'])
disp(['Two-tailed Wilcoxon signed-rank tests with Holm–Bonferroni correction:'])
disp(statistics)


%% Figure S4F - Drifting gratings - Ensemble rate correlation across blocks (neuropixels)
cell_cutoff = 15;

for area = 1:6
    valid_mice = neuropixels_cell_count(:,area,1) >= cell_cutoff & movie_repeats(:,1) == 10;
    current_area = neuropixels_drifting_gratings(valid_mice,area);
    elapsed_block_rate_corr = [];
    
    for mouse =1:length(current_area)
        clc;
        disp(['Calculating ensemble rate correlation between blocks:'])
        disp(['Dataset: Neuropixels | Stimulus: Drifting gratings | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse};
        current_mouse_sorted = [];
        for freq = 1:5
            current_mouse_sorted = [current_mouse_sorted,current_mouse(:,:,:,freq)];
        end
        
        mean_rate_corr = corr(squeeze(nanmean(current_mouse_sorted,2)));
        elapsed_block_rate_corr(mouse,:) = [nanmean(diag(mean_rate_corr,1)),diag(mean_rate_corr,2)];
    end
    elapsed_block_rate_corr_area{area} = elapsed_block_rate_corr;
end


pvalues = [];
zvalues = [];
plt = [];
figure('units','normalized','position',[0.35 0.35 0.2 0.3])
for area = 1:6
    current_area = elapsed_block_rate_corr_area{area};
    
    [pvalues(area),~,stats] = signrank(current_area(:,1),current_area(:,2));
    zvalues(area) = stats.zval;
    
    mean_stability = nanmean(current_area);
    ste_stability = nanstd(current_area)./sqrt(size(current_area,1));
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


corrected_pval = bonf_holm(pvalues);
VarNames = {'area','zvalue','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas(1:6)'),zvalues(:),pvalues(:),corrected_pval(:),'VariableNames',VarNames);

clc;
disp(['Ensemble rate correlation between proximal blocks compared to between distal blocks'])
disp(['Two-tailed Wilcoxon signed-rank tests with Holm–Bonferroni correction:'])
disp(statistics)

%% Figure S4G - Drifting gratings - Tuning curve correlation across blocks (neuropixels)
cell_cutoff = 15;

for area = 1:6
    valid_mice = neuropixels_cell_count(:,area,1) >= cell_cutoff & movie_repeats(:,1) == 10;
    current_area = neuropixels_drifting_gratings(valid_mice,area);
    
    elapsed_block_tuning_corr = [];
    for mouse =1:length(current_area)
        clc;
        disp(['Calculating tuning curve correlation between blocks:'])
        disp(['Dataset: Neuropixels | Stimulus: Drifting gratings | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse};
        
        current_mouse_sorted = [];
        for freq = 1:5
            current_mouse_sorted = [current_mouse_sorted,current_mouse(:,:,:,freq)];
        end
        
        mean_tuning_corr = [];
        for block1 = 1:3
            for block2 = 1:3
                mean_tuning_corr(block1,block2) = nanmedian(diag(corr(current_mouse_sorted(:,:,block1)',current_mouse_sorted(:,:,block2)','rows','pairwise')));
            end
        end
        
        elapsed_block_tuning_corr(mouse,:) = [nanmean(diag(mean_tuning_corr,1)),diag(mean_tuning_corr,2)];
    end
    elapsed_block_tuning_corr_area{area} = elapsed_block_tuning_corr;
end

pvalues = [];
zvalues = [];
figure('units','normalized','position',[0.35 0.35 0.2 0.3])
for area = 1:6
    current_area = elapsed_block_tuning_corr_area{area};
    
    [pvalues(area),~,stats] = signrank(current_area(:,1),current_area(:,2));
    zvalues(area) = stats.zval;
    
    mean_stability = nanmean(current_area);
    ste_stability = nanstd(current_area)./sqrt(size(current_area,1));
    hold on
    errorbar(mean_stability,ste_stability,'o','color',colors(area,:),...
        'markerfacecolor',colors(area,:),'capsize',0,'linestyle','-','linewidth',3,'markersize',3);
    
end

ylabel('Tuning curve correlation')
set(gca,'xtick',1:2,'xticklabels',[15,30])
xlabel('Elapsed time (min)')
xlim([0.5 2.5])


corrected_pval = bonf_holm(pvalues);
VarNames = {'area','zvalue','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas(1:6)'),zvalues(:),pvalues(:),corrected_pval(:),'VariableNames',VarNames);

clc;
disp(['Tuning curve correlation between proximal blocks compared to between distal blocks'])
disp(['Two-tailed Wilcoxon signed-rank tests with Holm–Bonferroni correction:'])
disp(statistics)

%% Figure S5A - Ensemble rate and tuning curve correlation difference between proximal and distal sessions
rate_tuning_corr_diff_across_areas = {};
nat_movie = 1;
cell_cutoff = 20;
for area = 1:6
    valid_mice = min(calcium_excitatory_cell_count{area},[],2)>= cell_cutoff;
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie);
    
    elapsed_rate_corr = [];
    elapsed_tuning_corr = [];
    for mouse = 1:length(current_area)
        clc;
        disp(['Calculating ensemble rate and tuning correlation between ssessions:'])
        disp(['Dataset: Calcium imaging | Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse};
        
        mean_activity_per_half = [];
        for half = 1:3
            current_half = [1:10] + 10*(half-1);
            mean_activity_per_half(:,:,half) = nanmean(current_mouse(:,:,current_half),3);
        end
        
        rate_corr = [];
        tuning_corr = [];
        for halfA = 1:size(mean_activity_per_half,3)
            halfA_activity = mean_activity_per_half(:,:,halfA);
            for halfB = 1:size(mean_activity_per_half,3)
                halfB_activity = mean_activity_per_half(:,:,halfB);
                
                valid_cells = [nanmean(halfA_activity,2) ~= 0] & [nanmean(halfB_activity,2) ~= 0];
                valid_halfA_activity = halfA_activity(valid_cells,:);
                valid_halfB_activity = halfB_activity(valid_cells,:);
                
                rate_corr(halfA,halfB) = corr(nanmean(valid_halfA_activity,2),nanmean(valid_halfB_activity,2));
                tuning_corr(halfA,halfB) = nanmedian(diag(corr(valid_halfA_activity',valid_halfB_activity')));
            end
        end
        
        elapsed_rate_corr(mouse,:) = [nanmean(diag(rate_corr,1)),diag(rate_corr,2)];
        elapsed_tuning_corr(mouse,:) = [nanmean(diag(tuning_corr,1)),diag(tuning_corr,2)];
    end
    
    rate_corr_diff = elapsed_rate_corr(:,1) - elapsed_rate_corr(:,2);
    tuning_corr_diff = elapsed_tuning_corr(:,1) - elapsed_tuning_corr(:,2);
    rate_tuning_corr_diff_across_areas{area} = [rate_corr_diff,tuning_corr_diff];
    
end


pvalue = [];
zvalue = [];
mean_stability = [];
ste_stability = [];

for area = 1:6
    current_area = rate_tuning_corr_diff_across_areas{area};
    mean_stability(area,:) = nanmean(current_area);
    std_stability = nanstd(current_area);
    ste_stability(area,:) = std_stability./sqrt(size(current_area,1));
    
    [pvalue(area,1),~,stats] = signrank(current_area(:,1));
    zvalue(area,1) = stats.zval;
    
    [pvalue(area,2),~,stats] = signrank(current_area(:,2));
    zvalue(area,2) = stats.zval;
end

plt = [];
figure('Units','Normalized','Position',[0.3 0.4 0.225 0.325])
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

pvalue = pvalue./2; % one tail
corrected_pvalue = [];
corrected_pvalue(:,1) = bonf_holm(pvalue(:,1));
corrected_pvalue(:,2) = bonf_holm(pvalue(:,2));

VarNames = {'area','Rate_zvalue','Rate_pvalue','Rate_bonf_holm','Tuning_zvalue','Tuning_pvalue','Tuning_bonf_holm'};
statistics = table(cell2mat(brain_areas(1:6)'),zvalue(:,1),pvalue(:,1),corrected_pvalue(:,1),...
    zvalue(:,2),pvalue(:,2),corrected_pvalue(:,2),'VariableNames',VarNames);

clc;
disp(['Ensemble rate and tuning curve correlation within blocks compared to between blocks'])
disp(['One-tailed Wilcoxon signed-rank tests with Holm–Bonferroni correction:'])
disp(statistics)

%% Figure S5B - Mean activity rate across sessions
cell_cutoff = 20;
mean_activity_all_mice_area = {};
for area = 1:6
    valid_mice = min(calcium_excitatory_cell_count{area},[],2)>= cell_cutoff;
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie);
    
    mean_activity_all_mice = [];
    for mouse = 1:length(current_area )
        current_mouse_activity = current_area{mouse}*30;
        sess1 = current_mouse_activity(:,:,1:10);
        valid_sess1 = nanmean(nanmean(sess1,2),3)>0;
        sess2 = current_mouse_activity(:,:,11:20);
        valid_sess2 = nanmean(nanmean(sess2,2),3)>0;
        sess3 = current_mouse_activity(:,:,21:30);
        valid_sess3 = nanmean(nanmean(sess3,2),3)>0;
        
        mean_sess1 = nanmean(squeeze(nanmean(sess1(valid_sess1,:,:),2)),2);
        mean_sess2 = nanmean(squeeze(nanmean(sess2(valid_sess2,:,:),2)),2);
        mean_sess3 = nanmean(squeeze(nanmean(sess3(valid_sess3,:,:),2)),2);
        
        mean_activity_all_mice(mouse,:) = [nanmean(mean_sess1),nanmean(mean_sess2),nanmean(mean_sess3)];
    end
    mean_activity_all_mice_area{area} = mean_activity_all_mice;
end

pvalues = [];
zvalues = [];
figure('units','normalized','position',[0.3 0.3 0.315 0.35])
for area = 1:6
    current_area = mean_activity_all_mice_area{area};
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
    
    sig_mat = [];
    for day1 = 1:3
        for day2 = 1:3
            if day1 ~= day2
                [sig_mat(day1,day2),~,stat] =signrank(current_area(:,day1),current_area(:,day2));
                z_mat(day1,day2) = stat.zval;
            end
        end
    end
    pvalues(area,:) = [sig_mat(1,2),sig_mat(2,3),sig_mat(1,3)];
    zvalues(area,:) = [z_mat(1,2),z_mat(2,3),z_mat(1,3)];
end

VarNames = {'area','zvalue_sess1_sess2','pvalue_sess1_sess2','zvalue_sess2_sess3','pvalue_sess2_sess3','zvalue_sess1_sess3','pvalue_sess1_sess3'};
statistics = table(cell2mat(brain_areas(1:6)'),zvalues(:,1),pvalues(:,1),zvalues(:,2),pvalues(:,2),zvalues(:,3),pvalues(:,3),'VariableNames',VarNames);

clc;
disp(['Difference in activity rates across sessions'])
disp(['Two-tailed Wilcoxon signed-rank tests with Holm–Bonferroni correction:'])
disp(statistics)

%% Figure S5C - Number of active cells across sessions
cell_cutoff = 20;
mean_num_active_cells_area = {};
for area = 1:6
    valid_mice = min(calcium_excitatory_cell_count{area},[],2)>= cell_cutoff;
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie);
    
    mean_num_active_cells = [];
    for mouse = 1:length(current_area )
        current_mouse_activity = current_area{mouse}*30;
        sess1 = current_mouse_activity(:,:,1:10);
        valid_sess1 = nanmean(nanmean(sess1,2),3)>0;
        sess2 = current_mouse_activity(:,:,11:20);
        valid_sess2 = nanmean(nanmean(sess2,2),3)>0;
        sess3 = current_mouse_activity(:,:,21:30);
        valid_sess3 = nanmean(nanmean(sess3,2),3)>0;
        
        mean_num_active_cells(mouse,:) = [sum(valid_sess1),sum(valid_sess2),sum(valid_sess3)];
    end
    mean_num_active_cells_area{area} = mean_num_active_cells;
end

pvalues = [];
zvalues = [];
figure('units','normalized','position',[0.3 0.3 0.315 0.35])
for area = 1:6
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
    
    sig_mat = [];
    for day1 = 1:3
        for day2 = 1:3
            if day1 ~= day2
                [sig_mat(day1,day2),~,stat] =signrank(current_area(:,day1),current_area(:,day2));
                z_mat(day1,day2) = stat.zval;
            end
        end
    end
    pvalues(area,:) = [sig_mat(1,2),sig_mat(2,3),sig_mat(1,3)];
    zvalues(area,:) = [z_mat(1,2),z_mat(2,3),z_mat(1,3)];
end

VarNames = {'area','zvalue_sess1_sess2','pvalue_sess1_sess2','zvalue_sess2_sess3','pvalue_sess2_sess3','zvalue_sess1_sess3','pvalue_sess1_sess3'};
statistics = table(cell2mat(brain_areas(1:6)'),zvalues(:,1),pvalues(:,1),zvalues(:,2),pvalues(:,2),zvalues(:,3),pvalues(:,3),'VariableNames',VarNames);

clc;
disp(['Difference in number of active cells across sessions'])
disp(['Two-tailed Wilcoxon signed-rank tests with Holm–Bonferroni correction:'])
disp(statistics)


%% Figure S5D - Mean running speed across sessions
cell_cutoff = 20;
mean_running_speed_area = {};
for area = 1:6
    valid_mice = min(calcium_excitatory_cell_count{area},[],2)>= cell_cutoff;
    current_area = calcium_excitatory_running_speed{area}(valid_mice,nat_movie);
    
    mean_running_speed = [];
    for mouse = 1:length(current_area )
        current_mouse = current_area{mouse};
        sess1 = nanmean(nanmean(current_mouse(:,:,1:10),2),3);
        sess2 = nanmean(nanmean(current_mouse(:,:,11:20),2),3);
        sess3 = nanmean(nanmean(current_mouse(:,:,21:30),2),3);
        
        mean_running_speed(mouse,:) = [sess1,sess2,sess3];
    end
    mean_running_speed_area{area,1} = mean_running_speed;
end

mean_running_speed_pooled = cell2mat(mean_running_speed_area);

figure('units','normalized','position',[0.35 0.35 0.2 0.3])
figure_boxplot(mean_running_speed_pooled);
xlim([0.1 3.9])
ylabel({'Mean running speed (cm/sec)'})
xlabel('Session')

sig_mat = [];
for day1 = 1:3
    for day2 = 1:3
        if day1 ~= day2
            [sig_mat(day1,day2),~,stat] =signrank(mean_running_speed_pooled(:,day1),mean_running_speed_pooled(:,day2));
            z_mat(day1,day2) = stat.zval;
        end
    end
end
pvalues = [sig_mat(1,2),sig_mat(2,3),sig_mat(1,3)];
zvalues = [z_mat(1,2),z_mat(2,3),z_mat(1,3)];
stats = [pvalues;zvalues];

VarNames = {'statistic','sess1_sess2','sess2_sess3','sess1_sess3'};
statistics = table({'pvalue';'zvalue'},stats(:,1),stats(:,2),stats(:,3),'VariableNames',VarNames);

clc;
disp(['Difference in mean running speed across sessions'])
disp(['Two-tailed Wilcoxon signed-rank tests with Holm–Bonferroni correction:'])
disp(statistics)

%% Figure S5E - Mean pupil size across sessions
cell_cutoff = 20;
mean_pupil_size_area = {};
for area = 1:6
    valid_mice = min(calcium_excitatory_cell_count{area},[],2)>= cell_cutoff;
    current_area = calcium_excitatory_pupil_size{area}(valid_mice,nat_movie);
    
    mean_pupil_size = [];
    for mouse = 1:length(current_area )
        current_mouse = current_area{mouse};
        sess1 = nanmean(nanmean(current_mouse(:,:,1:10),2),3);
        sess2 = nanmean(nanmean(current_mouse(:,:,11:20),2),3);
        sess3 = nanmean(nanmean(current_mouse(:,:,21:30),2),3);
        
        mean_pupil_size(mouse,:) = [sess1,sess2,sess3];
    end
    mean_pupil_size_area{area,1} = mean_pupil_size;
end

mean_pupil_size_pooled = cell2mat(mean_pupil_size_area);

figure('units','normalized','position',[0.35 0.35 0.2 0.3])
figure_boxplot(mean_pupil_size_pooled);
xlim([0.1 3.9])
ylabel({'Pupil size (a.u.)'})
xlabel('Session')

sig_mat = [];
for day1 = 1:3
    for day2 = 1:3
        if day1 ~= day2
            [sig_mat(day1,day2),~,stat] =signrank(mean_pupil_size_pooled(:,day1),mean_pupil_size_pooled(:,day2));
            z_mat(day1,day2) = stat.zval;
        end
    end
end
pvalues = [sig_mat(1,2),sig_mat(2,3),sig_mat(1,3)];
zvalues = [z_mat(1,2),z_mat(2,3),z_mat(1,3)];
stats = [pvalues;zvalues];

VarNames = {'statistic','sess1_sess2','sess2_sess3','sess1_sess3'};
statistics = table({'pvalue';'zvalue'},stats(:,1),stats(:,2),stats(:,3),'VariableNames',VarNames);

clc;
disp(['Difference in mean pupil size across sessions'])
disp(['Two-tailed Wilcoxon signed-rank tests with Holm–Bonferroni correction:'])
disp(statistics)

%% Figure S5F - Within day decoder across sessions
within_day_decoder_areas = {};
nat_movie = 1;
cell_cutoff = 20;
for area = 1:6
    valid_mice = min(calcium_excitatory_cell_count{area},[],2)>= cell_cutoff;
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie);
    within_day_decoder = [];
    for mouse = 1:length(current_area)
        clc;
        disp(['Performing time-lapse decoding within sessions:'])
        disp(['Dataset: Calcium imaging | Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse};
        for day = 1:3
            current_day = current_mouse(:,:,[1:10]+10*(day-1));
            valid_cells = nanmean(nanmean(current_day,3),2)>0;
            
            mdl = fitcknn(nanmean(current_day(valid_cells ,:,1:5),3)',[1:30]);
            ypred = predict(mdl,nanmean(current_day(valid_cells ,:,6:10),3)');
            
            within_day_decoder(mouse,day) = nanmean([(sum([ypred] == [1:30]')./30)*100]);
        end
    end
    within_day_decoder_areas{area} = within_day_decoder;
end

pvalues = [];
zvalues = [];
figure('units','normalized','position',[0.3 0.3 0.315 0.35])
for area = 1:6
    current_area = within_day_decoder_areas{area};
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
    
    sig_mat = [];
    for day1 = 1:3
        for day2 = 1:3
            if day1 ~= day2
                [sig_mat(day1,day2),~,stat] =signrank(current_area(:,day1),current_area(:,day2));
                z_mat(day1,day2) = stat.zval;
            end
        end
    end
    pvalues(area,:) = [sig_mat(1,2),sig_mat(2,3),sig_mat(1,3)];
    zvalues(area,:) = [z_mat(1,2),z_mat(2,3),z_mat(1,3)];
end

VarNames = {'area','zvalue_sess1_sess2','pvalue_sess1_sess2','zvalue_sess2_sess3','pvalue_sess2_sess3','zvalue_sess1_sess3','pvalue_sess1_sess3'};
statistics = table(cell2mat(brain_areas(1:6)'),zvalues(:,1),pvalues(:,1),zvalues(:,2),pvalues(:,2),zvalues(:,3),pvalues(:,3),'VariableNames',VarNames);

clc;
disp(['Difference in within day decoder performance across sessions'])
disp(['Two-tailed Wilcoxon signed-rank tests with Holm–Bonferroni correction:'])
disp(statistics)

%% Figure S5G - With day PV correlation across sessions
within_day_pv_areas = {};
nat_movie = 1;
cell_cutoff = 20;
for area = 1:6
    valid_mice = min(calcium_excitatory_cell_count{area},[],2)>= cell_cutoff;
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie);
    within_day_pv = [];
    for mouse = 1:length(current_area)
        clc;
        disp(['Calculation population vector correlation within sessions:'])
        disp(['Dataset: Calcium imaging | Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse};
        for day = 1:3
            current_day = current_mouse(:,:,[1:10]+10*(day-1));
            valid_cells = nanmean(nanmean(current_day,3),2)>0;
            
            within_day_pv(mouse,day) = nanmean(diag(corr(nanmean(current_day(valid_cells ,:,1:5),3),nanmean(current_day(valid_cells ,:,6:10),3))));
        end
    end
    within_day_pv_areas{area} = within_day_pv;
end

pvalues = [];
zvalues = [];
figure('units','normalized','position',[0.3 0.3 0.315 0.35])
for area = 1:6
    current_area = within_day_pv_areas{area};
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
    
    sig_mat = [];
    for day1 = 1:3
        for day2 = 1:3
            if day1 ~= day2
                [sig_mat(day1,day2),~,stat] =signrank(current_area(:,day1),current_area(:,day2));
                z_mat(day1,day2) = stat.zval;
            end
        end
    end
    pvalues(area,:) = [sig_mat(1,2),sig_mat(2,3),sig_mat(1,3)];
    zvalues(area,:) = [z_mat(1,2),z_mat(2,3),z_mat(1,3)];
end

VarNames = {'area','zvalue_sess1_sess2','pvalue_sess1_sess2','zvalue_sess2_sess3','pvalue_sess2_sess3','zvalue_sess1_sess3','pvalue_sess1_sess3'};
statistics = table(cell2mat(brain_areas(1:6)'),zvalues(:,1),pvalues(:,1),zvalues(:,2),pvalues(:,2),zvalues(:,3),pvalues(:,3),'VariableNames',VarNames);

clc;
disp(['Difference in within day PV correlation across sessions'])
disp(['Two-tailed Wilcoxon signed-rank tests with Holm–Bonferroni correction:'])
disp(statistics)

%% Figure S5H - PV correlation difference between pairs of subsequent sessions
pv_corr_diff_across_areas = nan(100,6);
nat_movie = 1;
cell_cutoff = 20;
pvalue = [];
zvalue = [];
for area = 1:6
    valid_mice = min(calcium_excitatory_cell_count{area},[],2)>= cell_cutoff;
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie);
    
    elapsed_pv_corr = [];
    for mouse = 1:length(current_area)
        clc;
        disp(['Calculating population vector correlation between ssessions:'])
        disp(['Dataset: Calcium imaging | Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse = current_area{mouse};
        
        mean_activity_per_half = [];
        for half = 1:3
            current_half = [1:10] + 10*(half-1);
            mean_activity_per_half(:,:,half) = nanmean(current_mouse(:,:,current_half),3);
        end
        
        pv_corr = [];
        for halfA = 1:size(mean_activity_per_half,3)
            halfA_activity = mean_activity_per_half(:,:,halfA);
            for halfB = 1:size(mean_activity_per_half,3)
                halfB_activity = mean_activity_per_half(:,:,halfB);
                
                valid_cells = [nanmean(halfA_activity,2) ~= 0] & [nanmean(halfB_activity,2) ~= 0];
                valid_halfA_activity = halfA_activity(valid_cells,:);
                valid_halfB_activity = halfB_activity(valid_cells,:);
                
                pv_corr(halfA,halfB) = nanmean(diag(corr(valid_halfA_activity,valid_halfB_activity)));
            end
        end
        
        elapsed_pv_corr(mouse,:) = [diag(pv_corr,1)];
    end
    
    pv_corr_diff = elapsed_pv_corr(:,1) - elapsed_pv_corr(:,2);
    pv_corr_diff_across_areas(1:mouse,area) = pv_corr_diff;
    [pvalue(area),~,stat] = signrank(elapsed_pv_corr(:,1),elapsed_pv_corr(:,2));
    zvalue(area) = stat.zval;
end

figure('units','normalized','position',[0.3 0.3 0.2 0.3])
xlim([0 7])
ylim([-0.6 0.6])
figure_boxplot(pv_corr_diff_across_areas);
ylabel('PV correlation difference')
set(gca,'xtick',1:6,'xticklabel',brain_areas(1:6))

VarNames = {'area','zvalue','pvalue'};
statistics = table(cell2mat(brain_areas(1:6)'),zvalue',pvalue','VariableNames',VarNames);

clc;
disp(['Difference in the PV correlation values between pairs of subsequent sessions'])
disp(['Two-tailed Wilcoxon signed-rank tests with Holm–Bonferroni correction:'])
disp(statistics)

%% Figure S5I - PV correlation across sessions with cells turnover
pv_corr_areas = {};
nat_movie = 1;
cell_cutoff = 20;
for area = 1:6
    valid_mice = min(calcium_excitatory_cell_count{area},[],2)>= cell_cutoff;
    current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie);
    for subset = 1:2
        within_between_session_stability = [];
        for mouse = 1:length(current_area)
            clc;
            disp(['Calculating population vector correlation between ssessions:'])
            disp(['Dataset: Calcium imaging | Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
            
            current_mouse = current_area{mouse};
            
            mean_activity_per_half = [];
            for half = 1:6
                current_half = [1:5] + 5*(half-1);
                mean_activity_per_half(:,:,half) = nanmean(current_mouse(:,:,current_half),3);
            end
            
            pv_corr = [];
            for halfA = 1:size(mean_activity_per_half,3)
                halfA_activity = mean_activity_per_half(:,:,halfA);
                for halfB = 1:size(mean_activity_per_half,3)
                    halfB_activity = mean_activity_per_half(:,:,halfB);
                    
                    if subset == 1
                        valid_cells = [nanmean(halfA_activity,2) ~= 0] & [nanmean(halfB_activity,2) ~= 0];
                    elseif subset == 2
                        valid_cells = [nanmean(halfA_activity,2) ~= 0] | [nanmean(halfB_activity,2) ~= 0];
                    end
                    
                    valid_halfA_activity = halfA_activity(valid_cells,:);
                    valid_halfB_activity = halfB_activity(valid_cells,:);
                    
                    pv_corr(halfA,halfB) = nanmean(diag(corr(valid_halfA_activity,valid_halfB_activity)));
                end
            end
            
            within_between_session_stability(mouse,1) = nanmean([pv_corr(1,2),pv_corr(3,4),pv_corr(5,6)]);
            within_between_session_stability(mouse,2) = nanmean([nanmean(nanmean(pv_corr(1:2,3:4))),...
                nanmean(nanmean(pv_corr(3:4,5:6)))]);
            within_between_session_stability(mouse,3) = nanmean(nanmean(pv_corr(1:2,5:6)));
        end
        
        pv_corr_areas{area,subset} =  within_between_session_stability;
    end
    
end

pvalues = [];
figure
for area = 1:6
    for subset = 1:2
        current_area = pv_corr_areas{area,subset};
        pvalues(subset,area) = signrank(current_area(:,2),current_area(:,3),'Tail','right');
        
        mean_stability = nanmean(current_area);
        ste_stability = nanstd(current_area)./sqrt(size(current_area,1));
        
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
    legend(plt,{'Both','Active?1'},'Location','southwest')
    legend('boxoff')
end

% TODO - statistics

%% Figure S5J - Ensemble rate correlation across sessions for spontaneous activity
cell_cutoff = 20;
nat_movie = 1;
rate_corr_areas = {};
for area = 1:6
    valid_mice = min(calcium_excitatory_cell_count{area},[],2)>= cell_cutoff;
    
    current_area_movie = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie);
    current_area_spont = calcium_excitatory_spont_population_vectors{area}(valid_mice,nat_movie);
    
    elapsed_sess_rate_corr_spont = [];
    elapsed_sess_rate_corr_movie = [];
    for mouse = 1:length(current_area_movie)
        clc;
        disp(['Calculating ensemble rate correlation between ssessions:'])
        disp(['Dataset: Calcium imaging | Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area_movie))])
        
        current_mouse_movie = current_area_movie{mouse};
        current_mouse_spont = current_area_spont{mouse};
       
        rate_corr_movie = [];
        rate_corr_spont= [];
        for half1 = 1:6
            current_half1_movie = nanmean(current_mouse_movie(:,:,[1:5]+5*(half1-1)),3);
            current_half1_spont = nanmean(current_mouse_spont(:,:,[1:5]+5*(half1-1)),3);
            for  half2 = 1:6
                current_half2_movie = nanmean(current_mouse_movie(:,:,[1:5]+5*(half2-1)),3);
                current_half2_spont = nanmean(current_mouse_spont(:,:,[1:5]+5*(half2-1)),3);
                
                valid_cells = [nanmean(current_half1_movie,2) > 0] & [nanmean(current_half2_movie,2) > 0];
                valid_half1_activity_movie = nanmean(current_half1_movie(valid_cells,:),2);
                valid_half2_activity_movie = nanmean(current_half2_movie(valid_cells,:),2);
                
                 valid_cells = [nanmean(current_half1_spont,2) > 0] & [nanmean(current_half2_spont,2) > 0];
                valid_half1_activity_spont = nanmean(current_half1_spont(valid_cells,:),2);
                valid_half2_activity_spont = nanmean(current_half2_spont(valid_cells,:),2);
                
                rate_corr_movie(half1,half2) = corr(valid_half1_activity_movie,valid_half2_activity_movie);
                rate_corr_spont(half1,half2) = corr(valid_half1_activity_spont,valid_half2_activity_spont);
            end
        end
        
        elapsed_sess_rate_corr_movie(mouse,1) = nanmean([rate_corr_movie(1,2),rate_corr_movie(3,4),rate_corr_movie(5,6)]);
        elapsed_sess_rate_corr_movie(mouse,2) = nanmean(nanmean([rate_corr_movie(1:2,3:4),rate_corr_movie(3:4,5:6)]));
        elapsed_sess_rate_corr_movie(mouse,3) = nanmean(nanmean(rate_corr_movie(1:2,5:6)));
        
        
        elapsed_sess_rate_corr_spont(mouse,1) = nanmean([rate_corr_spont(1,2),rate_corr_spont(3,4),rate_corr_spont(5,6)]);
        elapsed_sess_rate_corr_spont(mouse,2) = nanmean(nanmean([rate_corr_spont(1:2,3:4),rate_corr_spont(3:4,5:6)]));
        elapsed_sess_rate_corr_spont(mouse,3) = nanmean(nanmean(rate_corr_spont(1:2,5:6)));
        
    end
    
    
    rate_corr_areas{1,area} = elapsed_sess_rate_corr_movie;
    rate_corr_areas{2,area} = elapsed_sess_rate_corr_spont;
    
end

figure('units','normalized','position',[0.3 0.3 0.35 0.35])
pvalues = [];
zvalues = [];
plt = [];
for area = 1:6
    subplot(2,3,area)
    mean_rate_corr_movie = nanmean(rate_corr_areas{1,area});
    ste_rate_corr_movie = nanstd(rate_corr_areas{1,area})./sqrt(size(rate_corr_areas{1,area},1));
    
    mean_rate_corr_spont = nanmean(rate_corr_areas{2,area});
    ste_rate_corr_spont = nanstd(rate_corr_areas{2,area})./sqrt(size(rate_corr_areas{2,area},1));
    
    
    [pvalues(area),~,stats] = signrank(rate_corr_areas{2,area}(:,2),rate_corr_areas{2,area}(:,3),'tail','right');
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

corrent_pvalues = bonf_holm(pvalues);
VarNames = {'area','zvalue','pvalue','bonf_holm'};
statistics = table(cell2mat(brain_areas(1:6)'),zvalues',pvalues',corrent_pvalues','VariableNames',VarNames);

clc;
disp(['Difference in the ensemble rate correlation values between proximal'])
disp(['and distal sessions during blocks of spontaneous activity'])
disp(['One-tailed Wilcoxon signed-rank tests with Holm–Bonferroni correction:'])
disp(statistics)


%% Figure S5K - PV correlation across sessions using dF/F traces
cell_cutoff = 20;
nat_movie = 1;
pv_corr_areas = {};
for area = 1:6
    valid_mice = min(calcium_excitatory_cell_count{area},[],2)>= cell_cutoff;
    
    current_area_events = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie);
    current_area_raw = calcium_excitatory_population_vectors_raw{area}(valid_mice);
    
    elapsed_sess_pv_corr_raw = [];
    elapsed_sess_pv_corr_events = [];
    for mouse = 1:length(current_area_events)
        clc;
        disp(['Calculating population vector correlation between ssessions:'])
        disp(['Dataset: Calcium imaging | Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse_events = current_area_events{mouse};
        current_mouse_raw = current_area_raw{mouse};
        
        pv_corr_events = [];
        pv_corr_raw= [];
        for half1 = 1:6
            current_half1_events = nanmean(current_mouse_events(:,:,[1:5]+5*(half1-1)),3);
            current_half1_raw = nanmean(current_mouse_raw(:,:,[1:5]+5*(half1-1)),3);
            for  half2 = 1:6
                current_half2_events = nanmean(current_mouse_events(:,:,[1:5]+5*(half2-1)),3);
                current_half2_raw = nanmean(current_mouse_raw(:,:,[1:5]+5*(half2-1)),3);
                
                valid_cells = [nanmean(current_half1_events,2) > 0] & [nanmean(current_half2_events,2) > 0];
                
                valid_half1_activity_events = current_half1_events(valid_cells,:);
                valid_half2_activity_events = current_half2_events(valid_cells,:);
                
                valid_half1_activity_raw = current_half1_raw(valid_cells,:);
                valid_half2_activity_raw = current_half2_raw(valid_cells,:);
                
                pv_corr_events(half1,half2) = nanmean(diag(corr(valid_half1_activity_events,valid_half2_activity_events)));
                pv_corr_raw(half1,half2) = nanmean(diag(corr(valid_half1_activity_raw,valid_half2_activity_raw)));
            end
        end
        
        elapsed_sess_pv_corr_events(mouse,1) = nanmean([pv_corr_events(1,2),pv_corr_events(3,4),pv_corr_events(5,6)]);
        elapsed_sess_pv_corr_events(mouse,2) = nanmean(nanmean([pv_corr_events(1:2,3:4),pv_corr_events(3:4,5:6)]));
        elapsed_sess_pv_corr_events(mouse,3) = nanmean(nanmean(pv_corr_events(1:2,5:6)));
        
        
        elapsed_sess_pv_corr_raw(mouse,1) = nanmean([pv_corr_raw(1,2),pv_corr_raw(3,4),pv_corr_raw(5,6)]);
        elapsed_sess_pv_corr_raw(mouse,2) = nanmean(nanmean([pv_corr_raw(1:2,3:4),pv_corr_raw(3:4,5:6)]));
        elapsed_sess_pv_corr_raw(mouse,3) = nanmean(nanmean(pv_corr_raw(1:2,5:6)));
        
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
    mean_pv_corr_events = nanmean(pv_corr_areas{1,area});
    std_pv_corr_events = nanstd(pv_corr_areas{1,area})./sqrt(size(pv_corr_areas{1,area},1));
    
    mean_pv_corr_raw = nanmean(pv_corr_areas{2,area});
    std_pv_corr_raw = nanstd(pv_corr_areas{2,area})./sqrt(size(pv_corr_areas{2,area},1));
    
    
    [pvalues(area),~,stats] = signrank(pv_corr_areas{2,area}(:,2),pv_corr_areas{2,area}(:,3),'tail','right');
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

VarNames = {'area','zvalue','pvalue'};
statistics = table(cell2mat(brain_areas(1:6)'),zvalues',pvalues','VariableNames',VarNames);

clc;
disp(['Difference in the PV correlation values between proximal'])
disp(['and distal sessions when using dF(t)/F0 traces'])
disp(['One-tailed Wilcoxon signed-rank tests:'])
disp(statistics)


%%  - - PV correlation across sessions using dF/F traces SCATTER PLOT
cell_cutoff = 20;
nat_movie = 1;
pv_corr_areas = {};
for area = 1:6
    valid_mice = min(calcium_excitatory_cell_count{area},[],2)>= cell_cutoff;
    
    current_area_events = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie);
    current_area_raw = calcium_excitatory_population_vectors_raw{area}(valid_mice);
    
    elapsed_sess_pv_corr_raw = [];
    elapsed_sess_pv_corr_events = [];
    for mouse = 1:length(current_area_events)
        clc;
        disp(['Calculating population vector correlation between ssessions:'])
        disp(['Dataset: Calcium imaging | Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area))])
        
        current_mouse_events = current_area_events{mouse};
        current_mouse_raw = current_area_raw{mouse};
        
        pv_corr_events = [];
        pv_corr_raw= [];
        for half1 = 1:6
            current_half1_events = nanmean(current_mouse_events(:,:,[1:5]+5*(half1-1)),3);
            current_half1_raw = nanmean(current_mouse_raw(:,:,[1:5]+5*(half1-1)),3);
            for  half2 = 1:6
                current_half2_events = nanmean(current_mouse_events(:,:,[1:5]+5*(half2-1)),3);
                current_half2_raw = nanmean(current_mouse_raw(:,:,[1:5]+5*(half2-1)),3);
                
                valid_cells = [nanmean(current_half1_events,2) > 0] & [nanmean(current_half2_events,2) > 0];
                
                valid_half1_activity_events = current_half1_events(valid_cells,:);
                valid_half2_activity_events = current_half2_events(valid_cells,:);
                
                valid_half1_activity_raw = current_half1_raw(valid_cells,:);
                valid_half2_activity_raw = current_half2_raw(valid_cells,:);
                
                pv_corr_events(half1,half2) = nanmean(diag(corr(valid_half1_activity_events,valid_half2_activity_events)));
                pv_corr_raw(half1,half2) = nanmean(diag(corr(valid_half1_activity_raw,valid_half2_activity_raw)));
            end
        end
        
        elapsed_sess_pv_corr_events(mouse,1) = nanmean([pv_corr_events(1,2),pv_corr_events(3,4),pv_corr_events(5,6)]);
        elapsed_sess_pv_corr_events(mouse,2) = nanmean(nanmean([pv_corr_events(1:2,3:4),pv_corr_events(3:4,5:6)]));
        elapsed_sess_pv_corr_events(mouse,3) = nanmean(nanmean(pv_corr_events(1:2,5:6)));
        
        
        elapsed_sess_pv_corr_raw(mouse,1) = nanmean([pv_corr_raw(1,2),pv_corr_raw(3,4),pv_corr_raw(5,6)]);
        elapsed_sess_pv_corr_raw(mouse,2) = nanmean(nanmean([pv_corr_raw(1:2,3:4),pv_corr_raw(3:4,5:6)]));
        elapsed_sess_pv_corr_raw(mouse,3) = nanmean(nanmean(pv_corr_raw(1:2,5:6)));
        
    end
    
    
    pv_corr_areas{1,area} = elapsed_sess_pv_corr_events;
    pv_corr_areas{2,area} = elapsed_sess_pv_corr_raw;
    
end


figure('units','normalized','position',[0.3 0.3 0.35 0.35])
for area = 1:6
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
mouse = 14;
area = 2;

results_path = ['E:\daniel_master\AllenBrainObservatory\calcium_imaging\results_files\excitatory2\VISl\'];
mat_list = dir([results_path,'*.mat']);
mat_list = {mat_list.name};
registration_path = ['D:\daniel-master\AllenBrainObservatory\Advanced\experimental_data\VISl\'];

load([registration_path,mat_list{mouse}(1:end-4),'\registration\aligned_data_struct.mat'])

reg_file_name = dir([registration_path,mat_list{mouse}(1:end-4),'\registration\cellRegistered*.mat']);
quality_file_name = dir([registration_path,mat_list{mouse}(1:end-4),'\registration\','*.txt']);

load([registration_path,mat_list{mouse}(1:end-4),'\registration\aligned_data_struct.mat'])
load([registration_path,mat_list{mouse}(1:end-4),'\registration\modeled_data_struct.mat'])
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
    for cell = 1:size(current_sess,1)
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

results_path = ['E:\daniel_master\AllenBrainObservatory\calcium_imaging\results_files\excitatory2\VISl\'];
mat_list = dir([results_path,'*.mat']);
mat_list = {mat_list.name};
registration_path = ['D:\daniel-master\AllenBrainObservatory\Advanced\experimental_data\VISl\'];

load([registration_path,mat_list{mouse}(1:end-4),'\registration\aligned_data_struct.mat'])

reg_file_name = dir([registration_path,mat_list{mouse}(1:end-4),'\registration\cellRegistered*.mat']);
quality_file_name = dir([registration_path,mat_list{mouse}(1:end-4),'\registration\','*.txt']);

load([registration_path,mat_list{mouse}(1:end-4),'\registration\aligned_data_struct.mat'])
load([registration_path,mat_list{mouse}(1:end-4),'\registration\modeled_data_struct.mat'])
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
    for cell = 1:size(current_sess,1)
        clc;
        disp(['Registration of neuronal responses across sessions using Sheintuch et al. (2017) method:'])
        disp(['Stimulus: Natural movie 1 | Area: ',brain_areas{area},' | Mouse: ',num2str(mouse),'/',num2str(length(current_area)),...
            ' | Session: ',num2str(session),'/3 | Cell: ',num2str(cell),'/',num2str(size(current_sess,1))])
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

clearvars -except neuropixels_population_vectors neuropixels_drifting_gratings brain_areas...
    neuropixels_running_speed neuropixels_pupil_size neuropixels_cell_count movie_repeats...
    calcium_excitatory_population_vectors calcium_excitatory_drifting_gratings...
    calcium_excitatory_spont_population_vectors calcium_excitatory_cell_count...
    calcium_excitatory_imaging_depth calcium_excitatory_running_speed ...
    calcium_excitatory_pupil_size calcium_excitatory_population_vectors_raw...
    calcium_inhibitory_population_vectors calcium_inhibitory_cell_count...
    calcium_inhibitory_cre_line calcium_excitatory_sorted_mouse_age...
    neuropixels_population_vectors_tsne

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

area = 2;
nat_movie = 1;
cell_cutoff = 20;
valid_mice = min(calcium_excitatory_cell_count{area},[],2) >= cell_cutoff;
current_area = calcium_excitatory_population_vectors{area}(valid_mice,nat_movie);

results_path = ['E:\daniel_master\AllenBrainObservatory\calcium_imaging\results_files\excitatory2\VISl\'];
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
nat_movie = 1;
area = 1;
valid_mice = movie_repeats(:,nat_movie) == 30;
V1_pseudo = cell2mat(neuropixels_population_vectors_tsne(valid_mice,area,nat_movie));
V1_pseudo = V1_pseudo(:,:,31:60);


load('FigureS7A_neuropixels.mat','state')
rng(state)
cells_rand_ind = randperm(size(V1_pseudo,1));
cell_num_list = [75,100,125,150,200,250,500,1000,1500];

neuropixels_dim_reduceV1_cell_count = {};
for cell_included = 1:length(cell_num_list)
    clc;
     disp('Calculating structures:')
    disp(['Cells: ',num2str(cell_num_list(cell_included)), ' | ',...
        num2str(cell_included),'/',num2str(length(cell_num_list)),...
        ' | ', num2str([cell_included./length(cell_num_list)]*100),'%'])
  
    cell_num = cells_rand_ind(1:cell_num_list(cell_included));
    subset_V1_pseudo = V1_pseudo(cell_num,:,:);
    V1_pseudo_strcut = reshape(subset_V1_pseudo,[size(subset_V1_pseudo,1),90*30]);
    
    dim_reduceV1 = tsne(V1_pseudo_strcut','Algorithm','exact','Distance','cosine',...
        'NumDimensions',2,'NumPCAComponents',20,'Perplexity',200);
   neuropixels_dim_reduceV1_cell_count{cell_included} = dim_reduceV1;
end

xy_list = [2,1;1 2;1 2;2 1;1 2;1 2;1 2;1 2;1 2];
direction_list = [1 1;1 1;-1 1;1 -1; -1 -1;1 -1;-1 -1;1 -1; -1 -1];

figure('units','normalized','position',[0.3 0.2 0.3 0.525])
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

area = 1;
V1_pseudo = cell2mat(calcium_excitatory_population_vectors{area}(:,3));
valid_cells = (nanmean(nanmean(V1_pseudo(:,:,1:10),3),2) > 0 &...
    nanmean(nanmean(V1_pseudo(:,:,11:20),3),2) > 0 &...
    nanmean(nanmean(V1_pseudo(:,:,21:30),3),2) > 0);
V1_pseudo = V1_pseudo(valid_cells,:,:);
V1_pseudo(V1_pseudo==0) = 0.000001;
 
load('FigureS7A_calcium.mat','state')
rng(state)

cells_rand_ind = randperm(size(V1_pseudo,1));
cell_num_list = [75,100,125,150,200,250,500,1000,1500];


calcium_dim_reduceV1_cell_count = {};
for cell_included = 1:length(cell_num_list)
    clc;
     disp('Calculating structures:')
    disp(['Cells: ',num2str(cell_num_list(cell_included)), ' | ',...
        num2str(cell_included),'/',num2str(length(cell_num_list)),...
        ' | ', num2str([cell_included./length(cell_num_list)]*100),'%'])
  
    cell_num = cells_rand_ind(1:cell_num_list(cell_included));
    subset_V1_pseudo = V1_pseudo(cell_num,:,:);
    V1_pseudo_strcut = reshape(subset_V1_pseudo,[size(subset_V1_pseudo,1),90*30]);
    %,'Algorithm','exact'
    dim_reduceV1 = tsne(V1_pseudo_strcut','Distance','cosine',...
        'NumDimensions',2,'NumPCAComponents',20,'Perplexity',200);
   calcium_dim_reduceV1_cell_count{cell_included} = dim_reduceV1;
end

xy_list = [2 1;1 2;1 2;2 1;2 1;2 1;2 1;1 2;1 2];
direction_list = [-1 -1;1 -1;-1 1;1 -1; 1 -1;-1 1;1 1;1 -1; 1 1];

figure('units','normalized','position',[0.3 0.2 0.3 0.525])
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
nat_movie = 1;
area = 1;
valid_mice = movie_repeats(:,nat_movie) == 30;
V1_pseudo = cell2mat(neuropixels_population_vectors_tsne(valid_mice,area,nat_movie));
V1_pseudo = V1_pseudo(:,:,31:60);


load('FigureS7A_neuropixels.mat','state')
rng(state)
cells_rand_ind = randperm(size(V1_pseudo,1));
cell_num_list = [75,100,125,150,200,250,500,1000,1500];

neuropixels_dim_reduceV1_cell_count = {};
for cell_included = 1:length(cell_num_list)
    clc;
     disp('Calculating structures:')
    disp(['Cells: ',num2str(cell_num_list(cell_included)), ' | ',...
        num2str(cell_included),'/',num2str(length(cell_num_list)),...
        ' | ', num2str([cell_included./length(cell_num_list)]*100),'%'])
  
    cell_num = cells_rand_ind(1:cell_num_list(cell_included));
    subset_V1_pseudo = V1_pseudo(cell_num,:,:);
    V1_pseudo_strcut = reshape(subset_V1_pseudo,[size(subset_V1_pseudo,1),90*30]);
    
    [~,dim_reduceV1] = pca(V1_pseudo_strcut');
   neuropixels_dim_reduceV1_cell_count{cell_included} = dim_reduceV1;
end

xy_list = [1,2;1 2;1 2;2 1;1 2;1 2;1 2;1 2;1 2];
direction_list = [1 1;1 1;-1 1;1 -1; -1 -1;1 -1;-1 -1;1 -1; -1 -1];

figure('units','normalized','position',[0.3 0.2 0.3 0.525])
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


area = 1;
V1_pseudo = cell2mat(calcium_excitatory_population_vectors{area}(:,3));
valid_cells = (nanmean(nanmean(V1_pseudo(:,:,1:10),3),2) > 0 &...
    nanmean(nanmean(V1_pseudo(:,:,11:20),3),2) > 0 &...
    nanmean(nanmean(V1_pseudo(:,:,21:30),3),2) > 0);
V1_pseudo = V1_pseudo(valid_cells,:,:);
V1_pseudo(V1_pseudo==0) = 0.000001;
 
load('FigureS7A_calcium.mat','state')
rng(state)

cells_rand_ind = randperm(size(V1_pseudo,1));
cell_num_list = [75,100,125,150,200,250,500,1000,1500];

calcium_dim_reduceV1_cell_count = {};
for cell_included = 1:length(cell_num_list)
    clc;
     disp('Calculating structures:')
    disp(['Cells: ',num2str(cell_num_list(cell_included)), ' | ',...
        num2str(cell_included),'/',num2str(length(cell_num_list)),...
        ' | ', num2str([cell_included./length(cell_num_list)]*100),'%'])
  
    cell_num = cells_rand_ind(1:cell_num_list(cell_included));
    subset_V1_pseudo = V1_pseudo(cell_num,:,:);
    V1_pseudo_strcut = reshape(subset_V1_pseudo,[size(subset_V1_pseudo,1),90*30]);
    
    [~,dim_reduceV1] = pca(V1_pseudo_strcut');
   calcium_dim_reduceV1_cell_count{cell_included} = dim_reduceV1;
end

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
nat_movie = 1;
subset_population_vectors = {};
for area = 1:6
    for mouse = 1:size(calcium_excitatory_population_vectors{area},1)
        subset_population_vectors{area} = calcium_excitatory_population_vectors{area}(:,1);
    end
end

num_shuffles = 3000;
cell_cutoff_list = [5,10,25,50,75,100,150,200,250,300,350,400,500,600,700,800,1000,1200,1400];
between_similarity_acc_all_cutoffs = [];
between_similarity_acc_shuffled_all_cutoffs = [];
for cells_included = length(cell_cutoff_list)+1
    
    between_similarity_acc = nan(num_shuffles,6);
    between_similarity_acc_shuffled = nan(num_shuffles,6);
    
    for shuffle = 1:num_shuffles
        clc;[cells_included,shuffle]
        pseudo_area_cell_num = 0;
        if cells_included <= length(cell_cutoff_list)
            breaker = cell_cutoff_list(cells_included);
        else
            breaker = 1;
        end
        
        while   breaker > min(pseudo_area_cell_num(:))
            
            pseudo_mouseA = {};
            pseudo_mouseB = {};
            for area = 1:6
                pseudo_mouseA_ind = sort(randperm(size(subset_population_vectors{area},1),round(size(subset_population_vectors{area},1)./2)));
                pseudo_mouseB_ind = find(~ismember([1:size(subset_population_vectors{area},1)],pseudo_mouseA_ind));
                
                pseudo_mouseA{area} = cell2mat(subset_population_vectors{area}(pseudo_mouseA_ind));
                pseudo_mouseB{area} = cell2mat(subset_population_vectors{area}(pseudo_mouseB_ind));
                pseudo_area_cell_num(1,area) = size(pseudo_mouseA{area},1);
                pseudo_area_cell_num(2,area) = size(pseudo_mouseB{area},1);
            end
        end
        if cells_included <= length(cell_cutoff_list)
            min_cell_num = cell_cutoff_list(cells_included);
        else
            min_cell_num = min(pseudo_area_cell_num(:));
             min_cell_num_all_shuffles(shuffle) = min_cell_num;
        end
        
        pseudo_mouseA_subset = {};
        pseudo_mouseB_subset = {};
        for area = 1:6
            subset_cell_ids_mouseA = sort(randperm(pseudo_area_cell_num(1,area),min_cell_num));
            pseudo_mouseA_subset{area} = pseudo_mouseA{area}(subset_cell_ids_mouseA,:,:);
            
            subset_cell_ids_mouseB = sort(randperm(pseudo_area_cell_num(2,area),min_cell_num));
            pseudo_mouseB_subset{area} = pseudo_mouseB{area}(subset_cell_ids_mouseB,:,:);
        end
        
        all_cells_pseudo_mouseA = cell2mat(pseudo_mouseA');
        all_cells_pseudo_mouseB = cell2mat(pseudo_mouseB');
        rand_cells_id_mouseA = randperm(size(all_cells_pseudo_mouseA,1));
        rand_cells_id_mouseB = randperm(size(all_cells_pseudo_mouseB,1));
        
        shuffled_pseudo_mouseA_subset = {};
        shuffled_pseudo_mouseB_subset = {};
        for area = 1:6
            current_pseudo_area = [1:min_cell_num] + min_cell_num*(area-1);
            
            shuffled_pseudo_mouseA_subset{area} = all_cells_pseudo_mouseA(rand_cells_id_mouseA(current_pseudo_area),:,:);
            shuffled_pseudo_mouseB_subset{area} = all_cells_pseudo_mouseB(rand_cells_id_mouseB(current_pseudo_area),:,:);
        end
        
      
        internal_structures_pseudoA = [];
        internal_structures_pseudoB = [];
        internal_structures_pseudoA_shuffle = [];
        internal_structures_pseudoB_shuffle = [];
        internal_structures_labels = [];
        triu_ind = boolean(triu(ones(30),1));
        for area = 1:6
            current_structure_mouseA = [];
            current_structure_mouseA(:,:,1) = corr(nanmean(pseudo_mouseA_subset{area}(:,:,1:10),3));
            current_structure_mouseA(:,:,2) = corr(nanmean(pseudo_mouseA_subset{area}(:,:,11:20),3));
            current_structure_mouseA(:,:,3) = corr(nanmean(pseudo_mouseA_subset{area}(:,:,21:30),3));
            current_structure_mouseA = nanmean(current_structure_mouseA,3);
            internal_structures_pseudoA(:,area) = current_structure_mouseA(triu_ind);
            
            current_structure_mouseB = [];
            current_structure_mouseB(:,:,1) = corr(nanmean(pseudo_mouseB_subset{area}(:,:,1:10),3));
            current_structure_mouseB(:,:,2) = corr(nanmean(pseudo_mouseB_subset{area}(:,:,11:20),3));
            current_structure_mouseB(:,:,3) = corr(nanmean(pseudo_mouseB_subset{area}(:,:,21:30),3));
            current_structure_mouseB = nanmean(current_structure_mouseB,3);
            internal_structures_pseudoB(:,area) = current_structure_mouseB(triu_ind);
            
            current_structure_mouseA_shuffle = [];
            current_structure_mouseA_shuffle(:,:,1) = corr(nanmean(shuffled_pseudo_mouseA_subset{area}(:,:,1:10),3));
            current_structure_mouseA_shuffle(:,:,2) = corr(nanmean(shuffled_pseudo_mouseA_subset{area}(:,:,11:20),3));
            current_structure_mouseA_shuffle(:,:,3) = corr(nanmean(shuffled_pseudo_mouseA_subset{area}(:,:,21:30),3));
            current_structure_mouseA_shuffle = nanmean(current_structure_mouseA_shuffle ,3);
            internal_structures_pseudoA_shuffle(:,area) = current_structure_mouseA_shuffle(triu_ind);
            
            
            current_structure_mouseB_shuffle = [];
            current_structure_mouseB_shuffle(:,:,1) = corr(nanmean(shuffled_pseudo_mouseB_subset{area}(:,:,1:10),3));
            current_structure_mouseB_shuffle(:,:,2) = corr(nanmean(shuffled_pseudo_mouseB_subset{area}(:,:,11:20),3));
            current_structure_mouseB_shuffle(:,:,3) = corr(nanmean(shuffled_pseudo_mouseB_subset{area}(:,:,21:30),3));
            current_structure_mouseB_shuffle = nanmean(current_structure_mouseB_shuffle ,3);
            internal_structures_pseudoB_shuffle(:,area) = current_structure_mouseB_shuffle(triu_ind);
    
        end
        
        permutations = flipud(perms([1:6]));
        
        between_mouse_similarity = corr(internal_structures_pseudoA,internal_structures_pseudoB,'rows','complete');
        between_mouse_similarity_shuffled = corr(internal_structures_pseudoA_shuffle,internal_structures_pseudoB_shuffle,'rows','complete');
        
        similarity_sum = [];
        similarity_sum_shuffled = [];
        for perm = 1:size(permutations,1)
            similarity_sum(perm) = sum(diag(between_mouse_similarity(:,permutations(perm,:))));
            similarity_sum_shuffled(perm) = sum(diag(between_mouse_similarity_shuffled(:,permutations(perm,:))));
        end
    
        [B,I] = max(similarity_sum);
        between_similarity_acc(shuffle,:) = permutations(I,:) == [1:6];
        rank_between(shuffle) = I;
        
        [B,I] = max(similarity_sum_shuffled);
        between_similarity_acc_shuffled(shuffle,:) = permutations(I,:) == [1:6];
    end
    between_similarity_acc_all_cutoffs(cells_included,:) = sum(between_similarity_acc)./num_shuffles;
    between_similarity_acc_shuffled_all_cutoffs(cells_included,:) = sum(between_similarity_acc_shuffled)./num_shuffles;
end


plt = [];
xvalues = [cell_cutoff_list,round(nanmean(min_cell_num_all_shuffles))];
figure('units','normalized','position',[0.35 0.4 0.25 0.375])
for area = 1:6
    hold on
    plt(area) = plot(between_similarity_acc_all_cutoffs(:,area)*100,'color',colors(area,:),'linewidth',3);
    plot(between_similarity_acc_shuffled_all_cutoffs(:,area)*100,'color',[0.2 0.2 0.2]+0.1*(area-1),'linewidth',3)
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
nat_movie = 1;
subset_population_vectors = {};
for area = 1:6
    for mouse = 1:size(neuropixels_population_vectors,1)
        if ~isempty(neuropixels_population_vectors{mouse,area,nat_movie})
            subset_population_vectors{mouse,area} = neuropixels_population_vectors{mouse,area,nat_movie}(:,:,1:20);
        end
    end
end

num_shuffles = 3000;
cell_cutoff_list = [5,10,25,50,75,100,150,200,250,300,350,400,500,600,700];
between_similarity_acc_all_cutoffs = [];
between_similarity_acc_shuffled_all_cutoffs = [];

min_cell_num_all_shuffles = [];
for cells_included = 1:length(cell_cutoff_list)+1
    between_similarity_acc_shuffled = nan(num_shuffles,6);
    
    for shuffle = 1:num_shuffles
        clc;[cells_included,shuffle]
        pseudo_area_cell_num = 0;
        if cells_included <= length(cell_cutoff_list)
            breaker = cell_cutoff_list(cells_included);
        else
            breaker = 1;
        end
        
        while   breaker > min(pseudo_area_cell_num(:))
            pseudo_mouseA_ind = sort(randperm(size(subset_population_vectors,1),size(subset_population_vectors,1)./2));
            pseudo_mouseB_ind = find(~ismember([1:size(subset_population_vectors,1)],pseudo_mouseA_ind));
            
            
            pseudo_mouseA = {};
            pseudo_mouseB = {};
            for area = 1:6
                pseudo_mouseA{area} = cell2mat(subset_population_vectors(pseudo_mouseA_ind,area,nat_movie));
                pseudo_mouseB{area} = cell2mat(subset_population_vectors(pseudo_mouseB_ind,area,nat_movie));
                pseudo_area_cell_num(1,area) = size(pseudo_mouseA{area},1);
                pseudo_area_cell_num(2,area) = size(pseudo_mouseB{area},1);
            end
        end
        
        if cells_included <= length(cell_cutoff_list)
            min_cell_num = cell_cutoff_list(cells_included);
        else
            min_cell_num = min(pseudo_area_cell_num(:));
            min_cell_num_all_shuffles(shuffle) = min_cell_num;
        end
        
        pseudo_mouseA_subset = {};
        pseudo_mouseB_subset = {};
        for area = 1:6
            subset_cell_ids_mouseA = sort(randperm(pseudo_area_cell_num(1,area),min_cell_num));
            pseudo_mouseA_subset{area} = pseudo_mouseA{area}(subset_cell_ids_mouseA,:,:);
            
            subset_cell_ids_mouseB = sort(randperm(pseudo_area_cell_num(2,area),min_cell_num));
            pseudo_mouseB_subset{area} = pseudo_mouseB{area}(subset_cell_ids_mouseB,:,:);
        end

        internal_structures_pseudoA_shuffle = [];
        internal_structures_pseudoB_shuffle = [];
        
        internal_structures_labels = [];
        triu_ind = boolean(triu(ones(30),1));
        for area = 1:6
            rand_ind = circshift([1:30],randperm(30,1));
            current_structure_mouseA_shuffle = corr(nanmean(pseudo_mouseA_subset{area}(:,rand_ind,:),3));
            internal_structures_pseudoA_shuffle(:,area) = current_structure_mouseA_shuffle(triu_ind);
            
            rand_ind = circshift([1:30],randperm(30,1));
            current_structure_mouseB_shuffle = corr(nanmean(pseudo_mouseB_subset{area}(:,rand_ind,:),3));
            internal_structures_pseudoB_shuffle(:,area) = current_structure_mouseB_shuffle(triu_ind);
        end
        
        permutations = flipud(perms([1:6]));
        
        between_mouse_similarity_shuffled = corr(internal_structures_pseudoA_shuffle,internal_structures_pseudoB_shuffle);

        similarity_sum_shuffled = [];
        
        for perm = 1:size(permutations,1)
            similarity_sum_shuffled(perm) = sum(diag(between_mouse_similarity_shuffled(:,permutations(perm,:))));
        end
      
        [B,I] = max(similarity_sum_shuffled);
        between_similarity_acc_shuffled(shuffle,:) = permutations(I,:) == [1:6];
        
    end
     between_similarity_acc_shuffled_all_cutoffs(cells_included,:) = sum(between_similarity_acc_shuffled)./num_shuffles;
       
end

plt = [];
xvalues = [cell_cutoff_list,round(nanmean(min_cell_num_all_shuffles))];
figure('units','normalized','position',[0.35 0.4 0.25 0.375])
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
subset_population_vectors = cell(58,6,2);
repeat_num = 30;
for nat_movie = 1:2
    for area = 1:6
        for mouse = 1:size(neuropixels_population_vectors_tsne,1)
            if ~isempty(neuropixels_population_vectors_tsne{mouse,area,nat_movie}) && movie_repeats(mouse,1) == repeat_num
                subset_population_vectors{mouse,area,nat_movie} = neuropixels_population_vectors_tsne{mouse,area,nat_movie}(:,:,1:10);
            end
        end
    end
end

dim_reduce_all = {};
cells_ind = [];
load('figureS7G.mat','state')
rng(state)
  figure
for area = 1:6
    current_area = cell2mat(subset_population_vectors(:,area,1));
    cells_ind(area,:) = randperm(size(current_area,1),919);
    current_area = current_area(cells_ind(area,:),:,:);
    current_area_struct = reshape(current_area,[size(current_area,1),size(current_area,2)*size(current_area,3)]);
    clc;[area,1]
    dim_reduce = tsne(current_area_struct','Algorithm','exact','Distance','cosine',...
        'NumDimensions',3,'NumPCAComponents',20,'Perplexity',100);
    dim_reduce_all{1,area} = dim_reduce;

    current_area = cell2mat(subset_population_vectors(:,area,2));
    current_area_struct = reshape(current_area,[size(current_area,1),size(current_area,2)*size(current_area,3)]);
    clc;[area,2]
    dim_reduce = tsne(current_area_struct','Algorithm','exact','Distance','cosine',...
        'NumDimensions',3,'NumPCAComponents',20,'Perplexity',100);
    dim_reduce_all{2,area} = dim_reduce;
end

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
        title('Natural movie 1 (original frame sequance)')
    end
    subplot(2,6,area+6)
    scatter3(dim_reduce_all{2,area}(:,1),dim_reduce_all{2,area}(:,2),dim_reduce_all{2,area}(:,3),20,repmat([1:90],[1 10]),'filled')
    set(gca,'xtick',[],'ytick',[])
    view(view_list_shuffled(area,:))
    
    grid off
    if area == 3
        title('Shuffled natural movie 1 (fixed temporally shuffled frame sequance)')
    end

end
colormap(new_jet_colormap)

%% Figure S7F - Between pseudo-mice decoder using NM1 compared to SNM1

subset_population_vectors = cell(58,6,2);
repeat_num = 30;
for nat_movie = 1:2
    for area = 1:6
        for mouse = 1:size(neuropixels_population_vectors,1)
            if ~isempty(neuropixels_population_vectors{mouse,area,nat_movie}) && movie_repeats(mouse,1) == repeat_num
                if  nat_movie == 1
                    subset_population_vectors{mouse,area,nat_movie} = neuropixels_population_vectors{mouse,area,nat_movie}(:,:,21:40);
                    
                elseif nat_movie == 2
                    subset_population_vectors{mouse,area,nat_movie} = neuropixels_population_vectors{mouse,area,nat_movie}(:,:,1:20);
                end
            end
        end
    end
end

accurate_classification_movA = [];
accurate_classification_movB = [];
accurate_classification_across = [];
num_shuffles = 1000;
for shuffle = 1:num_shuffles
    clc;[shuffle]
    all_mice_ind = find(movie_repeats(:,1) == repeat_num);
    pseudo_mouseA_ind = all_mice_ind(sort(randperm(length(all_mice_ind),round(length(all_mice_ind)./2))));
    pseudo_mouseB_ind = all_mice_ind(~ismember(all_mice_ind,pseudo_mouseA_ind));
    
    cell_num_all_areas = [];
    pseudo_mouseA = {};
    pseudo_mouseB = {};
    for nat_movie = 1:2
        for area = 1:6
            pseudo_mouseA{nat_movie,area} = cell2mat(subset_population_vectors(pseudo_mouseA_ind,area,nat_movie));
            pseudo_mouseB{nat_movie,area} = cell2mat(subset_population_vectors(pseudo_mouseB_ind,area,nat_movie));
            cell_num_all_areas(1,area) = size(pseudo_mouseA{nat_movie,area},1);
            cell_num_all_areas(2,area) = size(pseudo_mouseB{nat_movie,area},1);
        end
    end
    
    min_num_cell = min(cell_num_all_areas(:));
    subset_pseudo_mouseA = {};
    subset_pseudo_mouseB = {};
    triu_ind = boolean(triu(ones(size(pseudo_mouseA{1,1},2)),1));
    movA_all_structures_pseudoA = [];
    movB_all_structures_pseudoA = [];
    movA_all_structures_pseudoB = [];
    movB_all_structures_pseudoB = [];
    for nat_movie = 1:2
        for area = 1:6
            rand_cellsA = randperm(size(pseudo_mouseA{nat_movie,area},1),min_num_cell);
            rand_cellsB = randperm(size(pseudo_mouseB{nat_movie,area},1),min_num_cell);
            subset_pseudo_mouseA{nat_movie,area} = pseudo_mouseA{nat_movie,area}(rand_cellsA,:,:);
            subset_pseudo_mouseB{nat_movie,area} = pseudo_mouseB{nat_movie,area}(rand_cellsB,:,:);
            
            pseudo_mouseA_struct = corr(nanmean(subset_pseudo_mouseA{nat_movie,area},3));
            pseudo_mouseB_struct = corr(nanmean(subset_pseudo_mouseB{nat_movie,area},3));
              
            if nat_movie ==1
                movA_all_structures_pseudoA(:,area) = pseudo_mouseA_struct(triu_ind);
                movA_all_structures_pseudoB(:,area) = pseudo_mouseB_struct(triu_ind);
            elseif nat_movie ==2
                movB_all_structures_pseudoA(:,area) = pseudo_mouseA_struct(triu_ind);
                movB_all_structures_pseudoB(:,area) = pseudo_mouseB_struct(triu_ind);
            end
        end
    end
    
    similarity_between_mice_movA = corr(movA_all_structures_pseudoA,movA_all_structures_pseudoB);
    similarity_between_mice_movB = corr(movB_all_structures_pseudoA,movB_all_structures_pseudoB);
    
    permutations = flipud(perms([1:6]));
    perms_similarity_movA = [];
    perms_similarity_movB = [];

    for perm = 1:size(permutations,1)
        current_perm = permutations(perm,:);
        perms_similarity_movA(perm) = sum(diag(corr(movA_all_structures_pseudoA,movA_all_structures_pseudoB(:,current_perm))));
        perms_similarity_movB(perm) = sum(diag(corr(movB_all_structures_pseudoA,movB_all_structures_pseudoB(:,current_perm))));
    end
    
    [~,best_perm_movA] = max(perms_similarity_movA);
    [~,best_perm_movB] = max(perms_similarity_movB);

    accurate_classification_movA(shuffle,:) = permutations(best_perm_movA,:);
    accurate_classification_movB(shuffle,:) = permutations(best_perm_movB,:);

end


acc_NM1 = 100*sum(accurate_classification_movA==[1:6])./shuffle;
acc_SNM1 = 100*sum(accurate_classification_movB==[1:6])./shuffle;

figure
hold on
plt = [];
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

%% Figure S7I - Internal structure stability VS PV stability when shuffling cells’ identities
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

