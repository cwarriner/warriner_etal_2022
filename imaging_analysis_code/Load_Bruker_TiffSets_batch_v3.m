%Load_Bruker_TiffSets_batch
%processes Bruker image time series from folders of individual tiffs, as they are stored on the acquisition machine

clear all

%%% Params
%First provide the folder containing the individual time series folders
raw_data_folder = 'G:\Data\CSMNs\Imaging Raw';

%Then provide the folder to which processed movies will be saved
save_folder = 'G:\Data\CSMNs\Imaging processed';

activity_ch = 3; %channel in which functional indicator is imaged
save_raw = 1; %set to 1 if you want to save the original movie
play_during_loading = 0; %set to 1 if you want to playback during loading
motion_correct = 1; %set to 1 if you want to motion correct
save_rigid_corrected = 0; %set to 1 if you want to save the rigid motion corrected movie
save_nonrigid_corrected = 1; %set to 1 if you want to save the non-rigid motion corrected movie
playback_frametime = 0.025; %in seconds - 0.025 works well

%Initialize
close all
num_movies = 0;
process = 0;

main_folder = dir(raw_data_folder);
cd(raw_data_folder)

for i = 1:numel(main_folder)
    dir_name = strsplit(main_folder(i).name,'-');
    if strcmp(dir_name(1),'TSeries')  % check if element in main_folder is a folder (that may have tiffs in it)
        sub_folder = dir(main_folder(i).name);
        process = 0; %reset process state

        %%% Check to see if folder has tiffs in it - this is needed if there is a chance there are other folders not with time series in them
        for j = 1:numel(sub_folder)
            [~,name,ext] = fileparts(sub_folder(j).name);
            if strcmp(ext,'.tif') && ~process
                process = 1;
            end
        end

        if process
            %%% Get params from the xml file - need to avoid the voltage recording xml files
            cd(sub_folder(1).folder)
            xml_file = dir('*.xml');
            for j = 1:numel(xml_file)
                xml_name = strsplit(xml_file(j).name,'.');
                if strcmp(xml_name(1),main_folder(i).name) % the xml file we want should have the same name as the tseries folder
                    xml_struct = xml2struct(xml_file(j).name);
                    mov_params.num_frames = size(xml_struct.PVScan.Sequence.Frame,2);
                    mov_params.frame_period = str2num(xml_struct.PVScan.PVStateShard.PVStateValue{6}.Attributes.value);
                    mov_params.lines_per_frame = str2num(xml_struct.PVScan.PVStateShard.PVStateValue{9}.Attributes.value);
                    mov_params.pixels_per_line = str2num(xml_struct.PVScan.PVStateShard.PVStateValue{17}.Attributes.value);
                    mov_params.microns_per_pixel_x = str2num(xml_struct.PVScan.PVStateShard.PVStateValue{11}.IndexedValue{1}.Attributes.value);
                    mov_params.microns_per_pixel_y = str2num(xml_struct.PVScan.PVStateShard.PVStateValue{11}.IndexedValue{2}.Attributes.value);
                    mov_params.x_coord = str2num(xml_struct.PVScan.PVStateShard.PVStateValue{19}.SubindexedValues{1}.SubindexedValue.Attributes.value);
                    mov_params.y_coord = str2num(xml_struct.PVScan.PVStateShard.PVStateValue{19}.SubindexedValues{2}.SubindexedValue.Attributes.value);
                    mov_params.z_coord = str2num(xml_struct.PVScan.PVStateShard.PVStateValue{19}.SubindexedValues{3}.SubindexedValue{1}.Attributes.value);
                end
            end

            %%% Determine which channels were recorded
            num_channels = size(xml_struct.PVScan.Sequence.Frame{1}.File,2);
            channels = zeros(1,num_channels);
            for j = 1:num_channels
                channels(j) = str2num(xml_struct.PVScan.Sequence.Frame{1}.File{j}.Attributes.channel);
            end

            %%% Load the tiffs
            for j = 1:num_channels
                num_movies = num_movies + 1; %just for figure indexing
                mov = zeros(mov_params.lines_per_frame,mov_params.pixels_per_line,mov_params.num_frames); %Initialize movie array
                for k = 1:mov_params.num_frames
                    mov(:,:,k) = loadtiff(xml_struct.PVScan.Sequence.Frame{k}.File{j}.Attributes.filename);
                    if play_during_loading
                        figure(num_movies)
                        colormap('bone')
                        imagesc(mov(:,:,k))
                        axis equal; axis tight;
                        pause(playback_frametime)
                    end
                end

                save_name = [main_folder(i).name '_mov_ch' num2str(channels(j))];
                if save_raw
                    cd(save_folder)
                    save(save_name,'mov','mov_params');
                    cd(sub_folder(1).folder)
                end
                figure(num_movies)
                colormap('bone')
                imagesc(mean(mov,3))
                axis equal; axis tight;
                title(['Channel ' num2str(channels(j)) ' Mean Image'])

                if motion_correct && channels(j) == activity_ch
                    Y = single(mov); % convert to single precision 
                    %Y = Y - min(Y(:));

                    %normcorre params - NB: current version of normcorre_batch giving output that sometimes doesn't appear to have the right bidirectional scanning correction
                    %batch parallelizes the routine for faster speed
                    options_rigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'bin_width',50,'max_shift',15,'us_fac',50);
                    options_nonrigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'grid_size',[32,32],'mot_uf',4,'bin_width',50,'max_shift',15,'max_dev',3,'us_fac',50);

                    %%% perform rigid motion correction
                    tic; [mov_rigid_corr,shifts1,template1] = normcorre(Y,options_rigid); toc

                    %%% now try non-rigid motion correction (also in parallel)
                    tic; [mov_nonrigid_corr,shifts2,template2] = normcorre(Y,options_nonrigid); toc %got rid of batch because shifts seem messed up

                    %%% compute metrics 
                    nnY = quantile(Y(:),0.005);
                    mmY = quantile(Y(:),0.995);

                    [corr_params.cY,corr_params.mY,corr_params.vY] = motion_metrics(Y,10);
                    [corr_params.cM1,corr_params.mM1,corr_params.vM1] = motion_metrics(mov_rigid_corr,10);
                    [corr_params.cM2,corr_params.mM2,corr_params.vM2] = motion_metrics(mov_nonrigid_corr,10);
                    T = length(corr_params.cY);

                    shifts_rigid = squeeze(cat(3,shifts1(:).shifts));
                    shifts_nr = cat(ndims(shifts2(1).shifts)+1,shifts2(:).shifts);
                    shifts_nr = reshape(shifts_nr,[],ndims(Y)-1,T);
                    shifts_x_nonrigid = squeeze(shifts_nr(:,1,:))';
                    shifts_y_nonrigid = squeeze(shifts_nr(:,2,:))';

                    %%% plot a movie with the results
                    figure(100+num_movies);
                    for t = 1:1:T
                        subplot(121);imagesc(Y(:,:,t),[nnY,mmY]); xlabel('raw data','fontsize',14,'fontweight','bold'); axis equal; axis tight;
                        title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); colormap('bone')
                        subplot(122);imagesc(mov_nonrigid_corr(:,:,t),[nnY,mmY]); xlabel('non-rigid corrected','fontsize',14,'fontweight','bold'); axis equal; axis tight;
                        title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); colormap('bone')
                        set(gca,'XTick',[],'YTick',[]);
                        drawnow;
                        pause(playback_frametime);
                    end
                    if save_rigid_corrected
                        cd(save_folder)
                        save(save_name,'mov_rigid_corr','corr_params','shifts_rigid','-append');
                        cd(sub_folder(1).folder)
                    end
                    if save_nonrigid_corrected
                        cd(save_folder)
                        save(save_name,'mov_nonrigid_corr','corr_params','shifts_x_nonrigid','shifts_y_nonrigid','-append');
                        cd(sub_folder(1).folder)
                    end

                end
            end
        end
    end
    cd(raw_data_folder)
end
