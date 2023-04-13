function xsectionsTransformData(file_folder, out_folder, time_interval, workers, varargin)
%{
    CrossSections

    Function to produce exploratory images of cross sections of transformed
    image volumes. 

    INPUTS
    file_folder: folder containing the raw images
    out_folder: folder to contain the transformed images
    time_interval: array specifying the first and last images
    workers: number of cores for the transformation

    OPTIONAL ARGUMENTS
    double_scan: bool. if the experiment is double scan
    scan_order: order of scans. Left:2, Right:1
    num_scan: number of exploratory slices

    Guillermo Serrano Najera 2019

%}
    %% Parse arguments
    opt_arguments = {
                     1,... % double_scan
                     [1;2],... % order of scans. Left:2, Right:1
                     5,... % number of scans
                     };
    
    opt_arguments(1:nargin-4) = varargin;
    double_scan = opt_arguments{1};
    scan_order = opt_arguments{2};
    num_scan = opt_arguments{3};
    
    % make output folders
    mkdir([out_folder filesep 'img']);

    % Read stack sizes from info file.
    Width=findMatchingNumber([file_folder filesep '_info.txt'],{'Width: ' '%d'},1);
    Height=findMatchingNumber([file_folder filesep '_info.txt'],{'Height: ' '%d'},1);
    Depth=findMatchingNumber([file_folder filesep '_info.txt'],{'Depth: ' '%d'},1);

    %% start parallel pool 
    pool = gcp('nocreate');
    if isempty(pool)
       parpool(workers);
    elseif pool.NumWorkers~=workers
       delete(pool);
       parpool(workers);
    end

    % split time points for the parallel workers for the images to be produced in order
    time_points_to_compute=cell(workers,1);
    if double_scan
        for i_ind=1:workers
            vector = (i_ind*2-1+time_interval(1)-1):2*workers:time_interval(2);
            time_points_to_compute{i_ind}=sort([vector,vector+1]);
        end
    else
        for i_ind=1:workers
            time_points_to_compute{i_ind}=(i_ind+time_interval(1)-1):workers:time_interval(2);
        end
    end
    
    % start work in all the workers
    parfor worker_ind=1:workers

        if double_scan
            workers_time_points = time_points_to_compute{worker_ind}(1:2:end);
        else
            workers_time_points = time_points_to_compute{worker_ind};            
        end
        
        for t_ind=workers_time_points
            
            if double_scan
                process_time_points = t_ind:t_ind+1;
            else
                process_time_points = t_ind;
            end
            
            scan_positions = round(linspace(1,Width,num_scan));

            if double_scan
                xsection = uint8(zeros(num_scan*Depth, 2*Height));
            else
                xsection = uint8(zeros(num_scan*Depth, Height));
            end
            
            scan_index = 1;
            for t_ind_ind = process_time_points
                disp(['image ' num2str(t_ind_ind)])
                
                % Check if the time point exists (time outs in some experiments).
                if exist([file_folder filesep 'image_' num2str(t_ind_ind,'%.04d') '_trans.raw'],'file')==0
                    continue;
                end

                %% Read in raw input file.
                % Preallocate image memory to speed up computations.
                img_vol=zeros(Width,Height,Depth,'uint8');
                rawInputImageHandle=fopen([file_folder filesep 'image_' num2str(t_ind_ind,'%.04d') '_trans.raw']);
                for iInputImage=1:Depth
                    img_vol(:,:,iInputImage)=...
                        fread(rawInputImageHandle,[Width Height],'uint8');
                end
                fclose(rawInputImageHandle);

                % permute image volume to have the dimensions in order [x,y,z]
                img_vol = permute(img_vol, [2,1,3]);

                if double_scan == 0 || scan_order(scan_index)==2 %if single scan or left part
                    scan_counter = 0;
                    for xs = scan_positions
                        origin = scan_counter*Depth+1;
                        xsection(origin:origin+Depth-1,1:Height) = squeeze(img_vol(:,xs,:))';
                        scan_counter = scan_counter + 1;
                    end
                else % if right part
                    scan_counter = 0;
                    for xs = scan_positions
                        origin = scan_counter*Depth+1;
                        xsection(origin:origin+Depth-1,Height+1:end) = squeeze(img_vol(:,xs,:))';
                        scan_counter = scan_counter + 1;
                    end

                end
                
                scan_index = scan_index+1;
            end
            
            if double_scan
                name_idx = (t_ind+1)/2;
            else
                name_idx = t_ind;
            end
            imwrite(xsection, [out_folder filesep 'img' filesep 'img_', num2str(name_idx, '%04d') '.tif'])
            
        end % t_ind    
    end % worker_ind
    
    eR = expReader([out_folder filesep 'img']);
    eR.ffmpeg;
end