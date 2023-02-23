function focusSharperstTransitionTransformData_dualChannel(file_folder, out_folder, time_points, workers, varargin)
%{
    Function to focus the surface of the embryo using 45 degrees transformed data. 
    This is a modified verision of the square gradient focusing algorithm. 
    Briefly, the volume is divided in columns. For each plane of the columns 
    it is calculated the average intensity of the image followed by the square 
    gradient in x of the image (mean(diff(subI, 1, 2).^2)) to obtain a 
    measure of sharpness. The next step is to find the plane of maximum
    change in sharpness (peak of the 1st derivative).
    
    This procedure is repeated n times, displacing the origin of the columns
    and the final heightmap is computed as the deepest position for each pixel
    among all the repetitions. Finally the 2D image is produced as the maximum 
    of the planes [-5:5] on the surface.

    Method developed to focus myosin in ACTM1 embryos

    INPUTS
    file_folder: folder containing the raw images
    out_folder: folder to contain the transformed images
    time_interval: array specifying the first and last images
    workers: number of cores for the transformation

    OPTIONAL ARGUMENTS
    kernel_size: size of the averaging kernel for the heightmaps (pixels)
    sq_side: side of the side of the square for focusing (pixels)
    width_lim: limit of the transition zone (pixels)
    overlaps: number of displaced repetitions of the algorithm

    Guillermo 2021

%}

    %% Parse arguments
    opt_arguments = {
                     16,... % kernel_size
                     25,... % sq_side
                     40,... % width_lim
                     3 % overlaps
                     };
    
    opt_arguments(1:nargin-4) = varargin;
    kernel_size = opt_arguments{1};
    sq_side = opt_arguments{2};
    width_lim = opt_arguments{3};
    overlaps = opt_arguments{4};

    %%%%%%%%%%%%%%%%%%%%%%%%%

    workers = 1; 
    time_points = 1:4;
    file_folder = 'A:\analysis\expkTest1\test_height\tmp_transform_data';
    rA.currentIdx = 1;
    results_folder=['A:\analysis\expkTest1\test_height' filesep num2str(rA.currentIdx,'%.4d') '_focusing_Sharperst_Transition'];
    mkdir(results_folder);

    out_folder = 'A:\analysis\expkTest1\test_height';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%     % make output folders
%     mkdir([out_folder filesep 'surface_00']);
%     mkdir([out_folder filesep 'h_map']);

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
    for i_ind=1:workers
%         time_points_to_compute{i_ind}=(i_ind+time_interval(1)-1):workers:time_interval(2);
          time_points_to_compute{i_ind}=time_points(i_ind:workers:length(time_points));
    end
    
    % start work in all the workers
    %parfor worker_ind=1:workers
    worker_ind=1:workers;


        %% run surface finding algorithm independetly for each time point
        for t_ind=time_points_to_compute{worker_ind}
            
            disp(['image ' num2str(t_ind)])
            % Check if the time point exists (time outs in some experiments).
            if exist([file_folder filesep 'image_' num2str(t_ind,'%.04d') '_trans.raw'],'file')==0
                continue;
            end
            
            %% Read in raw input file.
            % Preallocate image memory to speed up computations.
            img_vol=zeros(Width,Height,Depth,'uint8');
            rawInputImageHandle=fopen([file_folder filesep 'image_' num2str(t_ind,'%.04d') '_trans.raw']);
            for iInputImage=1:Depth
                img_vol(:,:,iInputImage)=...
                    fread(rawInputImageHandle,[Width Height],'uint8');
            end
            fclose(rawInputImageHandle);
           disp('read data');

            % permute image volume to have the dimensions in order [x,y,z]
            img_vol = permute(img_vol, [2,1,3]);

            % Store the heigh map of each repetition
            heightMap_stack = zeros(Height,Width,overlaps);

            for disIdx = 1:overlaps % repetions
                
                % every repetion dispace the columns. The mean of the maps
                % is the final height map
                displacement = (disIdx-1)*round(sq_side/overlaps)+1;
                if disIdx>1 % avoid 0s in the final matrix
                    heightMap_stack(:,:,disIdx) = heightMap_stack(:,:,disIdx-1); 
                end

                % crop columns in the stack
                for i = displacement:sq_side:size(img_vol,1)
                    
                    for j = displacement:sq_side:size(img_vol,2)

                        final_i = min([i+sq_side-1, size(img_vol,1)]);
                        final_j = min([j+sq_side-1, size(img_vol,2)]);

                        % initilize vectors
                        flist = zeros(1,size(img_vol,3)); %focus
                        ilist = zeros(1,size(img_vol,3)); %intensity
                        
                        % visit each slice in the column
                        for z = 1:size(img_vol,3)

                            subI = squeeze(img_vol(i:final_i, j:final_j, z));       %writes the elements in the pages of a 3D array into a single array

                            % Squared gradient (Eskicioglu95)
                            Ix = diff(subI, 1, 2);
                            f = Ix.^2;
                            f = sum(f(:)); % measure the focus of the plane
                            % f = mean2(f); % not really necesary, commented for speed

                            flist(z) = f; % measure the focus of the plane                    
                            ilist(z) = sum(subI(:)); % measure the average int. of the plane
                            
                        end

                        flist = movmean(flist,10); %smooth the vector
                        dfdz = movmean([0,diff(flist)],10);
                        [~, h] = max(dfdz);
                        heightMap_stack(i:final_i, j:final_j,disIdx) = h;
                    end
                end
            end

           disp('calculated height map');


            %% final heightmap
            heightMap = round(max(heightMap_stack,[],3));
            kernel = ones(kernel_size)./kernel_size^2;
            heightMap=round(imfilter(heightMap,kernel,'same'));

            %save(['A:\analysis\expkTest1\test_height' filesep 'surface4.mat'],'heightMap')
            save_Parfor_dualChannel(heightMap);
            
            for offset = 0%[-5, 0, 5];%[0, 40, 60, 80]
                final_folder = [out_folder filesep 'surface_' num2str(offset,'%04d')];
                mkdir(final_folder);
                mkdir([final_folder filesep 'surface_00'])

                
                %% produce the surface. The image is the average of 3 planes in the surface
                z_averaging = -5:2;%-1:1; %-5:5
                vol = zeros(size(heightMap,1), size(heightMap,2), length(z_averaging), 'uint8');
                counter = 1;
                for z=z_averaging
                    disp(z)
                    for i = 1:Height
                        for j = 1:Width
                            zlevel = min([Depth, max([1,heightMap(i,j)+offset+z])]);
                            vol(i,j,counter) = img_vol(i,j,zlevel);
                        end
                    end
                    counter = counter + 1;
                end
                surface_00 = uint8(max(vol,[],3));

                %%%%%%%%%%%%%%%%%%%%%%%%
                %output_path = [out_folder filesep 'height_map.mat'];
               % height_data = zeros(Height,Width,time_points)
                %save_Parfor_dualChannel('surface_00')         
                %surface_data = zeros(Height,Depth);
                %surface_data(t_ind) = surface_00;

                  %%%%%%%%%%%%%%%%%%%%%%%%  Shirooza

    %             figure;
    %             imshow(surface_00);

                %% save surface and heightmap
                heightMap_image = uint8(heightMap*(256/Depth)); % transform heightmap to 8 bit

                imwrite(surface_00,  [final_folder filesep 'surface_00' filesep 'img_' num2str(t_ind,'%.04d') '.jpeg']);
                if offset == 0
                    mkdir([final_folder filesep 'h_map'])
                    imwrite(heightMap_image, [final_folder filesep 'h_map' filesep 'img_' num2str(t_ind,'%.04d') '.jpeg']);
                end
            end
            
        end % t_ind    
        
   % end % worker_ind
   
end