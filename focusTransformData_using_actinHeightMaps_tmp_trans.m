function focusTransformData_using_actinHeightMaps_tmp_trans(rot_vol_actin, rot_vol_myosinORnuclei, Width, Height, Depth, results_folder_actin, results_folder_myosinORnuclei, t_ind, varargin)
%{
%     Function to focus the surface of the embryo using 45 degrees transformed data. 
%     This is a modified verision of the square gradient focusing algorithm. 
%     Briefly, the volume is divided in columns. For each plane of the columns 
%     it is calculated the average intensity of the image followed by the square 
%     gradient in x of the image (mean(diff(subI, 1, 2).^2)) to obtain a 
%     measure of sharpness. The next step is to find the plane of maximum
%     change in sharpness (peak of the 1st derivative).
%     
%     This procedure is repeated n times, displacing the origin of the columns
%     and the final heightmap is computed as the deepest position for each pixel
%     among all the repetitions. Finally the 2D image is produced as the maximum 
%     of the planes [-5:5] on the surface.
% 
%     Method developed to focus myosin in ACTM1 embryos
% 
%     INPUTS
%     file_folder_actin: folder containing the raw actin images
%     file_folder_myosinORnuclei: folder containing the raw myosin images
%     
%     out_folder_actin: folder to contain the actin surface images
%     out_folder_myosinORnuclei: folder to contain the myosin surface images
% 
%     time_interval: array specifying the first and last images
%     workers: number of cores for the transformation
% 
%     OPTIONAL ARGUMENTS
%     kernel_size: size of the averaging kernel for the heightmaps (pixels)
%     sq_side: side of the side of the square for focusing (pixels)
%     width_lim: limit of the transition zone (pixels)
%     overlaps: number of displaced repetitions of the algorithm
% 
%     Guillermo 2021

%}


    %% Parse arguments
%   opt_arguments = {
%                      16,... % kernel_size
%                      25,... % sq_side
%                      40,... % width_lim
%                      3 % overlaps
%                      };
%     
%     opt_arguments(1:nargin-4) = varargin;
%     kernel_size = opt_arguments{1};
%     sq_side = opt_arguments{2};
%     width_lim = opt_arguments{3};
%     overlaps = opt_arguments{4};
    
     kernel_size =16;
     sq_side = 25;
     width_lim = 40;
     overlaps = 3;

%     % make output folders
%     mkdir([out_folder_actin filesep 'surface_00']);
%     mkdir([out_folder_actin filesep 'h_map']);

%     % Read stack sizes from info file.
%      Width=findMatchingNumber([containerFolder filesep 'rotating_info.txt'],{'Width: ' '%d'},1);
%      Height=findMatchingNumber([containerFolder filesep 'rotating_info.txt'],{'Height: ' '%d'},1);
%      Depth=findMatchingNumber([containerFolder filesep 'rotating_info.txt'],{'Depth: ' '%d'},1);
%     
 
%  %%%%%%%% not using parallel pool for the current script
%      %% start parallel pool
%      pool = gcp('nocreate');
%      if isempty(pool)
%         parpool(workers);
%      elseif pool.NumWorkers~=workers
%         delete(pool);
%         parpool(workers);
%      end
%  
%      % split time points for the parallel workers for the images to be produced in order
%      time_points_to_compute=cell(workers,1);
%      for i_ind=1:workers
% %        time_points_to_compute{i_ind}=(i_ind+time_interval(1)-1):workers:time_interval(2);
%            time_points_to_compute{i_ind}=time_points(i_ind:workers:length(time_points));
%      end
% %     
% %     % start work in all the workers
%      parfor worker_ind=1:workers
%         for t_ind=time_points_to_compute{worker_ind}
% 
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        %for t_ind=1:time_points            

%  %% Read in raw input file.
%             % Preallocate image memory to speed up computations.
%             disp('image 1 start read');
%     tic
%   
%             img_vol=zeros(Width,Height,Depth,'uint8');
%             rawInputImageHandle=fopen([file_folder_actin filesep 'image_' num2str(t_ind,'%.04d') '_trans.raw']);
%             for iInputImage=1:Depth
%                 img_vol(:,:,iInputImage)=...
%                     fread(rawInputImageHandle,[Width Height],'uint8');
%             end
%             fclose(rawInputImageHandle);
%             
%             disp('image 1 read finished');
%    toc
            % permute image volume to have the dimensions in order [x,y,z]
twhole = tic; 
            
            tstart = tic; 
            rot_vol_actin = permute(rot_vol_actin, [2,1,3]);
            tend = toc(tstart);
             disp([size(rot_vol_actin)])
            disp(['permuting took ' num2str(tend) ' seconds.'])  %% 1.3 s

            % Store the heigh map of each repetion
            heightMap_stack = zeros(Height,Width,overlaps);

            tstart = tic; 
            for disIdx = 1:overlaps % repetions
                
                % every repetion dispace the columns. The mean of the maps
                % is the final height map
                displacement = (disIdx-1)*round(sq_side/overlaps)+1;
                if disIdx>1 % avoid 0s in the final matrix
                    heightMap_stack(:,:,disIdx) = heightMap_stack(:,:,disIdx-1); 
                end

                disp('crop columns in the stack')

                % crop columns in the stack
                 ti = tic; 
                for i = displacement:sq_side:size(rot_vol_actin,1)
                    
                    tj = tic; 
                    for j = displacement:sq_side:size(rot_vol_actin,2)

                        final_i = min([i+sq_side-1, size(rot_vol_actin,1)]);
                        final_j = min([j+sq_side-1, size(rot_vol_actin,2)]);

                        % initilize vectors
                        flist = zeros(1,size(rot_vol_actin,3)); %focus
%                         ilist = zeros(1,size(rot_vol_actin,3)); %intensity
                        
                        % visit each slice in the column
                         tz = tic; 
                        for z = 1:size(rot_vol_actin,3)

                            subI = squeeze(rot_vol_actin(i:final_i, j:final_j, z));

                            % Squared gradient (Eskicioglu95)
                            Ix = diff(subI, 1, 2);
                            f = Ix.^2;
                            f = sum(f(:)); % measure the focus of the plane
                            % f = mean2(f); % not really necesary, commented for speed

                            flist(z) = f; % measure the focus of the plate                    
%                             ilist(z) = sum(subI(:)); % measure the average int. of the plane
                            %%commented the ilist as it's not being used
                            
                        end
                        tendz = toc(tz);
                        disp(['z loop for the depth(283) took ' num2str(tendz) ' seconds.']) %% 0.005-0.01s -- time for 283 z loops

                        flist = movmean(flist,10); %smooth the vector
                        dfdz = movmean([0,diff(flist)],10);
                        [~, h] = max(dfdz);
                        heightMap_stack(i:final_i, j:final_j,disIdx) = h;
                    end
                    tendj = toc(tj);
                    disp(['***j loop took ' num2str(tendj) ' seconds.    j value' num2str(j) ]) %% 1.8-3s -- (time for ~269 j loops)
                   
                end
                disp('end of i loop')
                 tendi = toc(ti);
                 disp(['----i loop took ' num2str(tendi) ' seconds.---- i value' num2str(i)]) %% 242s (4min)  -- (time for ~103 i loops)
            end

           tend = toc(tstart);
           disp(['crop columns in the stack for all the 3 overlaps took ' num2str(tend) ' seconds.']) %% 690s (11.5 min)  -- time for 3 overlaps
            
    disp('calculate final height map');
            
        %% final heightmap
            heightMap = round(max(heightMap_stack,[],3));
            kernel = ones(kernel_size)./kernel_size^2;
            heightMap=round(imfilter(heightMap,kernel,'same'));         

            toffset = tic; 
            for offset = [0, 30] %[-5, 0, 5];%[0, 40, 60, 80]
                final_folder = [results_folder_actin filesep 'surface_' num2str(offset,'%04d')];
                mkdir(final_folder);
                mkdir([final_folder filesep 'surface_00'])
                
                %% produce the surface. The image is the average of 3 planes in the surface
                z_averaging = -1:1;%-1:1; %-5:5
                vol = zeros(size(heightMap,1), size(heightMap,2), length(z_averaging), 'uint8');
                counter = 1;

                tstart = tic; 
                for z=z_averaging
                    for i = 1:Height
                        for j = 1:Width
                            zlevel = min([Depth, max([1,heightMap(i,j)+offset+z])]);
                            vol(i,j,counter) = rot_vol_actin(i,j,zlevel);
                        end
                    end
                    counter = counter + 1;
                    disp('end of z loop')
                end
                surface_00 = uint8(max(vol,[],3));

                tend = toc(tstart);
                disp(['producing the surface took ' num2str(tend) ' seconds.   for offset value: ' num2str(offset)]) %% 8.9 s
            
    %             figure;
    %             imshow(surface_00);

                %% save surface and heightmap
%      disp('save surface 1');
                tstart = tic;
                heightMap_image = uint8(heightMap*(256/Depth)); % transform heightmap to 8 bit

                imwrite(surface_00,  [final_folder filesep 'surface_00' filesep 'img_' num2str(t_ind,'%.04d') '.jpeg']);
                disp('writing')
                if offset == 0
                    mkdir([final_folder filesep 'h_map'])
                    imwrite(heightMap_image, [final_folder filesep 'h_map' filesep 'img_' num2str(t_ind,'%.04d') '.jpeg']);
                end
                tend = toc(tstart);
                disp(['writing surface and height map took ' num2str(tend) ' seconds.']) %% 0.7 s
            
            end  %% for each offset value it takes nearly 10s for this loop
            tendA = toc(toffset);
            disp(['to produce a single surface actin ' num2str(tendA) ' seconds.']) %% 10.11 s  
            %%% if there are two offset values the above  time is doubled (10x2 = 20s))
            
            disp(['image ' num2str(t_ind)])


%%%%%%%%%%%%%%% does not have to read the images
%             % Check if the time point exists (time outs in some experiments).
%             if exist([filefolder_myosinORnuclei filesep 'image_' num2str(t_ind,'%.04d') '_trans.raw'],'file')==0
%                 continue;
%             end
                        
            
%             %% Read in raw input file.
%             % Preallocate image memory to speed up computations.

%             img_vol=zeros(Width,Height,Depth,'uint8');
%             rawInputImageHandle=fopen([filefolder_myosinORnuclei filesep 'image_' num2str(t_ind,'%.04d') '_trans.raw']);
%             for iInputImage=1:Depth
%                 img_vol(:,:,iInputImage)=...
%                     fread(rawInputImageHandle,[Width Height],'uint8');
%             end
%             fclose(rawInputImageHandle);
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % permute image volume to have the dimensions in order [x,y,z]
            rot_vol_myosinORnuclei = permute(rot_vol_myosinORnuclei, [2,1,3]);
            toffsetM = tic; 
            for offset = [0, 30]    %[-5, 0, 5];%[0, 20, 40, 60, 80]
                final_folder = [results_folder_myosinORnuclei filesep 'surface_' num2str(offset,'%04d')];
                mkdir(final_folder);
                mkdir([final_folder filesep 'surface_00'])
                
                %% produce the surface. The image is the average of 3 planes in the surface
                z_averaging = -1:1;%-1:1; %-5:5
                vol = zeros(size(heightMap,1), size(heightMap,2), length(z_averaging), 'uint8');
                counter = 1;
                for z=z_averaging
                    for i = 1:Height
                        for j = 1:Width
                            zlevel = min([Depth, max([1,heightMap(i,j)+offset+z])]);
                            vol(i,j,counter) = rot_vol_myosinORnuclei(i,j,zlevel);
                        end
                    end
                    counter = counter + 1;
                end
                surface_00 = uint8(max(vol,[],3));
    %             figure;
    %             imshow(surface_00);

                %% save surface %%   and heightmap

                
                heightMap_image = uint8(heightMap*(256/Depth)); % transform heightmap to 8 bit

                imwrite(surface_00,  [final_folder filesep 'surface_00' filesep 'img_' num2str(t_ind,'%.04d') '.jpeg']);


%                 if offset == 0
%                     mkdir([final_folder filesep 'h_map'])
%                     imwrite(heightMap_image, [final_folder filesep 'h_map' filesep 'img_' num2str(t_ind,'%.04d') '.jpeg']);
%                 end
            end
            tendM = toc(toffsetM);
            disp(['to produce a single surface myosin ' num2str(tendM) ' seconds.  for offset value: ' num2str(offset)]) %% 9.94 s
             %%% Now again if there are two offset values the above  time is doubled (10x2 = 20s))           
                     
%                     end % t_ind    
%    end % worker_ind

            tendW = toc(twhole);
            disp(['whole surface calculation for both channels took ' num2str(tendW) ' seconds.']) %% 715 s (12 mins)
            %%% for two off set values, 715+ 20 = 735s 
end