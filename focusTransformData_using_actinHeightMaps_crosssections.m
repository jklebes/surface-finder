function focusTransformData_using_actinHeightMaps_crosssections(file_folder_actin, filefolder_myosinORnuclei, out_folder_actin, outfolder_myosinORnuclei, xsection_out_folder_actin, xsection_outfolder_myosinORnuclei, time_points, workers, offsets, varargin)
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
%     of the planes [-5:5] on the surface. - now avg of [-1:1]
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
import utilities.expReader;

%% Parse arguments
%TODO parse arguments propoerly
kernel_size =16;
sq_side = 25;
%width_lim = 40;
overlaps = 3;

% start parallel pool
pool = gcp('nocreate');
if isempty(pool)
    parpool(workers);
elseif pool.NumWorkers~=workers
    delete(pool);
    parpool(workers);
end
%all_surfaces=cell(length(time_points),1); %collect in cell array by
%index
parfor t_ind=1:length(time_points) %parfor
    surfaces_=cell(2);
    if ismember(t_ind, time_points)
          
    %% Read in raw input file - actin
    info_filename=[file_folder_actin filesep '_info.txt'];
    rawInputImage=[file_folder_actin filesep 'image_' num2str(t_ind,'%.04d') '_trans.raw'];
    disp(['image ' num2str(t_ind) '  start read']);
    tstart = tic;
    img_vol_actin = read_3D_image(rawInputImage, info_filename);
    
    tend = toc(tstart);
    disp(['tmp_trans_ actin data reading ' num2str(t_ind) ' took ' num2str(tend) ' seconds.'])

    %% process to heightmaps
    heightmap=calculate_heightmap(img_vol_actin, overlaps, sq_side, kernel_size); %from actin channel
    %saved to

    %% save heightmap
    [~,~, Depth] = size(img_vol_actin);
    save_heightmap(heightmap, Depth, t_ind, out_folder_actin)

    %% extract and save actin surfaces
    surfaces_{1}=select_surfaces(img_vol_actin, t_ind, heightmap, offsets, out_folder_actin);
    
    %% read second channel data
    info_filename=[filefolder_myosinORnuclei filesep '_info.txt'];
    rawInputImage=[filefolder_myosinORnuclei filesep 'image_' num2str(t_ind,'%.04d') '_trans.raw'];
    img_vol_myosin = read_3D_image(rawInputImage, info_filename);

    heightmap_myosin=calculate_heightmap(img_vol_myosin, overlaps, sq_side, kernel_size); %from actin channel

    %% extract and save myosin surfaces
    surfaces_{2}=select_surfaces(img_vol_myosin ,  t_ind, heightmap_myosin, offsets, outfolder_myosinORnuclei);
    %fill array to return
    %all_surfaces{t_ind}=surfaces_;
    else 
    disp(["time point " num2str(t_ind) " not found"]);
    %placeholder for 3D data in cell array remains []
    end
    %% do cross sections, while the 3D data for tindex-1, tindex is in memory
    cross_sections(img_vol_actin, xsection_out_folder_actin, t_ind-1);
    cross_sections(img_vol_actin, xsection_outfolder_myosinORnuclei, t_ind-1);
end % parfor
%(so that parofr can slice by index)
%now we have all_surfaces, a cell array of t/2 time point pairs, each
%containing a 2x2 cell array in actin and myosin, pair upper and lower halves,
%each containing a Height x Widths x n_offsets array of images

%reshape all_surfaces before returning
%into a cell array of timepoints x channels
% all_surfaces_t = cell(length(time_points),1);
% first=@(x) x(1,:);
% all_surfaces_t(1:2:end,1) = cellfun(first,all_surfaces,'UniformOutput',false);
% second=@(x) x(2,:);
% all_surfaces_t(2:2:end,1) = cellfun(second,all_surfaces,'UniformOutput',false);
% vertcat(all_surfaces_t{:})
% %each entry all_surfaces_t(t_ind, channel) contains
% % a Height x Width x n_offsets stack of images
% 
end %function

function data=read_3D_image(filename, info_filename)
%% Read in raw input file. 
% Preallocate image memory to speed up computations.

% Read stack sizes from info file.
Width=findMatchingNumber(info_filename,{'Width: ' '%d'},1);
Height=findMatchingNumber(info_filename,{'Height: ' '%d'},1);
Depth=findMatchingNumber(info_filename,{'Depth: ' '%d'},1);

img_vol=zeros(Width,Height,Depth,'uint8');
rawInputImageHandle=fopen(filename);
for iInputImage=1:Depth
    img_vol(:,:,iInputImage)=...
        fread(rawInputImageHandle,[Width Height],'uint8');
end
fclose(rawInputImageHandle);

% permute image volume to have the dimensions in order [x,y,z]
data = permute(img_vol, [2,1,3]);
end

function heightMap=calculate_heightmap(img_vol, overlaps, sq_side, kernel_size)
[Height, Width, ~] = size(img_vol);
% Store the height map of each repetion
heightMap_stack = zeros(Height,Width,overlaps);

for disIdx = 1:overlaps % 3 planes of heightmap_stack are filled using
    % 3 different sq_size x sq_side griddings,
    % offset by sq_size/overlaps from each
    % other

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

                subI = squeeze(img_vol(i:final_i, j:final_j, z));

                % Squared gradient (Eskicioglu95)
                Ix = diff(subI, 1, 2);
                f = Ix.^2;
                f = sum(f(:)); % measure the focus of the plane
                % f = mean2(f); % not really necesary, commented for speed

                flist(z) = f; % measure the focus of the plate
                ilist(z) = sum(subI(:)); % measure the average int. of the plane

            end

            flist = movmean(flist,10); %smooth the vector
            dfdz = movmean([0,diff(flist)],10);
            [~, h] = max(dfdz);
            heightMap_stack(i:final_i, j:final_j,disIdx) = h;
        end
    end
end

disp('calculate final height map');

%% final heightmap
heightMap = round(max(heightMap_stack,[],3));
kernel = ones(kernel_size)./kernel_size^2;
heightMap=round(imfilter(heightMap,kernel,'same'));
end

function save_heightmap(heightMap, Depth, t_ind, out_folder_actin)
        h_map_folder = [out_folder_actin filesep 'h_map'];
        mkdir(h_map_folder);
        heightMap_image = uint8(heightMap*(256/Depth)); % transform heightmap to 8 bit
        imwrite(heightMap_image, [h_map_folder filesep 'img_' num2str(t_ind,'%.04d') '.jpeg']);
end

function surfaces=select_surfaces(img_vol, t_ind, heightMap, offsets, out_folder)
[Height, Width, Depth]=size(img_vol);
surfaces=zeros(Height, Width, size(offsets,2));
offset_index=1;
for offset = offsets
    final_folder = [out_folder filesep 'surface_' num2str(offset,'%02d')];
    mkdir(final_folder);

    % produce the surface. The image is the average of 3 planes in the surface
    z_averaging = -1:1;%-1:1; %-5:5
    vol = zeros(Height, Width, length(z_averaging), 'uint8');
    counter = 1;
    %TODO array operations?
    for z=z_averaging
        for i = 1:Height
            for j = 1:Width
                zlevel = min([Depth, max([1,heightMap(i,j)+offset+z])]);
                vol(i,j,counter) = img_vol(i,j,zlevel);
            end
        end
        counter = counter + 1;
    end
    surface = uint8(max(vol,[],3));
    %             figure;
    %             imshow(surface_00);

    imwrite(surface,  [final_folder filesep 'img_' num2str(t_ind,'%.04d') '.jpeg']);
    surfaces(:,:,offset_index) = surface;
    offset_index=offset_index+1;
end

end

function cross_sections(img_vol, out_folder, t_ind)
opt_arguments = {
    1,... % double_scan
    [1;2],... % order of scans. Left:2, Right:1
    5,... % number of scans
    };

%double_scan = opt_arguments{1};
scan_order = opt_arguments{2};
num_scan = opt_arguments{3};

%detect whether double scan by input format: call array of volumes or
%volume


[Height, Width, Depth]=size(img_vol);

scan_positions = round(linspace(1,Width,num_scan));

xsection = uint8(zeros(num_scan*Depth, Height));

scan_counter = 0;
for xs = scan_positions
    origin = scan_counter*Depth+1;
    xsection(origin:origin+Depth-1,1:Height) = squeeze(img_vol(:,xs,:))';
    scan_counter = scan_counter + 1;
end


name_idx = t_ind;

mkdir(out_folder)
imwrite(xsection, [out_folder filesep 'img_', num2str(name_idx, '%04d') '.jpeg'])
end
