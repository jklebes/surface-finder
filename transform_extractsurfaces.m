function transform_extractsurfaces(transform, file_folder, out_folder,  xsection_out_folder, expFolder, info_file, matrix_path, time_points, workers, offsets, varargin)
%{
%     Three-in-one: optionally 45deg transform, extract focussed surfaces,
%     and extract cross sections
%
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
%     file_folder: folder containing the raw images
%     
%     out_folder: folder to contain thesurface images
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
%TODO parse arguments propoerly
kernel_size =16;
sq_side = 25;
%width_lim = 40;
overlaps = 3;


if transform
    %% prepare affine transformation of raw data
    %read info files
    width = findMatchingNumber(info_file,{'Org Width:','%d'});
    height = findMatchingNumber(info_file,{'Org Height:','%d'});
    num_slices = findMatchingNumber(info_file,{'No Steps:','%d'});
    original_dims = [width; height; num_slices;0];
    %retrieve matrix
    matrix = load(matrix_path);
    %calculate new dimensions
    new_dims = round(matrix*original_dims); %matrix multiply
    %write new info file
    transf_info_file = [out_folder filesep '_info.txt'];
    fileID = fopen(transf_info_file, 'w');
    fprintf(fileID,['Width: ', num2str(new_dims(1)) '\n']); %it is unclear how imageJ does the rounding.
    fprintf(fileID,['Height: ',  num2str(new_dims(2)) '\n']);
    fprintf(fileID,['Depth: ',  num2str(new_dims(3)) '\n']);
    fclose(fileID);
    chunk_size=3*workers;
else
    info_file = [file_folder filesep '_info.txt']; %use post-transform dims
    %restarting pool not needed - everything in one "chunk"
    chunk_size=length(time_points);
    %these need to be defined going into parfor loop, even if not used
    matrix=[];
    original_dims=[];
    new_dims=[];
    
end

% Compute the number of images and the format
format=[];
if isempty(format)
    files = dir(expFolder);
    files = {files.name};
    files = {files{cellfun(@(x) length(x)>2, files)}};
    if isempty(files)
        format = [];
    else
        mainLength = cellfun(@(x) length(x),files);
        [repet,categ]=hist(mainLength,unique(mainLength));
        [~, maxRep] = max(repet);
        mainLength = mainLength == categ(maxRep);
        files = files(mainLength);

        constant = ones(size(files{1}));
        for f = 2:length(files)
            constant = constant.*(files{f-1}==files{f});
        end
        numlength = sum(~constant); % number of zeros

        format{1} = files{1}(1:find(~constant,1,'first')-1);
        format{2} = ['%' num2str(numlength,'%02d') 'd'];
        format{3} = files{1}(find(~constant,1,'last')+1:length(constant));
    end
end

% start parallel pool
pool = gcp('nocreate');
if isempty(pool)
    parpool(workers);
elseif pool.NumWorkers~=workers
    delete(pool);
    parpool(workers);
end

%prepare output directories 
h_map_folder = [out_folder filesep 'h_map'];
mkdir(h_map_folder);
final_folder=cell(size(offsets));
offset_index=1;
for offset=offsets
    final_folder{offset_index} = [out_folder filesep 'surface_' num2str(offset,'%02d')];
    mkdir(final_folder{offset_index});
    offset_index=offset_index+1;
end
mkdir(xsection_out_folder)
%calculate dividing timepoints into chunks
remainder_length=mod(length(time_points), chunk_size);
if length(time_points)-remainder_length>0
chunks=reshape(time_points(1:end-remainder_length), chunk_size,[]);
else
    chunks=[];
end
chunk_indices = cell([1 size(chunks,1)+1]);
for i=1:size(chunks,2)
chunk_indices{i} = chunks(:,i);
end
if remainder_length>0
chunk_indices{end} = time_points(end-remainder_length+1:end);
end
for chunk=1:size(chunks,2)+1 %chunk
    %the chunking is because of inadequate activity of java garbage collector: 
    %despite deallocating, something builds up in the
    %memory of each worker and causes 'Out of Memory' error farther into the
    %loop.  Workaround is to restart the pool periodically.  Takes extra
    %time.  Kicks in at transform=true route only because of chunk_size 
    % chosen above.
parfor t=1:length(chunk_indices{chunk}) %parfor
        t_ind=chunk_indices{chunk}(t);
        if ~transform
            %read in already-transformed 3D data
            rawInputImage=[file_folder filesep 'image_' num2str(t_ind,'%.04d') '_trans.raw'];
            disp(['image ' num2str(t_ind) '  start read']);
            tstart = tic;
            img_vol = read_3D_image(rawInputImage, info_file);
            tend = toc(tstart);
            disp(['tmp_trans_ data reading ' num2str(t_ind) ' took ' num2str(tend) ' seconds.'])
        else
            img_path = [expFolder filesep format{1} num2str(t_ind, format{2}) format{3}];
            %% call program that uses java to read in, affine transform raw 3d data
            java.lang.System.gc();
            disp(['image ' num2str(t_ind) ', start read and affine transform']);
            tstart = tic;
            img_vol =transform_45degrees_Java_single(img_path, matrix, original_dims(1:3), new_dims(1:3));   
            tend = toc(tstart);
            java.lang.System.gc(); %request java garbage collection, mostly works
            disp(['data read and transform ' num2str(t_ind) ' took ' num2str(tend) ' seconds.'])
        end

        %% extract heightmaps with surface finding algorithm
        heightmap=calculate_heightmap(img_vol, overlaps, sq_side, kernel_size); 

        %% save heightmap
        [~,~, Depth] = size(img_vol);
        save_heightmap(heightmap, Depth, t_ind, h_map_folder)

        %% extract and save surfaces
        select_surfaces(img_vol, t_ind, heightmap, offsets, final_folder); 

        %% do cross sections, while the 3D data is in memory
        cross_sections(img_vol, xsection_out_folder, t_ind);
end % parfor 
java.lang.System.gc();
%restarting the pool because nothing else works to clear memory
%accumulating in each worker,
%but this is very bad and slow.
pool = gcp('nocreate');
delete(pool);
parpool(workers);
end %chunks
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
    % 3 different sq_side x sq_side griddings,
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

function save_heightmap(heightMap, Depth, t_ind, h_map_folder)
heightMap_image = uint8(heightMap*(256/Depth)); % transform heightmap to 8 bit
my_imwrite(heightMap_image, [h_map_folder filesep 'img_' num2str(t_ind,'%.04d') '.jp2']);
end

function select_surfaces(img_vol, t_ind, heightMap, offsets, final_folder)
[Height, Width, Depth]=size(img_vol);
surfaces=zeros(Height, Width, size(offsets,2));
offset_index=1;
for offset = offsets

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

    my_imwrite(surface,  [final_folder{offset_index} filesep 'img_' num2str(t_ind,'%.04d') '.jp2']);
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
my_imwrite(xsection, [out_folder filesep 'img_', num2str(name_idx, '%04d') '.jp2'])
end
