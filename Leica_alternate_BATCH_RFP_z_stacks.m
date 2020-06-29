filename_originals = uigetfile('*RFP_t1.TIF','Multiselect','on');
mkdir('Tiff Files') 
filename_originals = filename_originals';
num_timelapses = length(filename_originals);
%% 

for j=1:num_timelapses %
    filename_replace = strrep(filename_originals {j}, 't1', '*');
    a = which(filename_originals{j}); %give it an initial file to start with 
    filelist = dir([fileparts(a) filesep filename_replace]);
    fileNames = {filelist.name}';
    filenames_sorted = natsortfiles(fileNames);
    number_of_files=length(fileNames); 
    num_files_expected = 30;
    if number_of_files <num_files_expected
         %msgbox('Timelapse incomplete you fool!','Error')
        continue
    end
    %mkdir('Tiff Files') 
    %% 
    filename = sprintf('Tiff Files/ZEY173_Z_Halo_%02d.tif', j-1) ;
    %colour_code = 3;
    %colours =3 ;
    %rd = bfopen('MCM4_Halo_.nd');
    %num_files = length(rd{2,1});
    for i = 1:number_of_files
        img_name = filenames_sorted{i};
        info_img = imfinfo(img_name);
        num_images = numel(info_img);
        IMG = [];
        for k = 1:num_images
            
            image_file = imread(img_name,k);
            IMG(:,:,k) = image_file;
            
        end
    max_Z_IMG = max(IMG,[],3);
    max_Z_16 = uint16(max_Z_IMG);
    imwrite(max_Z_16,filename,'tif','WriteMode','append');
    end
    disp('Success')
    disp(filename_originals{j})
end
