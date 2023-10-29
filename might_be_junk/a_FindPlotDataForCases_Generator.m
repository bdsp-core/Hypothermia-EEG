folder = '/Users/mbw/cdac Dropbox/Datasets/zz_Hypothermia/CodeAndData';

mfiles = dir(fullfile(folder,'*.m'));
mfilenames = {mfiles.name};

for i = 1:length(mfilenames)
  
  filename = fullfile(folder, mfilenames{i});
  
  fileContents = fileread(filename);
  
  if contains(fileContents, 'cs1_bursts')
    disp([mfilenames{i} ' generates the output files']);
  end
  
end