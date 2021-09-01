function v = chk_dir(folder_name)
%  checks if folder_name exists. if not makes one.
v = exist(folder_name, 'dir'); 
if v == 0
    mkdir(folder_name);
end


