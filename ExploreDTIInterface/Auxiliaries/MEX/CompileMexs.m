files2compile = dir('*.cpp');

for file_id=1:length(files2compile)
   disp(['Building ' files2compile(file_id).name]);
   cmd = ['mex ' files2compile(file_id).name ' -llapack'];
   try
       eval(cmd);
   catch
       disp('Failed');
   end
end

files2compile = dir('*.c');

for file_id=1:length(files2compile)
   disp(['Building ' files2compile(file_id).name]);
   cmd = ['mex ' files2compile(file_id).name ' -llapack'];
   try
       eval(cmd);
   catch
       disp('Failed');
   end
end
