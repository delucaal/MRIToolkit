% files2compile = dir('*.cpp');
files2compile = dir('mex_sh_v2.cpp');

for file_id=1:length(files2compile)
   disp(['Building ' files2compile(file_id).name]);
   cmd = ['mex -v LINKEXPORTCPP="" ' files2compile(file_id).name ' -llapack'];
   try
       eval(cmd);
   catch err
       disp('Failed');
       disp(err)
   end
end

return
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
