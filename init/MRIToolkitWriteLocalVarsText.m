function MRIToolkitWriteLocalVarsText(text)
    floc = which('MRIToolkitDefineLocalVars.m');
    if(isempty(floc))
        floc = fullfile(userpath,'MRIToolkitDefineLocalVars.m');
    end
    f = fopen(floc,'wt');
    for ij=1:length(text)
        fprintf(f,'%s%s',text{ij},newline);
    end
    fclose(f);
end