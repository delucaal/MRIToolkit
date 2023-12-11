function text = MRIToolkitReadLocalVarsText()
    text = {};
    floc = which('MRIToolkitDefineLocalVars.m');
    if(isempty(floc))
        floc = which('TemplateMRIToolkitDefineLocalVars.m');
    end
    if(isempty(floc))
        warning('I cannot find a template initialization file. Something may be wrong with this installation.')
        return
    end
    f = fopen(floc,'rt');
    while(~feof(f))
        l = fgetl(f);
        text = cat(1,text,{l});
    end
    fclose(f);
end