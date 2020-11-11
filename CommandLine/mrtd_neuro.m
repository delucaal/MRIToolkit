function mrtd_neuro(varargin)
    disp('mrtd_neuro');
    coptions = varargin;
    if(length(varargin{1}) > 1)
        coptions = varargin{1};
    end
    %     disp(varargin)

    if(isempty(coptions) || isempty(coptions{1}) || strcmpi(coptions{1},'-help'))
        help = 'This tool applies some popular neuroimaging pipelines';
        help = [help newline];
        help = [help newline 'FOD-based usage: mrtd_neuro -t1 t1_file.nii -flair flair.nii -apply cat12pipeline -out output_folder'];
        help = [help newline];
        help = [help newline '-t1: T1-weighted .nii file'];
        help = [help newline '-flair: FLAIR .nii file'];
        help = [help newline '-out: output folder'];
        help = [help newline '-apply: "cat12pipeline"'];
        help = [help newline];
        fprintf(help);

        return
    end

t1_file = '';
flair_file = '';
outfile = '';
method = '';

for input_id = 1:2:length(coptions)
    value = coptions{input_id};
    if(strcmp(value,'-t1'))
        t1_file = coptions{input_id+1};
    elseif(strcmp(value,'-flair'))
        flair_file = coptions{input_id+1};
    elseif(strcmp(value,'-out'))
        outfile = coptions{input_id+1};
    elseif(strcmp(value,'-apply'))
        method = (coptions{input_id+1});
    end
    
end

if(isempty(outfile))
    error('Missing mandatory argument -out');
end

try
%     cat_defaults
%     cat12
%     close all
catch err
    disp(err.message);
end

if(strcmpi(method,'cat12pipeline'))
    if(isempty(t1_file))
        error('Missing mandatory argument -t1');
    end
    should_remove_t1 = 0;
    if(contains(t1_file,'nii.gz'))
        t1_file = gunzip(t1_file);
        should_remove_t1 = 1;
    end
    should_remove_flair = 0;
    if(contains(flair_file,'nii.gz'))
        flair_file = gunzip(flair_file);      
        should_remove_flair = 1;
    end
    t1_file = fullfile(pwd,t1_file);
    flair_file = fullfile(pwd,flair_file);
    Neuro.CAT12Pipeline('nii_file',t1_file,'flair_file',flair_file,...
        'output',outfile);
    if(should_remove_t1 == 1)
        delete(t1_file);
    end
    if(should_remove_flair == 1)
        delete(flair_file);
    end
end
end
