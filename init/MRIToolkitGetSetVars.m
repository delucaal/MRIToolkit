function varargout = MRIToolkitGetSetVars(varargin)
    text = MRIToolkitReadLocalVarsText();
    if(nargout == 1)
        % expecting output -> get
        if(nargin == 1)
            % input 1 = variable
            found = -1;
            for ij = 1:length(text)
                if(contains(text{ij},varargin{1}))
                    % parameter found
                    found = ij;
                    break;
                end
            end
            if(found == -1)
                varargout{1} = '';
                warning('Text not found');
            else
                l = strfind(text{found},'''');
                varargout{1} = text{found}(l(1)+1:l(end)-1);                
            end
        else
            error('Invalid input/output combination (1 input - output expected)');
        end
    elseif(nargout == 0)
        % Set mode
        if(nargin == 2)
            % input 1 = variable, input 2 = new value
            found = -1;
            for ij = 1:length(text)
                if(contains(text{ij},varargin{1}))
                    % parameter found
                    found = ij;
                    break;
                end
            end
            if(found == -1)
                for ij=length(text):-1:1
                    if(contains(text{ij},'end')) % Insert a new statement before the end
                        break
                    end
                end
                new_text = text(1:ij-1);
                new_text(ij) = {[varargin{1} ' = ''' varargin{2} ''';']};
                new_text = cat(1,new_text,text(ij:end));
                text = new_text;
            else
                l = strfind(text{found},'''');
                text{found} = strrep(text{found},text{found}(l(1)+1:l(end)-1),varargin{2});                
            end
            MRIToolkitWriteLocalVarsText(text);
        end
    end
end