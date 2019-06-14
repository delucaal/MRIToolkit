%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%%%% Alberto De Luca - alberto@isi.uu.nl $%%%%%% Distributed under the terms of CC BY-NC-ND (https://creativecommons.org/licenses) %%%
classdef MRIToolkitLogger
    properties(Constant)
       LoggingOff = 0
       LoggingOn = 1
       LoggingDebug = 2
       MessageNormal = 0
       MessageDebug = 1
    end
    methods(Static)
        % Logging level = 0, 1, 2 (0 = off, 1 = normal, 2 = debug)
        function setupGlobalVars()
            global MRIKaleidoscope;
            MRIKaleidoscope.Logger.level = 1;
        end
        function setLoggingLevel(a)
            global MRIKaleidoscope;
            MRIKaleidoscope.Logger.level = a;
        end
        function out = getLoggingLevel()
            global MRIKaleidoscope;
            out = MRIKaleidoscope.Logger.level;
        end
        % level = 0 for normal log, = 1 for debug
        function Log(argument,level)
            global MRIKaleidoscope;
            if(MRIKaleidoscope.Logger.level > 0 && level < MRIKaleidoscope.Logger.level)
                disp([datestr(now) ' --- ' argument ' ---']);
            end
        end
    end
end