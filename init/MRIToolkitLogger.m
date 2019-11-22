%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%%%% Alberto De Luca - alberto@isi.uu.nl $%%%%%% Distributed under the terms of LGPLv3  %%%



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
            global MRIToolkit;
            MRIToolkit.Logger.level = 1;
        end
        function setLoggingLevel(a)
            global MRIToolkit;
            MRIToolkit.Logger.level = a;
        end
        function out = getLoggingLevel()
            global MRIToolkit;
            out = MRIToolkit.Logger.level;
        end
        % level = 0 for normal log, = 1 for debug
        function Log(argument,level)
            global MRIToolkit;
            if(MRIToolkit.Logger.level > 0 && level < MRIToolkit.Logger.level)
                disp([datestr(now) ' --- ' argument ' ---']);
            end
        end
    end
end