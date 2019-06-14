%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%%%% Alberto De Luca - alberto@isi.uu.nl $%%%%%% Distributed under the terms of CC BY-NC-ND (https://creativecommons.org/licenses) %%%
% MRIKaleidoscope RelaxometryMethods class - A. De Luca
% T1 and T2 quantification methods
% History: 
% v1: 23/09/2018
% 07/06/2019: Added absolute value in IR_T1_signal
%             Fixed T2->R2
classdef RelaxometryMethods < handle
    
    methods (Static)
    
        % The classic inversion recovery equation for T1 quantification
        % given multiple inversion times
        function [out,prediction] = IR_T1_signal(x0,params)
            Amplitude = x0(1);
            T1 = x0(2);
            prediction = abs(Amplitude*(1-2*exp(-params.TI/T1)+exp(-params.TR./T1)));
            out = params.S-prediction;
        end
    
        % Quantify T1 of a 1D normalized signal using a grid search method
        % Input:
        % T1_min_max = [min,max] in seconds;
        % T1_points = number of points to discretize the linear range
        % measured_signal: optimization target
        % params: needs params.TR = Repetition Time, params.TI = array of
        % inversion times
        % Output:
        % qT1: the quantified T1 value
        % cost: resnorm
        % The IR fit
        function [qT1,cost,fit] = IR_T1_fit_1D(T1_min_max,T1_points,measured_signal,params)
            T1Range = linspace(T1_min_max(1),T1_min_max(2),T1_points)';
            parameter_values = [ones(size(T1Range)) T1Range];
            [qT1,cost,fit] = OptimizationMethods.GridSearch1D(@RelaxometryMethods.IR_T1_signal,measured_signal,...
                parameter_values,params);
        end
        
        % Quantify T1 of an array of 1D normalized signal using a grid search method
        % Input:
        % T1_min_max = [min,max] in seconds;
        % T1_points = number of points to discretize the linear range
        % measured_signals: array of optimization target
        % params: needs params.TR = Repetition Time, params.TI = array of
        % inversion times
        % Output:
        % qT1: the quantified T1 value
        % cost: resnorm
        % The IR fit
        function [qT1,cost,fit] = IR_T1_fit_ND(T1_min_max,T1_points,measured_signals,params)
            T1Range = linspace(T1_min_max(1),T1_min_max(2),T1_points)';
            parameter_values = [ones(size(T1Range)) T1Range];
            [qT1,cost,fit] = OptimizationMethods.GridSearchND(@RelaxometryMethods.IR_T1_signal,...
                measured_signals,parameter_values,params);
        end
        
        % Quantify R2 (1/T2) of a signal using a LLS fit
        % Input:
        % measured_signal: optimization target
        % params: needs params.TE = Echo Time
        % Output:
        % qT2: the quantified T2 value (and intercept)
        % cost: resnorm
        % fit: The SE fit
        function [qR2,cost,fit] = SE_R2_fit_1D(measured_signal,params)
            X = [ones(size(params.TE)) -params.TE];
            qR2 = X\log(measured_signal);
            fit = exp(X*qR2);
            cost = sum((measured_signal - fit).^2);
        end
        
        % Quantify R2 (1/T2) of a signal array using a LLS fit
        % Input:
        % measured_signal: optimization target
        % params: needs params.TE = Echo Time
        % Output:
        % qT2: the quantified T2 value (and intercept)
        % cost: resnorm
        % fit: The SE fit
        function [qR2,cost,fit] = SE_R2_fit_ND(measured_signals,params,use_nnls)
            X = [ones(size(params.TE)) -params.TE];
            qR2 = zeros(size(measured_signals,1),2);
            cost = zeros(size(measured_signals,1),1);
            fit = zeros(size(measured_signals));
            if(nargin < 3)
                use_nnls = 0;
            end
            if(use_nnls == 0)
                for x=1:size(measured_signals,1)
                    S = measured_signals(x,:)';
                    if(sum(S == 0) > 0)
                        continue
                    end
                    qR2(x,:) = X\log(S);
                    fit(x,:) = exp(X*qR2(x,:)');
                    cost(x) = sum((S - fit(x,:)').^2);
                end
            else
                for x=1:size(measured_signals,1)
                    S = double(measured_signals(x,:)');
                    if(sum(S == 0) > 0)
                        continue
                    end
                    qR2(x,:) = lsqnonneg(X,log(S));
                    fit(x,:) = exp(X*qR2(x,:)');
                    cost(x) = sum((S - fit(x,:)').^2);
                end
            end
        end
        
    end  
end