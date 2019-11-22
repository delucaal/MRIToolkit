%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%%%% Alberto De Luca - alberto@isi.uu.nl $%%%%%% Distributed under the terms of LGPLv3  %%%



% MRIToolkit OptimizationMethods class - A. De Luca
% Optimization methods, as grid search
% History:
% v1: 23/09/2018
% 07/06/2019: Least-Squares optimization introduced in the grid-search

classdef OptimizationMethods < handle
    
    methods (Static)
        
        % Input:
        % function_handle: function(p,params) returning the signal vector
        % S: the target signal
        % parameter_values: possible solutions
        % params: structure passed to function_handle
        % Output:
        % best_p: parameter providing the best fit
        % cost: corresponding resnorm given best_p
        % fit: the best fit (corresponding to best_p)
        function [best_p,cost,fit] = GridSearch1D(function_handle,S,parameter_values,params)
            COST = zeros(size(parameter_values,1),1);
            params.S = S;
            for ij=1:size(parameter_values,1)
                [~,prediction] = function_handle(parameter_values(ij,:),params);
                X = [prediction ones(size(prediction))];
                p = X\S;
                prediction = X*p;
                COST(ij) = sum((S-prediction).^2);
            end
            if(any(~isfinite(COST)))
                error('Unexpected values in GridSearch1D');
            end
            
            [cost,IX] = min(COST);
            best_p = parameter_values(IX,:);
            [~,fit] = function_handle(best_p,params);
            X = [fit ones(size(fit))];
            p = X\S;
            fit = X*p;
        end
        
        % Input:
        % function_handle: function(p,params) returning the signal vector
        % S: the target signal
        % NS: Signal of the neighbours
        % parameter_values: possible solutions
        % params: structure passed to function_handle
        % reg_strength: weight for neighbouring voxels
        % Output:
        % best_p: parameter providing the best fit
        % cost: corresponding resnorm given best_p
        % fit: the best fit (corresponding to best_p)
        function [best_p,cost,fit] = GridSearch1DSpatReg(function_handle,S,NS,parameter_values,params,reg_strength)
            COST = zeros(size(parameter_values,1),1);
            params.S = S;
            for ij=1:size(parameter_values,1)
                [~,prediction] = function_handle(parameter_values(ij,:),params);
                X = [prediction ones(size(prediction))];
                p = X\S;
                prediction = X*p;
                COST(ij) = sum((S-prediction).^2);
            end
            
            for ns=1:size(NS,2)
                params.S = NS(:,ns);
                for ij=1:size(parameter_values,1)
                    [~,prediction] = function_handle(parameter_values(ij,:),params);
                    X = [prediction ones(size(prediction))];
                    p = X\S;
                    prediction = X*p;
                    COST(ij) = COST(ij)+reg_strength*sum((S-prediction).^2);
                end           
            end
            
            [cost,IX] = min(COST);
            best_p = parameter_values(IX,:);
            [~,fit] = function_handle(best_p,params);
            X = [fit ones(size(fit))];
            p = X\S;
            fit = X*p;
        end        
        
        % Input:
        % function_handle: function(p,params) returning the signal vector
        % S: the target signal array
        % parameter_values: possible solutions
        % params: structure passed to function_handle
        % Output:
        % best_p: parameter providing the best fit
        % cost: corresponding resnorm given best_p
        % fit: the best fit (corresponding to best_p)
        function [best_p,cost,fit] = GridSearchND(function_handle,S,parameter_values,params)
            Dictionary = zeros(size(S,2),size(parameter_values,1));
            params.S = S(1,:)';
            for ij=1:size(parameter_values,1)
                [~,Dictionary(:,ij)] = function_handle(parameter_values(ij,:),params);
            end
            best_p = zeros(size(S,1),size(parameter_values,2));
            cost = zeros(size(S,1),1);
            fit = zeros(size(S));
            
            Dictionary = Dictionary';
            
            for x=1:size(S,1)
                COST = zeros(size(parameter_values,1),1);
                K = S(x,:)';
                if(sum(K) == 0)
                    continue
                end
                for ij=1:size(parameter_values,1)
                    prediction = Dictionary(ij,:)';
                    X = [prediction ones(size(prediction))];
                    p = X\K;
                    prediction = X*p;
                    COST(ij) = sum((K-prediction).^2);
                end

                [cost(x),IX] = min(COST);
                best_p(x,:) = parameter_values(IX,:);
                [~,fit(x,:)] = function_handle(best_p(x,:),params);
            end
        end
        
    end
    
end