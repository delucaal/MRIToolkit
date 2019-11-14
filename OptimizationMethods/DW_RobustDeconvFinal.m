%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%%%% Alberto De Luca - alberto@isi.uu.nl $%%%%%% Distributed under the terms of LGPLv3  %%%
%%% Distributed under the terms of LGPLv3  %%%
function [LSQNONNEG_CL_DECONV,goodpoints] = DW_RobustDeconvFinal(Dictionary,NoisySimulatedSignal,kappa,options)

goodpoints = ones(size(NoisySimulatedSignal)) == 1;%squeeze(OUTLIERS(x,y,z,:))==0;
W = diag(NoisySimulatedSignal);
done = 0;
iters = 0;
if(nargin < 3)
    kappa = 3;
    options = optimset('TolX',1e-2);
end

while done == 0
    iters = iters+1;
    new_outliers_hes = 0;
    new_outliers_hos = 0;
    for iter=1:4
        LSQNONNEG_CL_DECONV = lsqnonneg(Dictionary(goodpoints,:),NoisySimulatedSignal(goodpoints),options);
        residuals = NoisySimulatedSignal-Dictionary*LSQNONNEG_CL_DECONV;
        sd = mad(residuals)*1.4826;
        bad_guys = abs(residuals-mean(residuals))>kappa*sd;
        if(sum(bad_guys(goodpoints == 1)==1) > 0)
            new_outliers_hes = 1;
        end
        goodpoints = goodpoints & bad_guys == 0;
        %             disp(['Outliers (HeS): ' num2str(length(find(goodpoints==0)))]);
    end
    for iter=1:4
        LSQNONNEG_CL_DECONV = lsqnonneg(W(goodpoints,goodpoints)*Dictionary(goodpoints,:),W(goodpoints,goodpoints)*NoisySimulatedSignal(goodpoints),options);
        W = diag((Dictionary*LSQNONNEG_CL_DECONV));
        residuals = NoisySimulatedSignal-Dictionary*LSQNONNEG_CL_DECONV;
        sd = mad(residuals)*1.4826;
        bad_guys = abs(residuals-mean(residuals))>kappa*sd;
        if(sum(bad_guys(goodpoints == 1)==1) > 0)
            new_outliers_hos = 1;
        end
        goodpoints = goodpoints & bad_guys == 0;
        %             disp(['Outliers (HoS): ' num2str(length(find(goodpoints==0)))]);
    end
    if(new_outliers_hos == 0 && new_outliers_hes == 0)
        done = 1;
    end
    if(iters >= 100)
        done = 1;
    end
end

LSQNONNEG_CL_DECONV = lsqnonneg(Dictionary(goodpoints,:),NoisySimulatedSignal(goodpoints),options);

end