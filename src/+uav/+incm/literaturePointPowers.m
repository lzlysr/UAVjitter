function [Rin, d] = literaturePointPowers(rule, Rx, A, desiredCount, noise, lambda, cfg)
%LITERATUREPOINTPOWERS 同一候选字典下比较论文功率估计核心，禁止读取真值。
%   lin_literal：Lin 2021式(23)-(24)，不加载、不额外加噪声；L<M由上层标记。
%   lin_loaded_noise：本研究显式变体，加载Capon点功率并补噪声，不冒称原式。
%   yue_ls：Yue 2023式(46)-(50)的实系数最小二乘，无非负或L2约束。
%   nnls：本工程非负L2拟合。所有版本最终仅保留干扰候选。
Rx = (Rx+Rx')/2;
M = size(Rx,1);
assert(size(A,1)==M && desiredCount>=0 && desiredCount<=size(A,2));
rule = string(rule);
d = struct('rule',rule,'negativePowerCount',0,'dictionaryRank',0, ...
    'gramCondition',NaN,'relativeResidual',NaN,'usedTruth',false);
switch rule
    case {"lin_literal","lin_loaded_noise"}
        Rc = Rx;
        addedNoise = 0;
        if rule=="lin_loaded_noise"
            Rc = Rc + cfg.incm.staticCovfit.loadingRatio*real(trace(Rx))/M*eye(M);
            addedNoise = noise;
        elseif rcond(Rc)<=eps
            error('uav:study:NotApplicable','Lin未加载Capon需要可逆SCM。');
        end
        p = 1./real(sum(conj(A).*(Rc\A),1)).';
        assert(all(isfinite(p)) && all(p>0),'Capon点功率必须为正且有限。');
    case "yue_ls"
        % 使用实值最小二乘而不显式求(B^H B)^-1，避免平方条件数。
        % 字典秩亏时求最小范数解；不裁剪负功率，不暗中修复PSD。
        B = zeros(M*M,size(A,2));
        for q=1:size(A,2), B(:,q)=reshape(A(:,q)*A(:,q)',[],1); end
        target = reshape(Rx-noise*eye(M),[],1);
        D = [real(B);imag(B)]; y=[real(target);imag(target)];
        p = lsqminnorm(D,y);
        d.dictionaryRank=rank(D); d.gramCondition=cond(D'*D);
        d.relativeResidual=norm(D*p-y)/max(norm(y),eps);
        addedNoise = noise;
    case "nnls"
        [p,fit] = uav.incm.fitNonnegativePowers(Rx,A,noise,lambda);
        d.dictionaryRank=fit.dictionaryRank;
        d.relativeResidual=fit.relativeResidual;
        addedNoise = noise;
    otherwise
        error('uav:study:UnknownPowerRule','未知的论文功率估计核心。');
end
Ai = A(:,desiredCount+1:end); pi = p(desiredCount+1:end);
Rin = (Ai.*pi.')*Ai'+addedNoise*eye(M);
Rin = (Rin+Rin')/2;
d.powers=p; d.interfererPowers=pi;
d.negativePowerCount=sum(p<0);
d.minimumEigenvalue=min(real(eig(Rin)));
d.noisePower=addedNoise;
end
