function obsShort = truncatePeriodObservation(obsLong, L)
%TRUNCATEPERIODOBSERVATION 截取配对观测的前L个快拍并重算样本协方差。
%   除X、pilot和真值干扰符号外，其余观测与真值字段保持不变，因此Oracle
%   与最终指标仍可读取完整的同周期物理真值；输入结构不会被修改。

validateInputs(obsLong, L);
obsShort = obsLong;
obsShort.X = obsLong.X(:, 1:L);
obsShort.pilot = obsLong.pilot(:, 1:L);
obsShort.truth.interfererSymbols = obsLong.truth.interfererSymbols(:, 1:L);
obsShort.Rx = obsShort.X * obsShort.X' / L;
end

function validateInputs(obsLong, L)
%VALIDATEINPUTS 保证截取的三类快拍字段具有一致的原始长度。
required = {'X', 'Rx', 'pilot', 'truth'};
if ~isstruct(obsLong) || ~all(isfield(obsLong, required)) || ...
        ~isstruct(obsLong.truth) || ...
        ~isfield(obsLong.truth, 'interfererSymbols')
    error('uav:scenario:InvalidPeriodObservation', ...
        'obsLong must be a generated period observation with isolated truth.');
end
snapshotCount = size(obsLong.X, 2);
if ~isnumeric(L) || ~isscalar(L) || ~isfinite(L) || L < 1 || ...
        L ~= floor(L) || L > snapshotCount
    error('uav:scenario:InvalidTruncationLength', ...
        'L must be a positive integer no greater than the original snapshot count.');
end
if ~ismatrix(obsLong.X) || isempty(obsLong.X) || ...
        ~isequal(size(obsLong.Rx), [size(obsLong.X, 1), size(obsLong.X, 1)]) || ...
        size(obsLong.pilot, 2) ~= snapshotCount || ...
        size(obsLong.truth.interfererSymbols, 2) ~= snapshotCount
    error('uav:scenario:InvalidPeriodObservation', ...
        'Observation snapshots, covariance, pilot, and interferer symbols must agree.');
end
end
