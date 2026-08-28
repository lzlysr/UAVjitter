function stream = makeStream(baseSeed, substream)
%MAKESTREAM 创建确定性的局部Threefry随机流。
%   baseSeed固定整个研究工程，substream区分实验层级。返回独立RandStream，
%   不调用rng，因此不会污染用户MATLAB会话的全局随机状态。

arguments
    baseSeed (1,1) double {mustBeInteger, mustBePositive}
    substream (1,1) double {mustBeInteger, mustBePositive} = 1
end

stream = RandStream('Threefry', 'Seed', baseSeed);
stream.Substream = substream;
end
