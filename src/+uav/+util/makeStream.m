function stream = makeStream(baseSeed, substream)
%MAKESTREAM Create a deterministic local random stream without touching rng.

arguments
    baseSeed (1,1) double {mustBeInteger, mustBePositive}
    substream (1,1) double {mustBeInteger, mustBePositive} = 1
end

stream = RandStream('Threefry', 'Seed', baseSeed);
stream.Substream = substream;
end
