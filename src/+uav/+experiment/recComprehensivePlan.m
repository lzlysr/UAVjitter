function [geometries,cases]=recComprehensivePlan(cfg)
%RECCOMPREHENSIVEPLAN 在任何算法结果产生前确定几何库和条件清单。
g=cfg.experiments.recComprehensive;
geometries=struct('id',{},'family',{},'worldAzElDeg',{},'randomIndex',{});
stream=uav.util.makeStream(cfg.random.baseSeed,g.geometrySeed);
for index=1:g.randomGeometryCount
    % 独立均匀下半球方向，不用算法性能筛选或拒绝不利几何。
    count=1+max(g.interfererCounts);
    directions=[360*rand(stream,count,1)-180,-asind(rand(stream,count,1))];
    addGeometry("random_world",directions,index);
end
ulaIndices=[]; upaIndices=[];
for desired=g.ulaDesiredGridDeg
    for first=g.ulaInterferer1GridDeg
        for second=g.ulaInterferer2GridDeg
            theta=[desired;first;second];
            elevation=cfg.experiments.recUlaGeometry.fixtureElevationDeg;
            directions=[asind(sind(theta)/cosd(elevation)),repmat(elevation,3,1)];
            ulaIndices(end+1)=addGeometry("ula_cartesian",directions,0); %#ok<AGROW>
        end
    end
end
for az=g.upaDesiredAzGridDeg
    for el=g.upaDesiredElGridDeg
        for iaz=g.upaInterfererAzGridDeg
            for iel=g.upaInterfererElGridDeg
                upaIndices(end+1)=addGeometry("upa_cartesian", ...
                    [az,el;iaz,iel;g.upaOtherInterfererDeg],0); %#ok<AGROW>
            end
        end
    end
end
% 旧反例只作锚点，不混入随机几何总体推断。
upaIndices(end+1)=addGeometry("anchor_default", ...
    [cfg.sources.desiredAzElDeg;cfg.sources.interfererAzElDeg],0);
upaIndices(end+1)=addGeometry("anchor_orthogonal", ...
    [cfg.experiments.recRootCause.orthogonalDesiredAzElDeg; ...
    cfg.experiments.recRootCause.orthogonalInterfererAzElDeg],0);
for angles={[80,-50,-20],[5,-31,-29]}
    elevation=cfg.experiments.recUlaGeometry.fixtureElevationDeg;
    directions=[asind(sind(angles{1}(:))/cosd(elevation)),repmat(elevation,3,1)];
    ulaIndices(end+1)=addGeometry("anchor_ula",directions,0);
end
cases=struct('id',{},'block',{},'arrayType',{},'geometryIndex',{},'snrDb',{}, ...
    'inrDb',{},'halfWidthDeg',{},'jitterStdDeg',{},'snapshotCounts',{}, ...
    'trials',{},'gridSize',{},'population',{});
for arrayType=g.arrayTypes
    if arrayType=="ULA", native=ulaIndices; grid=g.ulaGridSize; else, native=upaIndices; grid=g.upaGridSize; end
    if cfg.profile=="smoke", native=native(1:min(numel(native),g.randomGeometryCount)); end
    for geo=[1:g.randomGeometryCount,native]
        addCase("population",arrayType,geo,g.baseSnrDb,g.baseInrDb,g.baseHalfWidthDeg,0,Inf,1,grid,true);
    end
    for geo=1:g.randomGeometryCount
        addCase("snapshots",arrayType,geo,g.baseSnrDb,g.baseInrDb,g.baseHalfWidthDeg,0, ...
            g.snapshotCounts,g.trialCount,grid,false);
    end
    for geo=1:min(g.factorGeometryCount,g.randomGeometryCount)
        if isfield(g,'inrPairsDb')
            for pair=1:size(g.inrPairsDb,1)
                addCase("imbalance",arrayType,geo,g.baseSnrDb,g.inrPairsDb(pair,:),g.baseHalfWidthDeg,0, ...
                    g.factorSnapshotCounts,g.trialCount,grid,false);
            end
        end
        for snr=g.powerSnrDb
            for inr=g.powerInrDb
                addCase("power",arrayType,geo,snr,[inr,inr],g.baseHalfWidthDeg,0, ...
                    g.factorSnapshotCounts,g.trialCount,grid,false);
            end
        end
        for half=g.halfWidthsDeg
            addCase("sector",arrayType,geo,g.baseSnrDb,g.baseInrDb,half,0, ...
                g.factorSnapshotCounts,g.trialCount,grid,false);
        end
        for jitter=g.jitterStdDeg
            addCase("jitter",arrayType,geo,g.baseSnrDb,g.baseInrDb,g.baseHalfWidthDeg,jitter, ...
                g.factorSnapshotCounts,g.jitterTrials,grid,false);
        end
    end
    for geo=1:min(g.countGeometryCount,g.randomGeometryCount)
        for count=g.interfererCounts
            addCase("source_count",arrayType,geo,g.baseSnrDb,repmat(g.baseInrDb(1),1,count), ...
                g.baseHalfWidthDeg,0,Inf,1,grid,true);
        end
    end
    if arrayType=="ULA", grids=g.numericalUlaGrids; else, grids=g.numericalUpaGrids; end
    for geo=1:min(g.numericalGeometryCount,g.randomGeometryCount)
        for value=grids
            addCase("numerical",arrayType,geo,g.baseSnrDb,g.baseInrDb,g.baseHalfWidthDeg,0,Inf,1,value,true);
        end
    end
end
    function id=addGeometry(family,directions,randomIndex)
        id=numel(geometries)+1;
        geometries(id)=struct('id',id,'family',family,'worldAzElDeg',directions,'randomIndex',randomIndex);
    end
    function addCase(block,arrayType,geo,snr,inr,half,jitter,counts,trials,grid,population)
        id=numel(cases)+1;
        cases(id)=struct('id',id,'block',block,'arrayType',arrayType,'geometryIndex',geo, ...
            'snrDb',snr,'inrDb',inr,'halfWidthDeg',half,'jitterStdDeg',jitter, ...
            'snapshotCounts',counts,'trials',trials,'gridSize',grid,'population',population);
    end
end
