function result = exp08_gu_leshem2012_upa_transfer(cfg)
%EXP08_GU_LESHEM2012_UPA_TRANSFER 诊断原论文Rec向当前UPA场景的迁移。
%   本实验不声称是原论文复现，不修改固定排除区，也不加入时序或结构增强。
result = uav.experiment.runGuLeshem2012UpaTransfer(cfg);
end
