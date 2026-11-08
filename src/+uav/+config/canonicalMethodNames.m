function names = canonicalMethodNames(names)
%CANONICALMETHODNAMES 将历史方法标识转换为当前统一名称，保持顺序和维度。
%   兼容旧MAT/脚本；不在此处去重，调用方必须拒绝新旧同义名重复选择。
names = string(names);
aliases = ["gu_leshem2012_rec","static","CovFit","CovFit_local","static_incm_mvdr"];
canonical = ["incm_integral_origin","incm_integral_processed", ...
    "incm_covfit","incm_covfit_refined","incm_covfit"];
[found,index] = ismember(lower(names),lower(aliases));
names(found) = canonical(index(found));
end
