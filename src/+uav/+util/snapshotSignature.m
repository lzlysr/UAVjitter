function signature = snapshotSignature(X)
%SNAPSHOTSIGNATURE Return a deterministic observation identity diagnostic.

if ~isnumeric(X) || ~ismatrix(X) || isempty(X) || any(~isfinite(X), 'all')
    error('uav:util:InvalidSnapshot', 'X must be a finite nonempty matrix.');
end
signature = [sum(real(X), 'all'), sum(imag(X), 'all'), ...
    norm(X, 'fro'), real(sum(abs(X).^2, 'all'))];
end
