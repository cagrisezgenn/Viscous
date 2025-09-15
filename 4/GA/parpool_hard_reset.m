function parpool_hard_reset(nWorkers)
% PARPOOL_HARD_RESET parpool'u sıfırlayıp iş parçacıklarını kısıtlar.
% Eski işleri temizleyerek güvenli bir havuz açılışı sağlar.

if nargin<1 || isempty(nWorkers)
    nWorkers = feature('numcores');
else
    validateattributes(nWorkers, {'numeric'}, {'scalar','integer','positive'});
end

c = parcluster('Processes');

if ~isempty(c.Jobs)
    delete(c.Jobs);  % çökmüş işleri temizle
end

p = gcp('nocreate');
if isempty(p) || ~isvalid(p)
    parpool(c, min(nWorkers, c.NumWorkers));
end

pctRunOnAll maxNumCompThreads(1);          % CPU aşırı kullanımını engelle
pctRunOnAll set(0,'DefaultFigureVisible','off');
pctRunOnAll setenv('OMP_NUM_THREADS','1');
pctRunOnAll setenv('MKL_NUM_THREADS','1');

end
