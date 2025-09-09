function parpool_hard_reset(nWorkers)
% PARPOOL_HARD_RESET parpool'u sıfırlayıp iş parçacıklarını kısıtlar.
% Eski işleri temizleyerek güvenli bir havuz açılışı sağlar.

    if nargin<1 || isempty(nWorkers), nWorkers = feature('numcores'); end

    %% Havuz Temizliği
    try
        c = parcluster('Processes');
        if ~isempty(c.Jobs), delete(c.Jobs); end     % crash dump'lı işleri temizle
        p = gcp('nocreate');
        if isempty(p) || ~isvalid(p)
            parpool(c, min(nWorkers, c.NumWorkers));
        end

        %% Thread Sınırlandırma
        % oversubscription önlemleri
        pctRunOnAll maxNumCompThreads(1);
        pctRunOnAll set(0,'DefaultFigureVisible','off');
        try
            pctRunOnAll setenv('OMP_NUM_THREADS','1');
            pctRunOnAll setenv('MKL_NUM_THREADS','1');
        catch, end
    catch ME
        warning('[parpool_hard_reset] %s', ME.message);
        try
            p = gcp('nocreate'); if ~isempty(p), return; end
            parpool('Threads');
        catch, end
    end
end

