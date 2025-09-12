function parpool_hard_reset(nWorkers)
% PARPOOL_HARD_RESET parpool'u sıfırlayıp iş parçacıklarını kısıtlar.
% Eski işleri temizleyerek güvenli bir havuz açılışı sağlar.

    if nargin<1 || isempty(nWorkers), nWorkers = feature('numcores'); end

    try
        c = parcluster('Processes');
    catch ME
        warning('[parpool_hard_reset] küme alınamadı: %s', ME.message);
        try
            parpool('Threads');
        catch ME2
            warning('[parpool_hard_reset] yedek havuz açılamadı: %s', ME2.message);
        end
        return;
    end

    %% İş Temizliği
    try
        if ~isempty(c.Jobs)
            delete(c.Jobs);  % çökmüş işleri temizle
        end
    catch ME
        warning('[parpool_hard_reset] iş temizliği başarısız: %s', ME.message);
    end

    %% Havuz Açma
    try
        p = gcp('nocreate');
        if isempty(p) || ~isvalid(p)
            parpool(c, min(nWorkers, c.NumWorkers));
        end
    catch ME
        warning('[parpool_hard_reset] havuz açma başarısız: %s', ME.message);
        try
            parpool('Threads');
        catch ME2
            warning('[parpool_hard_reset] yedek havuz açılamadı: %s', ME2.message);
        end
    end

    %% Thread Sınırlandırma
    try
        pctRunOnAll maxNumCompThreads(1);          % CPU aşırı kullanımını engelle
        pctRunOnAll set(0,'DefaultFigureVisible','off');
        try
            pctRunOnAll setenv('OMP_NUM_THREADS','1');
            pctRunOnAll setenv('MKL_NUM_THREADS','1');
        catch ME
            warning('[parpool_hard_reset] ortam değişkenleri ayarlanamadı: %s', ME.message);
        end
    catch ME
        warning('[parpool_hard_reset] thread sınırlandırma başarısız: %s', ME.message);
    end
end

