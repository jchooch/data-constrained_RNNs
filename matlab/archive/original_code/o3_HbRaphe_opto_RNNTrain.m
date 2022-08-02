function  o3_HbRaphe_opto_RNNTrain(simfilename, outputname, learn, nRunTot, nFree, nFreePre)
    load(simfilename);

    if (exist('nFreePre') == 0)
      nFreePre = 0;
    end

    N = double(N);

    tData_ind = (data_start_0ndx+1):(data_end_0ndx+1);
    tData = vol_starttimes(tData_ind);

    t = data_start_time:dt:data_end_time;
    lastchi2s = zeros(length(tData), 1);

    iModelSample = zeros(length(tData), 1);
    for i=1:length(tData)
        [tModelSample, iModelSample(i)] = min(abs(tData(i)-t));
    end
    
    Adata = datadyn./max(max(abs(datadyn)));
    Adata = Adata(:, tData_ind);
    Adata = min(Adata, 0.999);
    Adata = max(Adata, -0.999);

    stdData = std(reshape(Adata, N*length(tData_ind), 1));

    ampWN = sqrt(tauWN/dt);
    iWN = ampWN*randn(N, length(t));
    inputN = ones(N, length(t));
    inputC = zeros(N, length(t));
    inputB = zeros(N, length(t));
    
    weights_chr_l_LHb = wlHb0;
    weights_chr_r_LHb = wrHb0;
    weights_blue = wBlue0;
    
    stim_lHb_when = stim_mk_vec(stim_ndx_lHb, dtData, dt, data_start_time, data_end_time);
    stim_rHb_when = stim_mk_vec(stim_ndx_rHb, dtData, dt, data_start_time, data_end_time);
    stim_blue_when = or(stim_lHb_when, stim_rHb_when);
    
    for tt = 2: length(t)
        inputN(:, tt) = iWN(:, tt) + (inputN(:, tt - 1) - iWN(:, tt))*exp(-(dt/tauWN));
    end
    inputN = ampIn*inputN;
    
    Ttrial = size(Adata, 2);
    J = J0;
    R = NaN(N, length(t));
    AR = NaN(N+3, length(t)); %augmented activity for training input wts
    JR = zeros(N, 1);

    if (learn)
        PJ = P0 * eye(N+3, N+3);
    end

    [data_units, data_steps] = size(Adata);
    fprintf('data_steps: %d\n', data_steps);
    fprintf('dt: %f\n', dt);
    fprintf('length(t): %d\n', length(t));
    assert(data_steps >= length(t) * dt);

    chi2 = zeros(1, nRunTot);
    pVars = zeros(1, nRunTot);
    for nRun = 1:nRunTot
        H = Adata(:, 1);
        R(:, 1) = tanh(H);
        tLearn = 0;
        iLearn = 1;
        epoch_blue = 0;
        epoch_left = 0;
        epoch_right = 0;
        for tt = 2:length(t)
            tLearn = tLearn + dt;
            R(:, tt) = tanh(H);
            if stim_blue_when(tt) == 1
              epoch_blue = 1;
            end
            if stim_lHb_when(tt) == 1
              epoch_left = 1;
            end
            if stim_rHb_when(tt) == 1
              epoch_right = 1;
            end
            inputB(:, tt) = ampBlue*stim_blue_when(tt).*weights_blue;
            JR = J*R(:, tt) + inputN(:, tt) + inputB(:, tt);
            if ampChr2 ~= 0
                inputC(:, tt) = ampChr2 * stim_lHb_when(tt) * units_l_LHb' .* weights_chr_l_LHb;
                inputC(:, tt) = inputC(:, tt) + ampChr2 * stim_rHb_when(tt) * units_r_LHb' .* weights_chr_r_LHb;
                JR = JR + inputC(:, tt);
            end
            H = H + dt*(-H + JR )/tau;
            if (tLearn >= dtData)
                tLearn = 0;
                err = JR - Adata(:, iLearn);
                meanerr2 = mean(err.^2);
                chi2(nRun) = chi2(nRun) + meanerr2;
                lastchi2s(iLearn) = meanerr2;
                if ((learn) && (nRun <= nRunTot - nFree) && (nRun > nFreePre))
                    AR = [R(:, tt); epoch_left; epoch_right; epoch_blue]; %augmented Dyn variable
                    k = PJ*AR;
                    rPr = AR'*k;
                    c = 1.0/(1.0 + rPr);
                    PJ = PJ - c*(k*k');
                    if epoch_blue
                        weights_blue = weights_blue - c*err*k(end);
                    end
                    if ampChr2 ~= 0
                        if epoch_left
                            weights_chr_l_LHb = weights_chr_l_LHb - c*err*k(end-2);
                        end
                        if epoch_right
                            weights_chr_r_LHb = weights_chr_r_LHb - c*err*k(end-1);
                        end
                    end
                    J = J - c*err*k(1:N)';
                end
                iLearn = iLearn + 1;
                epoch_blue = 0;
                epoch_left = 0;
                epoch_right = 0;
            end
        end
        rModelSample = R(:, iModelSample);
        pVar = 1 - (norm(Adata - rModelSample, 'fro')/(sqrt(N*length(tData))*stdData)).^2;
        pVars(nRun) = pVar;
        fprintf('trial=%d pVar=%f chi2=%f\n', nRun, pVar, chi2(nRun));
    end

    varData = var(reshape(Adata, N*length(tData_ind), 1));
    chi2 = chi2/(sqrt(N*length(tData_ind))*varData);
    lastchi2s = lastchi2s/(sqrt(N*length(tData_ind))*varData);


    save(outputname, 'R', 'N', 't', 'chi2', 'Adata', 'tData', 'tData_ind', 'nRunTot', 'nFree', 'nFreePre', 'ampChr2', 'ampBlue', 'ampShock', 'ampIn', 'data_start_time', 'data_end_time', 'inputN', 'inputC', 'inputB', 'J', 'pVars', 'lastchi2s', 'datafile', 'simfilename', 'weights_blue', 'weights_chr_l_LHb', 'weights_chr_r_LHb', '-v7.3');
end
