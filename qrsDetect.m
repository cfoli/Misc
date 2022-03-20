function [pqrst, pqrstIdx] = qrsDetect(ecg,Rthresh,t,Fs,tperiod)

    %%%% qrsDetect locates the PQRST complexes in the given ECG recording

    %%%%   Inputs:
    %       *ecg--------------- ECG signal
    %
    %       *Rthresh----------- Amplitude threshold to use in detecting R peaks in ECG signal. Must have
    %                           the same units as sig.
    %       *t----------------- Time vector corresponding to sampling times (in sec)
    %       *Fs---------------- Sampling frequency of ECG (in Hz)
    %
    %%%%   Outputs:
    %       *pqrst------------- A cell vector, where cells 1 through 5 correspond to the amplitude
    %                           values for the detected P,Q,R,S and T peaks respectively.
    %                           e.g. pqrst{4} contains the detected S peaks.
    %       *pqrstIdx---------- A cell vector where cells 1 through 5 correspond to the time
    %                           indexes of the detected P,Q,R,S and T peaks respectively.
    %                           e.g. pqrst{2} contains the time locations (in sec) of the detected Q peaks.

    %% Parameters

    arguments
        ecg;
        Rthresh;
        t;
        Fs;
        tperiod; % period of one heart beat (in sec)
    end

    %% Detect peaks

    % Intuition:
    % 1. First find the R peaks and use that to extract the corresponding ECG waveforms
    % 2. The Q trough is the minimum value b/n the P and the R peaks
    % 3. The P peak is the maximum value b/n the beginning of the waveform and the Q trough
    % 4. The S trough is the minimum value within the entire waveform
    % 5. The T peak is the maximum value b/n the S trough and the end of the waveform:

    %     bpm2Hz = (tperiod/60); % Heart rate in Hz
    %     tperiod = 1/bpm2Hz; % Sampling period (in sec)

    % Find R peaks
    pks     = find(islocalmax(ecg,'MinProminence',floor(Rthresh), ...
        "MinSeparation",floor(0.8*tperiod*Fs))); % find threshold crossings of interest

    pqrst{3}       = ecg(pks); % extract detected peaks (amplitudes) from original signal
    pqrstIdx{3}    = t(pks); % Extract temporal locations of peaks

    %     N       = 0.75*Fs; % Limit width of ECG waveform to 0.75sec
    N       = tperiod*Fs; % Number of samples in per heart beat cycle
    w       = -ceil(N/3):ceil(2*N/3); % Generic ECG window (around R peak)

    ecgIdxs =  bsxfun(@plus, w', pks); % Find indexes of all ECG waveforms
    ecgs    =  ecg(ecgIdxs); % Extract ECG waveforms

    % Find Q troughs
    [pqrst{2},pqrstIdx{2}] = min(ecgs(w<0,:)); forP = pqrstIdx{2}(1);
    pqrstIdx{2}            = t(ecgIdxs(sub2ind(size(ecgIdxs),pqrstIdx{2},1:size(ecgIdxs,2))));

    % Find S troughs
    [pqrst{4},pqrstIdx{4}] = min(ecgs); forT = pqrstIdx{4}(1);
    pqrstIdx{4}            = t(ecgIdxs(sub2ind(size(ecgIdxs),pqrstIdx{4},1:size(ecgIdxs,2))));

    % Find P peaks
    [pqrst{1},pqrstIdx{1}] = max(ecgs(1:forP,:)); forQ = pqrstIdx{1}(1);
    pqrstIdx{1}            = t(ecgIdxs(sub2ind(size(ecgIdxs),pqrstIdx{1},1:size(ecgIdxs,2))));

    % Update Q troughs (using detected P peaks)
    [pqrst{2},pqrstIdx{2}] = min(ecgs(forQ:find(w==0),:)); 
    pqrstIdx{2}            = t(forQ-1+ecgIdxs(sub2ind(size(ecgIdxs),pqrstIdx{2},1:size(ecgIdxs,2))));

    % Find T peaks
    [pqrst{5},pqrstIdx{5}] = max(ecgs(forT:end,:));
    pqrstIdx{5}            = t(forT-1+ecgIdxs(sub2ind(size(ecgIdxs),pqrstIdx{5},1:size(ecgIdxs,2))));

end