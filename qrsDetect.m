function [pqrst, pqrstIdx] = qrsDetect(ecg,args)

    %%%% qrsDetect locates the P-QRS-T complexes in the given ECG recording

    %%%%   Inputs:
    %       *ecg--------------- ECG signal
    %
    %       *Rthresh----------- Amplitude threshold to use in detecting R peaks in ECG signal. Must have
    %                           the same units as sig.
    %       *t----------------- Time vector corresponding to sampling times (in sec)
    %       *fh----------------  [Approximate] Heart rate (in Hz)
    %
    %%%%   Outputs:
    %       *pqrst------------- A cell vector, where cells 1 through 5 correspond to the amplitude
    %                           values for the detected P,Q,R,S and T peaks respectively.
    %                           e.g. pqrst{4} contains the detected S peaks.
    %       *pqrstIdx---------- A cell vector where cells 1 through 5 correspond to the time
    %                           indexes of the detected P,Q,R,S and T peaks respectively.
    %                           e.g. pqrst{2} contains the time locations (in sec) of the detected Q peaks.


    %%% AUTHOR     : Crispin Foli
    %%% DATE       : Mar 10.2022.Su
    %%% Revisions  : 
    %%% FILENAME   : qrsDetect.m


    %% Parameter validation

    arguments
        ecg;
        args.t             =  1:numel(ecg); % Default to 'regular indexes' if time vector is not provided
        args.Rthresh (1,1) = 3*std(ecg); % Default to 3 times the std of the signal if R threshold is not provided
        args.fh (1,1)      = 1; % Default to 1Hz if heart rate is not provided
    end

    t       = args.t;
    Rthresh = args.Rthresh;
    tperiod = 0.8*(1/args.fh); % period of one heart beat cycle (in sec)..'subtract' a 20% 'buffer'
    Fs      = numel(ecg)/t(end);
    N       = floor(tperiod*Fs); % Number of samples in per heart beat cycle

    if size(ecg,1)>size(ecg,2), ecg=ecg'; end % Ensure ecg array is a row vector

     %% Detect peaks

     % Intuition:
     % 1. First find the R peaks and use that to extract the corresponding ECG waveforms
     % 2. The Q trough is the minimum value b/n the P and the R peaks
     % 3. The P peak is the maximum value b/n the beginning of the waveform and the Q trough
     % 4. The S trough is the minimum value within the entire waveform
     % 5. The T peak is the maximum value b/n the S trough and the end of the waveform:

    % Find R peaks
    pks     = find(islocalmax(ecg,'MinProminence',floor(Rthresh),"MinSeparation",N)); 

    pqrst{3}       = ecg(pks); % Extract detected R peaks (amplitudes) from original signal
    pqrstIdx{3}    = t(pks); % Extract temporal locations of peaks

    w       = -ceil(N/3):ceil(2*N/3); % Generic ECG window (around R peak)

    ecgIdxs =  bsxfun(@plus, w', pks); % Find indexes of all ECG waveforms
    ecgs    =  ecg(ecgIdxs); % Extract ECG waveforms

    % Find Q troughs
    [pqrst{2},pqrstIdx{2}] = min(ecgs(w<0,:)); forP = pqrstIdx{2}(1); % Q troughs
    pqrstIdx{2}            = t(ecgIdxs(sub2ind(size(ecgIdxs),pqrstIdx{2},1:size(ecgIdxs,2)))); % Corresponding temporal locations/indexes

    % Find S troughs
    [pqrst{4},pqrstIdx{4}] = min(ecgs); forT = pqrstIdx{4}(1); % S troughs
    pqrstIdx{4}            = t(ecgIdxs(sub2ind(size(ecgIdxs),pqrstIdx{4},1:size(ecgIdxs,2)))); % Corresponding temporal locations/indexes

    % Find P peaks
    [pqrst{1},pqrstIdx{1}] = max(ecgs(1:forP,:)); forQ = pqrstIdx{1}(1); % P peaks
    pqrstIdx{1}            = t(ecgIdxs(sub2ind(size(ecgIdxs),pqrstIdx{1},1:size(ecgIdxs,2)))); % Corresponding temporal locations/indexes

    % Update Q troughs (using detected P peaks)
    [pqrst{2},pqrstIdx{2}] = min(ecgs(forQ:find(w==0),:)); % Q troughs
    pqrstIdx{2}            = t(forQ-1+ecgIdxs(sub2ind(size(ecgIdxs),pqrstIdx{2},1:size(ecgIdxs,2))));  % Corresponding temporal locations/indexes

    % Find T peaks
    [pqrst{5},pqrstIdx{5}] = max(ecgs(forT:end,:)); % T peaks
    pqrstIdx{5}            = t(forT-1+ecgIdxs(sub2ind(size(ecgIdxs),pqrstIdx{5},1:size(ecgIdxs,2)))); % Corresponding temporal locations/indexes

end
