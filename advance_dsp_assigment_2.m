clear;close all;clc;
lowFreqs=[697,770,852,941];highFreqs=[1209,1336,1477];
dtmfMap={'1','2','3';'4','5','6';'7','8','9';'*','0','#'};
filename='dialed_number.mp3';
[signal,fs]=audioread(filename);signal=signal(:,1);
signal=signal/max(abs(signal));
lowCutoff=600;highCutoff=1700;order=4;
[b,a]=butter(order,[lowCutoff,highCutoff]/(fs/2),'bandpass');
filtered_signal=filtfilt(b,a,signal);
f1=figure; t=(0:length(signal)-1)/fs;
plot(t,signal,'b-',t,filtered_signal,'r--'); xlabel('Time (s)'); ylabel('Amplitude'); legend('Original','Filtered'); title('Raw vs. Filtered Dialing Sound Signal');
print(f1,'Raw_vs_Filtered.jpg','-djpeg','-r300');
signal=filtered_signal;
windowLength=round(0.03*fs); stepSize=round(0.01*fs);
numSteps=floor((length(signal)-windowLength)/stepSize)+1; energy=zeros(1,numSteps);
for k=1:numSteps
    idx=(1:windowLength)+(k-1)*stepSize;
    energy(k)=sum(signal(idx).^2);
end
f2=figure; tEnergy=((1:numSteps)*stepSize+windowLength/2)/fs;
plot(tEnergy,energy); xlabel('Time (s)'); ylabel('Energy'); title('Short-Time Energy Envelope');
print(f2,'Energy_Envelope.jpg','-djpeg','-r300');
energyThreshold=0.1*max(energy); toneIdx=energy>energyThreshold;
regions=[]; inRegion=false;
for k=1:length(toneIdx)
    if toneIdx(k)&&~inRegion
        regionStart=(k-1)*stepSize+1; inRegion=true;
    elseif ~toneIdx(k)&&inRegion
        regionEnd=(k-1)*stepSize+windowLength;
        regions=[regions; regionStart,regionEnd]; inRegion=false;
    end
end
if inRegion, regionEnd=length(signal); regions=[regions; regionStart,regionEnd]; end
fprintf('Detected %d tone regions.\n',size(regions,1));
detectedDigits='';
for r=1:size(regions,1)
    seg=signal(regions(r,1):regions(r,2));
    seg=seg.*hamming(length(seg));
    N=2^nextpow2(length(seg)); Y=fft(seg,N);
    f=fs*(0:(N/2))/N; Y_mag=abs(Y(1:N/2+1));
    f3=figure; plot(f,Y_mag); xlabel('Frequency (Hz)'); ylabel('Magnitude');
    title(sprintf('Spectrum of tone region %d',r));
    print(f3,sprintf('Spectrum_Tone_Region_%d.jpg',r),'-djpeg','-r300');
    searchRange=f>600 & f<1600; f_search=f(searchRange); Y_search=Y_mag(searchRange);
    [pks,locs]=findpeaks(Y_search,f_search,'SortStr','descend');
    if length(pks)<2, fprintf('Warning: less than 2 peaks detected in region %d.\n',r); continue; end
    dominantFreqs=sort(locs(1:2));
    fprintf('Region %d: Detected frequencies: %.1f Hz, %.1f Hz\n',r,dominantFreqs(1),dominantFreqs(2));
    [~,idxLow]=min(abs(lowFreqs-dominantFreqs(1)));
    [~,idxHigh]=min(abs(highFreqs-dominantFreqs(2)));
    if abs(lowFreqs(idxLow)-dominantFreqs(1))>50 || abs(highFreqs(idxHigh)-dominantFreqs(2))>50
        [~,idxLow_swap]=min(abs(lowFreqs-dominantFreqs(2)));
        [~,idxHigh_swap]=min(abs(highFreqs-dominantFreqs(1)));
        if abs(lowFreqs(idxLow_swap)-dominantFreqs(2))<50 && abs(highFreqs(idxHigh_swap)-dominantFreqs(1))<50
            idxLow=idxLow_swap; idxHigh=idxHigh_swap;
        else
            fprintf('Region %d: Frequency mapping error.\n',r); continue;
        end
    end
    detectedDigit=dtmfMap{idxLow,idxHigh};
    detectedDigits=[detectedDigits detectedDigit];
end
fprintf('\nDetected Phone Number: %s\n',detectedDigits);
