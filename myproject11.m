% chapter 11 project MATLAB for neuroscientists book
%hadi km shahrivar 1400


% load audio files
clear ;close all;
num_peak = 5 %choosing 5 peaks (formants)
figure('Position', [100 100 1200 800])
LIMIT = 1000 %5000 is good
smoothing = 20
show = true

for audio_num=1:8
    t = sprintf('voice-sample-%d.wav', audio_num);
    [amp, fs] = audioread(t);
    % sound(amp,fs);
    
    L = length(amp); %length of audio
    Y = fft(amp); % power(amplitude?) of audio
    
%     NyLimit = 5000; %manually changed (not true Nyquist)
%     P = abs(Y/L); P = P(1:NyLimit); %replaced NyLimit instead of L/2+1
%     plot(P(1:NyLimit))
%     pause
    
    NyLimit = fs/ 2;
    F = linspace(0,1,LIMIT)*NyLimit; %??????
    P = (Y(1:LIMIT).*conj(Y(1:LIMIT)));
    plot(F, P);
    pause
    
    %smoothing
    smoothed = smooth(P,smoothing);
    plot(F,smoothed)
    hold on
    % [pks,locs] = findpeaks(smoothed);
    % findpeaks(smoothed)
    % text(locs+.02,pks,num2str((1:numel(pks))'))
    %     pause(2)
    [psor,lsor] = findpeaks(smoothed,'SortStr','descend');
    freq(1:length(lsor)) = lsor;
    lsor = F(lsor);
    plot(lsor(1:num_peak),psor(1:num_peak), 'v')
    text(lsor(1:num_peak)+.02,psor(1:num_peak),num2str((1:numel(psor(1:num_peak)))'))
    if show
        disp([num2str(num_peak),' Freq with most powers: ', num2str(freq(1:num_peak)), '  '])
    end
    pause
    hold off
end





%for peak analysis see this link
...https://www.mathworks.com/help/signal/ug/peak-analysis.html;jsessionid=6a83275df8ed7526c0c18cac537c
