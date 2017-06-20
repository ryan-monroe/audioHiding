% A quick demo of the work demonstrated in https://pdfs.semanticscholar.org/b964/4dbf393117b4c4bb91122cf7b08f71b71fbd.pdf

%BER=  0 @ 1bps
%BER=  1.26% @ 4bps
%BER=  8.4 % @ 10bps



useExponential=1;  %0: use delta functions ("straight delayed copy").  1: use exponential kernel
exponentialAlpha=4;  %higher numbers: exponential tapers off faster
exponentialSamples=5; %choose how long the filter should be here
delayMagnitudeScale=0.7; %what is the amplitude scaling on the echo?
inputIsWhiteNoise=0;  %if 1, use white noise instead of the dulcet tones of the Jackson 5
secondsToProcess=[]; %if empty ("[]"), choose a time that's useful.  otherwise, use seconds
windowTime=0.1; %seconds
delaySamples=[140, 200]; %how many samples delayed is the echo?  These should not be multiples of each other.  Optimal choice of delays will need more work, but shouldn't be hard


if inputIsWhiteNoise
    fs=44000; %make up a sampling rate
    if isempty(secondsToProcess) secondsToProcess=10;end
    nIn=10*fs
    nIn=round(fs*windowTime);
    din=wgn(nIn,1,0); %make white noise
else
    [audioThis,fs]=audioread('/home/rmonroe/Dropbox/Projects/audioStamping/ABC-Jackson5.ogg'); %ead file
    indFirst=find(audioThis(:,1),1,'first'); %the first bit was all zeros
    if isempty(secondsToProcess)
    din=audioThis(2*indFirst:end-2*indFirst,1); %only using one channel.  the left one?
    else
        indsIn=2*indFirst+(1:fs*secondsToProcess);
        din=audioThis(indsIn,1); %only using one channel.  the left one?
    end
end
windowLenSamples=round(windowTime*fs); %how many samples per bit?

nFrames=floor(numel(din)/windowLenSamples);
% dinTrunc=din(1:nFrames*windowLenSamples);
clear kernelIn din_bits

if useExponential
    kernelBase=exp(-1*(0:(exponentialSamples-1))*exponentialAlpha).'; %exponential kernel
    kernelIn=zeros(max(delaySamples)+exponentialSamples,2); %make the workspace we stamp it onto.  kernelIn(:,1) is for a zero, kernelIn(:,2) is for a one
    
    for i=1:2;kernelIn(1:exponentialSamples,i)=kernelIn(1:exponentialSamples,i)+kernelBase;end %the original copy
    for i=1:2;kernelIn((1:exponentialSamples)+delaySamples(i),i)=kernelIn((1:exponentialSamples)+delaySamples(i),i)+delayMagnitudeScale*kernelBase;end %the delayed copy
else
    for i=1:2;kernelIn(:,i)=zeros(1+max(delaySamples),1);kernelIn([1,delaySamples(i)+1],i)=[1, delayMagnitudeScale];end %delta functions instead
end

for i=1:2;din_bits(:,i)=conv(din,kernelIn(:,i)); end%convolution to encode the data

din_bits_trunc=din_bits(1:nFrames*windowLenSamples,:); %truncate to a even number of bits

din_bits_rs=reshape(din_bits_trunc,windowLenSamples,nFrames,2);

bitsIn=unidrnd(2,nFrames,1)   %generate some random bits

din_encoded=0*din_bits_rs(:,:,1); %sloppy but it gets the same type and size easily

for i=1:nFrames;din_encoded(:,i)=din_bits_rs(:,i,bitsIn(i));end  %choosing the third dimension to be 1 when the bit is '0', and choose it to be 2 when bit is '1'

assert(windowLenSamples>8*max(delaySamples)+2,'your delays are too long for your window length!  make the delays shorter or the window longer...');

binsAdd=[];
for mult=1:4; for addVal=-1:1; binsAdd(end+1,:)=delaySamples*mult+addVal;end;end %choose what bins of the cepstrum to add.  This is a less efficient version of using a matched filter.  but much faster and easier to make
    
%now decode
scoresOut=zeros(nFrames,2);
for i=1:nFrames
    dinThis=din_encoded(:,i);
    cepstrumOut=abs(ifft(log(fft(dinThis)).^2));
    for bitVal=1:2;scoresOut(i,bitVal)=sum(cepstrumOut(binsAdd(:,bitVal)));end  %primitive matched filter to produce scores.  a bit better than just choosing the single delay bin.
    i
end
[~,bitsChosen]=max(scoresOut,[],2); %choose the highest score as the winner!

bitsRight=bitsChosen==bitsIn;
bitErrorRate=1-mean(bitsRight)
%%
windowRmsPower=rms(din_encoded(:,:,1)); %rough estimate of power in window.  really this should be done on the non-convolved data, but whatever

binRanges=linspace(min(windowRmsPower),max(windowRmsPower),20);
histWrongs=hist(windowRmsPower(bitsRight==0),binRanges);
histRights=hist(windowRmsPower(bitsRight==1),binRanges);

figure
s=subplot(3,1,1);
scatter(windowRmsPower,bitsRight);
xlabel('windows RMS power (approx)')
ylabel('bit correct');
s(2)=subplot(3,1,2);
bar(binRanges,histWrongs);
xlabel('windows RMS power (approx)')
ylabel('number of events');
title('histogram of window powers, for bits which were guessed incorrectly')

s(3)=subplot(3,1,3);
bar(binRanges,histRights);
xlabel('windows RMS power (approx)')
ylabel('number of events');
title('histogram of window powers, for bits which were guessed correctly')

audioRaw=din_encoded(:);
ap=audioplayer(audioRaw,fs);
ap.playblocking();



