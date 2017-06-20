% A quick demo of the work demonstrated in https://pdfs.semanticscholar.org/b964/4dbf393117b4c4bb91122cf7b08f71b71fbd.pdf


useExponential=1;  %0: use delta functions ("straight delayed copy").  1: use exponential kernel
exponentialAlpha=4;  %higher numbers: exponential tapers off faster
exponentialSamples=5; %choose how long the filter should be here
delayMagnitudeScale=0.7; %what is the amplitude scaling on the echo?
inputIsWhiteNoise=0;  %if 1, use white noise instead of the dulcet tones of the Jackson 5

windowTime=0.1; %seconds
delaySamples=140; %how many samples delayed is the echo?


if inputIsWhiteNoise
    fs=44000; %make up a sampling rate
    nIn=round(fs*windowTime);
    din=wgn(nIn,1,0); %make white noise
else
    [audioThis,fs]=audioread('/home/rmonroe/Dropbox/Projects/audioStamping/ABC-Jackson5.ogg'); %ead file
    indFirst=find(audioThis(:,1),1,'first'); %the first bit was all zeros
    nIn=round(fs*windowTime); %how many samples are we processing
    
    din=audioThis(2*indFirst+(1:nIn)-1,1); %only using one channel.  the left one?
end


if useExponential
    kernel0=exp(-1*(0:(exponentialSamples-1))*exponentialAlpha).'; %exponential kernel
    kernelIn=zeros(delaySamples+exponentialSamples-1,1); %make the workspace we stamp it onto
    
    kernelIn(1:exponentialSamples)=kernelIn(1:exponentialSamples)+kernel0; %the original copy
    kernelIn((1:exponentialSamples)+delaySamples)=kernelIn(1:exponentialSamples)+delayMagnitudeScale*kernel0; %the delayed copy
else
    kernelIn=zeros(1+delaySamples,1);kernelIn([1,delaySamples+1])=[1, delayMagnitudeScale]; %delta functions instead
end

din_encoded=conv(din,kernelIn); %convolution to encode the data



figure;
subplot(3,1,1)
semilogy(abs(fft(din_encoded)))
title('output spectrum -- just a curiousity, not expected to extract information')

subplot(3,1,2)
plot(abs(ifft((abs(fft(din_encoded))).^2)))
title('autocorrelation spectrum');

subplot(3,1,3)
plot(abs(ifft(log(fft(din_encoded)).^2)))
title('cepstrum autocorrelation -- actual useful result comes out here');


