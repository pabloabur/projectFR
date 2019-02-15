clear all
close all
pkg load ltfat
pkg load signal

%*******************
%**** Main variables
%*******************
tf = 9000;
dt = 0.05;
pulseWidth = 10; % ms
f = 10; % Hz
T = 100; % ms
t=[0:dt:tf];
tCut = tf-T/2; % First element of a interval to be placed at
               % the beginning of wave
sampleCut = tCut/dt;
nMN = 300; % Number of MNs to be activated
noiseMeanO = 25;
noiseMeanC = 60;
noiseStd = 100;
status = 0; % variable used to test whiteness
csiAmplitude = 100; % Half of desired value because max(square)=2

%*******************
%**** CSI input
%*******************
csiInput = (square(2*pi*f*t*1e-3, pulseWidth/T)+1)*csiAmplitude;
waveTail=csiInput(sampleCut:end);
csiInput=[waveTail,csiInput(1:sampleCut-1)];

% Saving data
data=csiInput';
save csi2.txt data

%*******************
%**** Independent noise
%*******************
tempINInput = (square(2*pi*f*t*1e-3, pulseWidth/(T))+1)/2;  
waveTail = tempINInput(sampleCut:end);
tempINInput = [waveTail,tempINInput(1:sampleCut-1)];
  
for k = 1:2 % 1 for open loop and 2 for closed loop
    inInput = [];
    for j = 1:nMN
      while (status==0)
        if (k==1)
            noiseInput = noiseStd*noise(length(tempINInput), 1,'white')+noiseMeanO;
        else
            noiseInput = noiseStd*noise(length(tempINInput), 1,'white')+noiseMeanC;
        endif
      
        % Test whiteness
        N = length(tempINInput);
        y=0;
        E=2*(N-2)/3;
        var=(16*N-29)/90;
        x=noiseInput;
        for i=2:N-1
          if ((x(i-1)<x(i) && x(i)>x(i+1)) || (x(i-1)>x(i) && x(i)<x(i+1)))
            y++;
          endif
        end
        z = abs((y-E)/sqrt(var));
        if (z<1.4)
          status = 1;
        endif
      endwhile
      
      inInput = vertcat(inInput, tempINInput.*noiseInput');
      status = 0;

      if (k==1)
        save ino2.txt inInput
      else
        save inc2.txt inInput
      endif
      % Removed headers manually, but I could just use csvwrite
    endfor
endfor
  
%*******************
%**** Plots
%*******************
%plot(t, csiInput)
%hold on
%plot(t, tempINInput)
