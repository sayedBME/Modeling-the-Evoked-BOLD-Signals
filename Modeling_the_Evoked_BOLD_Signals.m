clc
clear all
close all
TR = 1; % REPETITION TIME
t = 1:TR:20; % MEASUREMENTS
h = gampdf(t,6) + -.5*gampdf(t,10); % HRF MODEL
h = h/max(h); % SCALE HRF TO HAVE MAX AMPLITUDE OF 1
trPerStim = 30; % # TR PER STIMULUS
nRepeat = 2; % # OF STIMULUS REPEATES
nTRs = trPerStim*nRepeat + length(h);
impulseTrain0 = zeros(1,nTRs);
% VISUAL STIMULUS
impulseTrainLight = impulseTrain0;
impulseTrainLight(1:trPerStim:trPerStim*nRepeat) = 1;
% AUDITORY STIMULUS
impulseTrainTone = impulseTrain0;
impulseTrainTone(5:trPerStim:trPerStim*nRepeat) = 1;
% SOMATOSENSORY STIMULUS
impulseTrainHeat = impulseTrain0;
impulseTrainHeat(9:trPerStim:trPerStim*nRepeat) = 1;
% COMBINATION OF ALL STIMULI
impulseTrainAll = impulseTrainLight + impulseTrainTone + impulseTrainHeat;
%% SIMULATE VOXELS WITH VARIOUS SELECTIVITIES
visualTuning = [4 0 0]; % VISUAL VOXEL TUNING
auditoryTuning = [0 2 0]; % AUDITORY VOXEL TUNING
somatoTuning = [0 0 3]; % SOMATOSENSORY VOXEL TUNING
noTuning = [1 1 1]; % NON-SELECTIVE
beta = [visualTuning', ...
 auditoryTuning', ...
 somatoTuning', ...
 noTuning'];
%% EXPERIMENT DESIGN / STIMULUS SEQUENCE
D = [impulseTrainLight',impulseTrainTone',impulseTrainHeat'];
% CREATE DESIGN MATRIX FOR THE THREE STIMULI
X = conv2(D,h'); % X = D * h
X(nTRs+1:end,:) = []; % REMOVE EXCESS FROM CONVOLUTION
% DISPLAY STIMULUS AND DESIGN MATRICES

subplot(121);
imagesc(D);
colormap gray;
xlabel('Stimulus Condition')
ylabel('Time (TRs)');
title('Stimulus Train, D');
set(gca,'XTick',1:3); set(gca,'XTickLabel',{'Light','Tone','Heat'});

subplot(122);
imagesc(X);
xlabel('Stimulus Condition')
ylabel('Time (TRs)');
title('Design Matrix, X = D * h')
set(gca,'XTick',1:3); set(gca,'XTickLabel',{'Light','Tone','Heat'});
y0 = X*beta;

figure;
subplot(211);
imagesc(beta); colormap hot;
axis tight
ylabel('Condition')
set(gca,'YTickLabel',{'Visual','Auditory','Somato.'})
xlabel('Voxel');
set(gca,'XTick',1:4)
title('Voxel Selectivity, \beta')

subplot(212);
plot(y0,'Linewidth',2);
legend({'Visual Voxel','Auditory Voxel','Somato.Voxel','Unselective'});
xlabel('Time (TRs)'); ylabel('BOLD Signal');
title('Activity for Voxels with Different Stimulus Tuning')
set(gcf,'Position',[100 100 750 540])
subplot(211);
colorbar

%% SIMULATE NOISY VOXELS & ESTIMATE TUNIN
SNR = 5; % (APPROX.) SIGNAL-TO-NOISE RATIO
noiseSTD = max(y0(:))./SNR; % NOISE LEVEL FOR EACH VOXEL
noise = bsxfun(@times,randn(size(y0)),noiseSTD);
y = y0 + noise;
betaHat = inv(X'*X)*X'*y; % OLS
yHat = X*betaHat; % GLM PREDICTION

figure
subplot(211);
plot(y,'Linewidth',3);
xlabel('Time (s)'); ylabel('BOLD Signal');
legend({'Visual Voxel','Auditory Voxel','Somato. Voxel','Unselective'});
title('Noisy Voxel Responses');

subplot(212)
h1 = plot(y0,'Linewidth',3); 
hold on
h2 = plot(yHat,'-o');
legend([h1(end),h2(end)],{'Actual Responses','Predicted Responses'})
xlabel('Time (s)'); ylabel('BOLD Signal');
title('Model Predictions')
set(gcf,'Position',[100 100 750 540])

figure
subplot(211);
imagesc(beta);
colormap hot(5);
axis tight
ylabel('Condition')
set(gca,'YTickLabel',{'Visual','Auditory','Somato.'})
xlabel('Voxel');
set(gca,'XTick',1:4)
title('Actual Selectivity, \beta')

subplot(212)
imagesc(betaHat); 
colormap hot(5);
axis tight
ylabel('Condition')
set(gca,'YTickLabel',{'Visual','Auditory','Somato.'})
xlabel('Voxel');
set(gca,'XTick',1:4)
title('Noisy Estimated Selectivity')

drawnow
%% MODEL OF THE HRF
t = 0:.1:20;
hrfModel = gampdf(t,6) + -.5*gampdf(t,10);
% DISPLAY

figure
plot(t,hrfModel,'r','Linewidth',2);
hold on;
hb = plot(t,zeros(size(t)),'k--');
title('Model HRF');
xlabel('Time From Activity Onset (s)');
ylabel('BOLD Signal');
legend(hb,'Baseline')
%% FIR MODEL OF BOLD RESPONSE
% HRF AS MEASURED BY MRI SCANNER
TR = 1; % REPETITION TIME
t = 1:TR:20; % MEASUREMENTS
h = gampdf(t,6) + -.5*gampdf(t,10);
% CREATE A STIMULUS IMPULSE TRAIN
% FOR THREE SEPARATE STIMULUS CONDITIONS
trPerStim = 30; % # TR PER STIMULUS
nRepeat = 2;
nTRs = trPerStim*nRepeat + length(h);
impulseTrain0 = zeros(1,nTRs);
impulseTrainModel = impulseTrain0;
impulseTrainModel([1,24,28]) = 1;
boldModel = conv(impulseTrainModel,h);
boldModel(nTRs+1:end) = [];
% DISPLAY AN EXAMPLE OF FIR RESPONSE

figure
stem(impulseTrainModel*max(h),'k');
hold on;
plot(boldModel,'r','Linewidth',2); xlim([0,nTRs]);
title('HRF Convolved with Impulse Stimulus Train');
xlabel('Time (TRs)'); ylabel('BOLD Signal');
legend({'Activity Onset', 'Voxel Response'})
set(gcf,'Position',[100 100 750 380])
%% SIMULATE A BLOCK-DESIGN EXPERIMENT
% CREATE BLOCK STIMULUS TRAIN
blocks = repmat([ones(8,1);zeros(8,1)], round((nTRs-10)/16),1);
blockImpulseTrain = impulseTrain0;
blockImpulseTrain(1:numel(blocks)) = blocks;
boldBlock = conv(blockImpulseTrain,h);

% DISPLAY BOLD RESPONSES FROM BLOCK DESIGN
figure
stem(blockImpulseTrain*max(h),'k');
hold on;
plot(boldBlock(1:nTRs),'r','Linewidth',2); xlim([0,nTRs]);
title('Simulated Block Design');
xlabel('Time (TRs)');
ylabel('BOLD Signal');
legend({'Stimulus Train', 'Voxel Response'})
set(gcf,'Position',[100 100 750 380])
%% SIMULATE AN EVENT-RELATED EXPERIMENT
% VISUAL STIMULUS
impulseTrainLight = impulseTrain0;
impulseTrainLight(1:trPerStim:trPerStim*nRepeat) = 1;
% AUDITORY STIMULUS
impulseTrainTone = impulseTrain0;
impulseTrainTone(5:trPerStim:trPerStim*nRepeat) = 1;
% SOMATOSENSORY STIMULUS
impulseTrainPress = impulseTrain0;
impulseTrainPress(9:trPerStim:trPerStim*nRepeat) = 1;
% COMBINATION OF ALL STIMULI
impulseTrainAll = impulseTrainLight + impulseTrainTone + impulseTrainPress;

%% SIMULATE BOLD SIGNAL EVOKED BY EACH CONDITION
boldLight = conv(impulseTrainLight,h);
boldTone = conv(impulseTrainTone,h);
boldPress = conv(impulseTrainPress,h);
boldAll = conv(impulseTrainAll,h);
% DISPLAY STIMULUS ONSETS FOR EACH CONDITION

figure
subplot(211)
hold on
stem(impulseTrainLight,'k');
stem(impulseTrainTone,'b');
stem(impulseTrainPress,'g');
xlim([0,nTRs]);
xlabel('Time (TRs)'); ylabel('BOLD Signal');
legend({'Light Stimulus Train', 'Tone Stimulus Train','Press Stimulus Train'})
title('Impulse Trains for 3 Different Stimuli');

% DISPLAY COMBINATION OF BOLD RESPONSES FROM EACH CONDITION
subplot(212)
hold on;
plot(boldLight(1:nTRs),'k');
plot(boldTone(1:nTRs),'b');
plot(boldPress(1:nTRs),'g');
plot(boldAll(1:nTRs),'r','Linewidth',2);
xlim([0,nTRs]);
xlabel('Time (TRs)'); ylabel('BOLD Signal');
legend({'Response to Light','Response to Tone','Response to Press','Total Response'});
title('Simulation of BOLD Signal from Overlapping Stimuli');
set(gcf,'Position',[100 100 750 540])


