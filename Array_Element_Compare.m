clear
clc
close all

%look in CostFunctions folder in same directory as this .m file for
%implementations of the Problem.m class
addpath ./CostFunctions 
load('patch_patterns.mat')

el_pattern = 10.^(uStripPatch27GHzElementPattern.dBGainTotalFreq27GHzPhi900000000000002deg./10);
theta = uStripPatch27GHzElementPattern.Thetadeg;
el4_pattern = 10.^(uStripPatch27GHz1x4ArrayPattern.dBGainTotal_1./10);
el8_pattern = 10.^(uStripPatch27GHz1x8ArrayPattern.dBGainTotal_Phi_90Freq27GHzPhi900000000000002deg./10);


c = physconst('LightSpeed');
Dk = 2.33; %relative permittivity
freq0 = 27e9; %Operating frequency of 28 GHz
lambda = c/freq0/sqrt(Dk);

%% 
for n = [3,7]
    % n = (# of array elements) - 1
    pos = (1:n).*5.56e-3; %uniform element spacing of 5.56mm
    
    %construct SLL class
    SLLFunc = SLLCost(n);
    SLLFunc.theta = theta;
    SLLFunc.Dk = Dk;
    SLLFunc.freq0 = freq0;
    
    %calculate array factor
    AF = SLLFunc.getAF(pos,1);
    
    %use element pattern from HFSS to calculate the array pattern
    pattern_calc = (abs(AF).^2).*el_pattern;
    
    %generate plots to compare with HFSS array simulation
    if n == 3
        hfss_pattern = el4_pattern;
    else
        hfss_pattern = el8_pattern;
    end
    figure()
    subplot(1,2,1)
    sgtitle(sprintf('%d Element Array',n+1))
    hold on 
    plot(theta, db(pattern_calc,'power'), 'DisplayName', 'Calculated')
    plot(theta, db(hfss_pattern,'power'), 'DisplayName', 'HFSS')
    xlabel('Angle (\circ)', 'Interpreter','tex')
    xlim([0 360])
    ylabel('Magnitude (dB)')
    legend
    hold off
    
    subplot(1,2,2)
    p = polarplot(deg2rad(theta), db(pattern_calc,'power'), deg2rad(theta), db(hfss_pattern,'power'));
    ax = gca;
    ax.ThetaZeroLocation = 'top';
    rlim([-35 25])
    rholabel
    legend('Caclulated','HFSS','Location','SouthOutside')
end