% Alex's KAP Simulation
% Simulates KAP diffusion with IFT delivery to tip and sink at base

function [] = KAPdiffusionsimulation()
    clear
    clc
    clf
    close all

    D = 1.7; %um2/sec
    L = 10; %um
    dt = 0.005; %sec
    MaxParticles = 500;
    IFTfreq = 1.3; %Train/sec, from meuller et al KAP-GFP paper
    %TrainLoad = 6; %KAP/train, from engel et al. 2009
    LArray = [];
    InFlagMeans = [];
    InFlagSTDs = [];
    
    %Account for train loading at different lengths
    Counter = 1;
    TrainLoad = [13, 12, 11, 10, 9, 8, 7, 6];
    
    for L = 5:12 % for flagellar lengths
        InFlagellaArray = [];
        for X = 1:10 %repeat each simulation this many times
            [TotalInjects, TotalExits, InFlagella] = KAPdiffusionfunc(D, L, dt, MaxParticles, IFTfreq, TrainLoad(Counter)); %main simulation
            InFlagellaArray = [InFlagellaArray, InFlagella]; %adding the number of KAP in the flagella to an array
            disp(['The number is ' num2str(X)]) % the number of the simulation trial
        end
        InFlagMeans = [InFlagMeans, mean(InFlagellaArray)];
        InFlagSTDs = [InFlagSTDs, std(InFlagellaArray)];
        LArray = [LArray, L];
        disp(L)
        Counter = Counter + 1;
    end
    
    
    figure(2)
    clf
    NoDiffTrains = LArray./2.*IFTfreq.*2.*TrainLoad; %The number of KAP in active transport, assuming 2 um/s velocity, and doubling number for both directions
    InFlagMeans = InFlagMeans + NoDiffTrains./2; %Adding to the number of KAP in the flagellum
    
    errorbar(LArray, InFlagMeans, InFlagSTDs)
    
    hold on
    plot(LArray, NoDiffTrains, 'r-') %plotting the comparison to active transport in both directions
    title('KAP-GFP accumulation in flagella by length, simulation data')
    legend('Returning by diffusion', 'Returning by directed transport')
    ylabel('Number of molecules in flagellum')
    xlabel('Flagellum length (um)')
    axis([0 inf 0 inf])
    disp(NoDiffTrains)
end

function [TotalInjects, TotalExits, InFlagella] = KAPdiffusionfunc(D, L, dt, MaxParticles, IFTfreq, TrainLoad)




Gridsize = sqrt(2.*D.*dt); %The grid step size is the mean square displacement in dt
MaxPosition = floor(L./Gridsize); %How big the grid is
TheGrid = zeros(1, MaxPosition); % Keeps track of number of molecules in each grid spot



Positions = zeros(1, MaxParticles); % A list of the positions of all active molecules.


IFTperiod = 1./IFTfreq; %How often to add new molecules
Timer = 9999; % Running timer, is set high at first to introduct new molecules at the start

changesbase = zeros(1, MaxParticles) - 1; 
Injections = 0; %number of injections so far
Exits = 0; % number of exits so far
for AA = 1:100000 %the number of time points in the simulation
    Timer = Timer + dt; %incrementing time
    Positions = sort(Positions); % Ordering positions, will be a bunch of zeros and then the grid positions of all active molecules
    NumZeros = nnz(~Positions); % number of zeros
    NumNonZeros = MaxParticles - NumZeros; %Keeping track of the number of active molecules
    Positions = [Positions(NumZeros+1:end), Positions(1:NumZeros)]; % Putting the positions of all active molecules first in the array, followed by zeros

    if Timer >= IFTperiod %When timer goes over IFTperiod:
        Timer = 0; %reset timer
        % Add a number of new molecules (TrainLoad) to the end of the
        % flagellum (MaxPosition), increasing the number of nonzeros in the
        % positions array
        Positions(NumNonZeros+1:NumNonZeros+TrainLoad) = MaxPosition;
        Injections = Injections + 1; %track the number of injections
    end

    Changes = changesbase + 2.*round(rand(1, length(changesbase))); %Create a random vector of -1 and +1



    Positions(1:NumNonZeros) = Positions(1:NumNonZeros)+Changes(1:NumNonZeros); %For each non-zero in the position array, add or subtract one (move in grid)
    Positions = min(Positions, MaxPosition); %Prevents molecules from moving out of the length of the flagellum
    Exits = Exits + nnz(~Positions(1:NumNonZeros)); %See how many values in the array moved from 1 to 0 (molecules exiting flagellum)
    


end
    A = nonzeros(Positions); 
    TheGrid = zeros(1, MaxPosition);
    for x = 1:length(A)
        TheGrid(A(x))= TheGrid(A(x))+1; % Tally the number of molecules in each grid position
    end
    figure(1)
    %plot(TheGrid, 'k')
    
    TotalInjects = TrainLoad.*Injections; % Tally the number of molecules injected
    TotalExits = Exits; % Tally the number of molecules that have exited. Theoretically, this should be very close to the number of injects.
    InFlagella = sum(TheGrid); % The number of molecules currently inside the flagellum.
    
end