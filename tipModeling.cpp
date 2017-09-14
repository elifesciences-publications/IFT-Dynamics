// tip1: Tip modeling.
// Given all the parameters. Calculate the average number of returned particles
// 		and the average return time for 1st, 2nd, ... particles

// tip3: Scan over number of particle in the tip
// tip4: Alex's attempt to include a photobleaching rate parameter.
// tip5: Alex's change to how new bright anterograde trains are introduced.
// tip6: Alex's change to loop over number of anterograde trains allowed through gate.
// tip7: Alex's change to give output of distribution of number of returning trains.
// tipModeling: Version 20170831.

#include <iostream>
#include <fstream>

#define MAXNP 100 // maximum number of particle returned


int main()
{
	double incomeFrequency = 1.15; //Anterograde frequency from our calculations, 1.15
	double returnFrequency = 1.31; //Retrograde frequency from our calculations, 1.31
	double timeResolution = 0.01;
	double bleachFrequency = 0.073; // 0.073; // bleaching events per molecule per second
	double travelTime = 5.0; // travel time for anterograde train to reach tip, in seconds, plus tip waiting time	
	
	

	

	// For number of particles on trains, used Wallace Marshall's KAP-GFP bleaching paper and his signal comparison.
	// He saw an average of 6 KAP-GFP per train. Since IFT27-GFP is slightly higher, we estimate 7 for anterograde and 6 for retrograde
	// to account for the difference in frequencies.
	int numberParticleAnt = 6; 
	int numberParticleRet = 6;
	int numberParticleTip = 30; //Estimated from the tip intensity data I took, as compared to IFT27 trains. This changes later in code...
	
	double incomeFreqPerFrame = incomeFrequency*timeResolution;
	double returnFreqPerFrame = returnFrequency*timeResolution;
	double bleachFreqPerFrame = bleachFrequency*timeResolution;

	int numberRuns =10000;
	
	int countReturnTrain[MAXNP];
	double timeReturnTrain[MAXNP];
	double numParticleReturnTrain[MAXNP];
	int returnNumberDist[MAXNP];
	
	int countIncomeTrain;
	double numberReturnTrain;
	std::ofstream outFile;
	outFile.open("Tip7output.txt");

	outFile << "Anterograde train frequency (per sec): " << incomeFrequency << std::endl;
	outFile << "Retrograde train frequency (per sec): " << returnFrequency << std::endl;
	outFile << "Bleaching rate (per molecule per second): " << bleachFrequency << std::endl;
	outFile << "Anterograde travel time (sec): " << travelTime << std::endl;
	outFile << "Number of particles on anterograde train: " << numberParticleAnt << std::endl;
	outFile << "Number of particles on anterograde train: " << numberParticleRet << std::endl;
	outFile << std::endl;

	for (int trainsThruGate = 1; trainsThruGate <= 3; trainsThruGate++) // number of trains allowed through gate
	{
		outFile << "Allowing " << trainsThruGate << " trains through gate:" << std::endl;
		std::cout << "Allowing " << trainsThruGate << " trains through gate:" << std::endl;
		for (int i = 0; i < MAXNP; i++)
		{
			countReturnTrain[i] = 0;
			timeReturnTrain[i] = 0;
			numParticleReturnTrain[i] = 0;
			returnNumberDist[i] = 0;
		}
		countIncomeTrain = 0;
		numberReturnTrain = 0;



		srand(10000);

		for (int k = 8; k <= 13; k = k + 2) //For different values of number of tip particles...
		{
			numberParticleTip = k;

			for (int i = 0; i < numberRuns; i++) // For the number of repeats of this simulation...
			{

				int numberBrightParticle = numberParticleAnt; // number of bright particles arriving at the tip
				numberBrightParticle -= int(numberParticleAnt*bleachFrequency*travelTime); //Subtracts bleached molecules before they get to the tip
				int numberAllParticle = numberParticleTip + numberParticleAnt;
				int brightTrainsLeft = trainsThruGate - 1;

				int tempNumberReturnTrain = 0; //number of bright trains returning


				for (int i = 0; i < 10000; i++) //For 10,000 instants, each 0.01 seconds
				{
					// bleaching at the tip
					int tempNumberBleachedTip = 0;
					for (int m = 0; m < numberBrightParticle; m++)
					{
						if (static_cast<double>(rand()) / RAND_MAX < bleachFreqPerFrame)
						{
							tempNumberBleachedTip += 1;
						}
					}
					numberBrightParticle -= tempNumberBleachedTip;

					// retrograd particle
					if (static_cast<double>(rand()) / RAND_MAX < returnFreqPerFrame) // If a retrograde train leaves at this instant
					{
						int numberBrightParticleReturn = 0;

						for (int j = 0; j < numberParticleRet; j++) // For every particle that makes up a retrograde train...
						{
							if (static_cast<double>(rand()) / RAND_MAX*numberAllParticle < numberBrightParticle) //If it's a bright particle...
							{
								numberBrightParticleReturn += 1;
								numberBrightParticle -= 1;
								//std::cout << "yes ";
							}
							//else std::cout << "no ";
							if (numberAllParticle > 0) numberAllParticle -= 1; // As long as there are particles left at tip...
						}

						if (numberBrightParticleReturn > 0) // If this retrograde train has any bright particles...
						{
							countReturnTrain[tempNumberReturnTrain] += 1;
							timeReturnTrain[tempNumberReturnTrain] += (i + 1)*timeResolution;
							numParticleReturnTrain[tempNumberReturnTrain] += numberBrightParticleReturn;

							tempNumberReturnTrain += 1;
						}
					}

					//anterograde particle
					if (static_cast<double>(rand()) / RAND_MAX < incomeFreqPerFrame) // If Anterograde train arrives this instant
					{
						if (brightTrainsLeft > 0) // If there are more bright anterograde trains left
						{
							numberBrightParticle += numberParticleAnt - int(numberParticleAnt*bleachFrequency*travelTime);
							brightTrainsLeft = brightTrainsLeft - 1;
						}
						numberAllParticle += numberParticleAnt;
					}

					if (numberBrightParticle == 0) break; //When there are no more bright particles at tip, stop simulation
				}

				countIncomeTrain += 1;
				numberReturnTrain += tempNumberReturnTrain;
				returnNumberDist[tempNumberReturnTrain]++;
			}

			outFile << k << "\t" << numberReturnTrain / countIncomeTrain; // prints the number of returning trains

			for (int i = 0; i < 10; i++) // For the first 10 returning trains...
			{
				outFile << "\t" << timeReturnTrain[i] / countReturnTrain[i]; //Denominator always 1?, prints average time of returning trains
			}
			outFile << std::endl;


			for (int i = 0; i < 10; i++) // For the first 10 returning trains...
			{
				outFile << "\t" << returnNumberDist[i];//Denominator always 1?, prints distribution of number of returning trains
			}
			outFile << std::endl;


			std::cout << k << std::endl; //prints the value of k

		}
	}
	//for( int i=0; i<100; i++ )std::cout << static_cast<double>(rand())/*RAND_MAX*/ << std::endl;
	
	outFile.close();
	
	std::cout << "Simulation Done" << std::endl;
	std::cin.get(); // Waits to get user input before closing window
}
