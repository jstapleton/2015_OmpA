
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <string>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace std;

//prototypes
char detWeightedAA (double);
int countAAsInSeq (string, char);
vector<int> countAAs_Vector (string);
void printVector (vector<int>); //for debugging
int convertAAtoIndex (char);
double getEnergyOfInsertion (char, double);
double getEnergyOfSequence (string);
double Kelvin(double);
double calcStepSize (double, double, int);

double getSequenceEntropy (vector<int>);

//start list of fixed variables and arrays
//
//list of amino acids
static char AA_Array[20] = {'A','R','N','D','C','Q','E','G','H','I',
							'L','K','M','F','P','S','T','W','Y','V'};

//z-Coords of target residues for mutation


static double z_OmpA[43] =

{				// strand 1
-9.415,  
-4.713, 
0.581, 
5.215, 
10.033, 
				// strand 2
12.552,  
7.184, 
2.558, 
-2.633, 
-7.820, 
-11.438, 
				// strand 3
-9.335, 
-4.231, 
1.164, 
6.528, 
12.107, 
				// strand 4
8.404, 
3.253, 
-2.187, 
-5.440, 
-10.771,
				// strand 5 
-10.952, 
-6.520, 
-1.613, 
2.217, 
7.358, 
11.871,
				// strand 6
9.508, // P
5.628, 
0.770, // G
-4.818, 
-8.923, 
-11.668, 
				// strand 7
-10.812,
-5.602, 
-0.907, 
4.896, 
9.420, 
				// strand 8
11.435, // G
6.920, 
1.568,
-4.141, 
-8.352 // Y
};

//cumulative probability of picking the corresponding AA
static double cuPbPick[21] =
{0,
 0.179876465, 0.179876465, 0.190810611, 0.198628225, 0.198628225,
 0.206162973, 0.206162973, 0.206162973, 0.209617443, 0.328851702,
 0.543269104, 0.545872473, 0.571341673, 0.626695788, 0.626695788,
 0.652938489, 0.718035208, 0.737610534, 0.811937945, 1
};

/*statistically derived energy gained by inserting corresponding AA
in "core region" (from water) of transmembrane beta-barrel*/
static double deltaE[20] = 
{-0.77,   1.42,   1.15,   1.26,     0,   0.64,   1.15,   0.65,   1.21,   -0.99,
 -1.98,   1.28,  -0.6,  -3.0,   0.74,   0.42,   0.44,  -2.08,  -1.33,   -1.48
};

//z-Coord that gives corresponding AA max energy
static double zMid[20] =
{5.96,  11.68,  14.41,  15.09,      0, 11.14,  15.07,   9.86,  7.99,   14.49,
 17.0,  14.27,  6.17,   9.5,  11.13, 5.75,  3.16,  11.49,  9.34,   14.91
};

//best rate of transferring corresponding AA to membrane layer
static double sTrans[20] =
{ 7.06,   1.73,   5.41,     7.22,   0,     11.14,   15.07,   2.86,   13.68,   17.97,
  2.93,   3.73,   30.0,   10.18,   6.53,  3.39,    3.29,   11.49,    9.34,   21.98
};

//frequencies at which corresponding AA is chosen to be mutated to * seqLength
static double absFreq[20] = 
{8.813946794, 0,           0.535773157, 0.383063088, 0,
 0.369202656, 0,           0,           0.169268992, 5.842478707,
10.50645272,  0.127565037, 1.247990838, 2.712351605, 0,
 1.285892373, 3.189739225, 0.959190953, 3.642043143, 9.215040709
};


int main()
{
	// initialize string of surface residues
    string currSeq = "HLLMSFVAAIWFLLINVALMTWILLDSHVLVFKWVLMTHYLMI";


    // redesign #1: IILSVQVAAWVYLLVVTLLMFTFLAGVVILMYHWILMVSILLV
    // redesign #2: FMLMVYIAAVHYLLVWIAASVQLLLTTVLGMISVLLLFWILLI
    // redesign #3: HLLMSFVAAIWFLLINVALMTTVLAIWYILLYQVLLMVSILLV


	string prevSeq;
	string bestSeq = currSeq;

	int currSeqLength = currSeq.length();
	int randLocNumber;

	double energyOfCurrSeq, energyOfPrevSeq, energyOfBestSeq;
	
	double eModOfCurrSeq, eModOfPrevSeq, eModOfBestSeq;
	double penOfBestSeq;
    
	energyOfCurrSeq = getEnergyOfSequence(currSeq);
	
	eModOfCurrSeq = energyOfCurrSeq + 5.06*10E-70*(3.23E35 - getSequenceEntropy(countAAs_Vector(currSeq))) * (3.23E35 - getSequenceEntropy(countAAs_Vector(currSeq)));
	
	eModOfBestSeq = eModOfCurrSeq;

	cout << currSeq << " " << energyOfCurrSeq << " " << eModOfCurrSeq << endl;

	double boltzmannRandom, boltzmannCriteria;

    double numberOfRounds = 500000000;

    for(double j = 1; j<=numberOfRounds; j++)     
        {

        
        srand((unsigned)time(0)); //initiating random seed generator
        int stepNumber = 15;
        
        double T_Celsius_High = 100;
        double T_Celsius_Start = (numberOfRounds - j)/numberOfRounds;
        double T_Celsius_End = .2*(numberOfRounds - j)/numberOfRounds;
        
        double T_Celsius = T_Celsius_Start;

        for(int i = 1; i<= stepNumber; i++)
        {
            prevSeq = currSeq; //copying currSeq to prevSeq: need seq_i-1, where changing seq => seq_i
            randLocNumber = (rand()%(currSeqLength+1)); //between 1 and 49
            double randPropDeterm = (rand()/(static_cast<double>(RAND_MAX)+1.0)); //between 0 and 1
            char dwAA = detWeightedAA(randPropDeterm); //determined AA based on probability weight
            currSeq[randLocNumber] = dwAA; //mutate


		//  case: self-mutation
            if(currSeq == prevSeq)
            {
                i=i-1;
//			    
                continue;
            }

            eModOfCurrSeq = getEnergyOfSequence(currSeq) + 5.06*10E-70*(3.23E35 - getSequenceEntropy(countAAs_Vector(currSeq))) * (3.23E35 - getSequenceEntropy(countAAs_Vector(currSeq)));
            eModOfPrevSeq = getEnergyOfSequence(prevSeq) + 5.06*10E-70*(3.23E35 - getSequenceEntropy(countAAs_Vector(prevSeq))) * (3.23E35 - getSequenceEntropy(countAAs_Vector(prevSeq)));


            //case: energy criteria not met -> Metropolis
            if(eModOfCurrSeq >= eModOfPrevSeq)
            {
                boltzmannRandom = (rand()/(static_cast<double>(RAND_MAX)+1.0));
                boltzmannCriteria = exp(-(eModOfCurrSeq - eModOfPrevSeq)*(T_Celsius/T_Celsius_High));
                if(boltzmannRandom <= boltzmannCriteria)
                {
                    i=i-1;
                    currSeq = prevSeq;
//				    
                    continue;
                }
                
            }

        
            
            if(eModOfBestSeq > eModOfCurrSeq)
            {
                bestSeq = currSeq;
                eModOfBestSeq = eModOfCurrSeq;
                
                energyOfBestSeq = getEnergyOfSequence(bestSeq);
                
                penOfBestSeq = energyOfBestSeq - eModOfBestSeq;

                cout << "Best sequence so far is: " << endl
                    << bestSeq << " with best energy " << energyOfBestSeq << " - penalty " << penOfBestSeq
                    
                    << " = " << eModOfBestSeq << endl << "Round: " << j << ", "
                    
                    << "Cycle: " << i << " with temperature " << T_Celsius << endl;
            }

            T_Celsius = T_Celsius - calcStepSize(T_Celsius_End, T_Celsius_Start, stepNumber);
        }


        
        currSeq = bestSeq;
        
        energyOfCurrSeq = energyOfBestSeq;
    
    }
    
    cout << "Best sequence after " << numberOfRounds << " rounds is: " << endl
         << bestSeq << " with best energy " << energyOfBestSeq << endl;

    
    printVector(countAAs_Vector(bestSeq));
	return(0);
}

char detWeightedAA(double prob)
{
	char output;

	for(int i = 0; i <= 20; i++)
	{
		if(prob >= cuPbPick[i] && prob < cuPbPick[i+1])
		{
			output = AA_Array[i];
			break;
		}
	}

	return output;
}

int countAAsInSeq(string sequence, char amino)
{
	int AA_counter = 0;
	for (int i = 0; i< sequence.length(); i++)
	{
		if(sequence[i] == amino)
			AA_counter++;
	}
	return AA_counter;
}

vector<int> countAAs_Vector(string sequence)
{
	vector<int> vectorOfAAs;
	for (int i = 0; i< 20; i++)
	{
		vectorOfAAs.push_back(countAAsInSeq(sequence, AA_Array[i]));
	}
	return vectorOfAAs;
}

void printVector(vector<int> countAAs_Vector)
{
	for(int i = 0; i< 20; i++)
	{
		cout << "current AA_Array" << " has " << setw(2)
			 << countAAs_Vector[i] << " instances of the letter " << AA_Array[i] << endl << endl;
	}
}

int convertAAtoIndex(char letter)
{
	int indexNum;
	switch(letter)
	{
		case 'A':
			indexNum = 0; break;
		case 'R':
			indexNum = 1; break;
		case 'N':
			indexNum = 2; break;
		case 'D':
			indexNum = 3; break;
		case 'C':
			indexNum = 4; break;
		case 'Q':
			indexNum = 5; break;
		case 'E':
			indexNum = 6; break;
		case 'G':
			indexNum = 7; break;
		case 'H':
			indexNum = 8; break;
		case 'I':
			indexNum = 9; break;
		case 'L':
			indexNum = 10; break;
		case 'K':
			indexNum = 11; break;
		case 'M':
			indexNum = 12; break;
		case 'F':
			indexNum = 13; break;
		case 'P':
			indexNum = 14; break;
		case 'S':
			indexNum = 15; break;
		case 'T':
			indexNum = 16; break;
		case 'W':
			indexNum = 17; break;
		case 'Y':
			indexNum = 18; break;
		case 'V':
			indexNum = 19; break;
		default:
			indexNum = -1; break;
	}
	return indexNum;
}
double getEnergyOfInsertion (char letter, double z)
{
	double energyOfInsertion;
	int iNum = convertAAtoIndex(letter);
	switch(letter)
	{
		case 'A':
		case 'M':
		case 'V':
		case 'L':
		case 'I':
		case 'N':
		case 'D':
		case 'Q':
		case 'E':
		case 'K':
		case 'R':
		case 'H':
		case 'P':
			energyOfInsertion = deltaE[iNum]/(1+(pow((fabs(z) / zMid[iNum]),sTrans[iNum])));
			break;
		case 'Y':
		case 'W':
        case 'F':
        case 'G':
			energyOfInsertion = deltaE[iNum] * (exp(-1 * (pow((fabs(z) - zMid[iNum]),2)) / (2 * pow(sTrans[iNum],2))));
//
//  **** NOTE: when the first set of redesigns were created, this code contained the 
//        following typo that affected the energy calculations:
//			energyOfInsertion = deltaE[iNum] * (exp(-1 * (pow((fabs(z) - zMid[iNum]),2)))) / (2 * pow(sTrans[iNum],2));
//
			break;
		default:
			energyOfInsertion = 0;
			break;
	}

	return energyOfInsertion;
}

double getEnergyOfSequence (string sequence)
{
	double energyOfSequence = 0;
    for(int i = 0; i< sequence.length(); i++)
	{
		energyOfSequence += getEnergyOfInsertion(sequence[i],z_OmpA[i]);
	}
	return energyOfSequence;
}

double Kelvin(double T_Celsius)
{
	double T_Kelvin = T_Celsius + 273.15;
	return T_Kelvin;
}
double calcStepSize (double numFinal, double numInit, int stepNum)
{
	double stepSize = (numInit - numFinal) / (stepNum-1);
	return stepSize;
}



// Calculates sequence complexity: (factorial of number of residues)/(product of factorials of number of each residue)

double getSequenceEntropy (vector<int> countVect)

{
    
    double numerator = 0;
    
    double denominator = 1;
    
    for (int i = 0; i < 20; i++)
    
    {
        
        int num = countVect[i];
        
        numerator += num;
        
        for (int j = num; j > 1; j--)
        
        {
            
            denominator *= j;
        
        }
    
    }

    length = numerator

    for (int i = length - 1; i > 1; i--)
        
    {
            
        numerator *= i;
        
    }
    
    double seqEnt = numerator/denominator;
    
    return seqEnt;

}
