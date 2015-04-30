/*
This script calculates the EzBeta energy of a given string of surface residues.
This was used to generate the EzBeta energy values in Table 1 and Supplemental Table 1.
*/

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <string>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace std;

//prototypes
int convertAAtoIndex (char);
double getEnergyOfInsertion (char, double);
double getEnergyOfSequence (string);

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


int main()
{
	// initialize string of surface residues
    string currSeq = "WTALWLAAGYVVFMYWVLALYLITLGVPFGVYIITLYWGLLVY"; // Wild type
//    string currSeq = "LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL"; // All leucine
//    string currSeq = "IILSVQVAAWVYLLVVTLLMFTFLAGVVILMYHWILMVSILLV"; // Redesign #1
//    string currSeq = "FMLMVYIAAVHYLLVWIAASVQLLLTTVLGMISVLLLFWILLI"; // Redesign #2
//    string currSeq = "HLLMSFVAAIWFLLINVALMTTVLAIWYILLYQVLLMVSILLV"; // Redesign #3
//    string currSeq = "FLLMIQIAAIWHLLVVYALMVVVLLIQFMLSHVVMLTHWVLAY"; // Redesign #4
	
//    string currSeq = "WTALWQIAAIVVFMYWVLALYVVLLIVPFGVYIIMLTHWVLAY"; // H1346
//    string currSeq = "WTALWLAAGYVVLLVVYALMVLITLGVPFGVYIIMLTHWVLAY"; // H1256
//    string currSeq = "FLLMILAAGYVVFMYWVLALYLITLGVFMLSHIIMLTHWVLAY"; // H2345
//    string currSeq = "FLLMIQIAAIVVFMYWYALMVLITLGVPFGVYIITLYWWVLAY"; // H3567
//    string currSeq = "FLLMILAAGYVVFMYWYALMVVVLLIVPFGVYIIMLTHGLLVY"; // H2368
//    string currSeq = "WTALWQIAAIVVLLVVYALMVVVLLIVPFGVYIITLYWGLLVY"; // H1678
//    string currSeq = "WTALWLAAGYVVFMYWYALMVVVLLIVFMLSHIITLYWWVLAY"; // H1237
//    string currSeq = "FLLMILAAGYVVLLVVYALMVLITLGVFMLSHIITLYWGLLVY"; // H2578
//    string currSeq = "FLLMIQIAAIVVFMYWVLALYVVLLIVFMLSHIITLYWGLLVY"; // H3478
//    string currSeq = "WTALWLAAGYVVLLVVVLALYVVLLIVFMLSHIIMLTHGLLVY"; // H1248
//    string currSeq = "FLLMILAAGYVVLLVVVLALYVVLLIVPFGVYIITLYWWVLAY"; // H2467
//    string currSeq = "WTALWQIAAIVVLLVVVLALYLITLGVFMLSHIITLYWWVLAY"; // H1457
//    string currSeq = "FLLMIQIAAIVVLLVVVLALYLITLGVPFGVYIIMLTHGLLVY"; // H4568
//    string currSeq = "WTALWQIAAIVVFMYWYALMVLITLGVFMLSHIIMLTHGLLVY"; // H1358

//    string currSeq = "WTALWQIAAIVVLLVVYALMVVVLLIVPFGVYIIMLTHWVLAY"; // H16
//    string currSeq = "FLLMIQIAAIVVFMYWYALMVVVLLIVPFGVYIIMLTHWVLAY"; // H36
//    string currSeq = "WTALWQIAAIVVFMYWYALMVVVLLIVPFGVYIIMLTHWVLAY"; // H136
//    string currSeq = "WTALWQIAAIVVLLVVVLALYVVLLIVPFGVYIIMLTHWVLAY"; // H146
//    string currSeq = "FLLMIQIAAIVVFMYWVLALYVVLLIVPFGVYIIMLTHWVLAY"; // H346

//    string currSeq = "FLLMILIAGIVVLLVVYALMVVVLLGVPMGSYIIMLTWLLLLY"; // OR4cons
//    string currSeq = "FLLMILIAAIVVLLVVYALMVVVLLGVPMGSYIIMLTWLLLLY"; // OR4cons_G41A
//    string currSeq = "FLLMILIAGIVVLLVVYALMVVVLLIVPMGSYIIMLTWLLLLY"; // OR4cons_G99I
//    string currSeq = "FLLMILIAGIVVLLVVYALMVVVLLGVFMGSYIIMLTWLLLLY"; // OR4cons_P121F
//    string currSeq = "FLLMILIAGIVVLLVVYALMVVVLLGVPMLSYIIMLTWLLLLY"; // OR4cons_G125L
//    string currSeq = "FLLMILIAGIVVLLVVYALMVVVLLGVPMGSHIIMLTWLLLLY"; // OR4cons_Y129H
//    string currSeq = "FLLMILIAGIVVLLVVYALMVVVLLGVPMGSLIIMLTWLLLLY"; // OR4cons_Y129L



	double energyOfCurrSeq;
	    
	energyOfCurrSeq = getEnergyOfSequence(currSeq);
			
	cout << currSeq << " " << energyOfCurrSeq << endl;

	return(0);
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
