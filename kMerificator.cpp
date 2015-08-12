//============================================================================
// Name        : kMerificator.cpp
// Author      : 
// Version     :
// Copyright   : 
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <vector>

#include "Utilities.h"
#include "Validator.h"

using namespace std;

int main(int argc, char *argv[]) {

	vector<string> arguments (argv + 1, argv + argc + !argc);

	int cortex_height = 26;
	int cortex_width = 140;
	int kMer_size = 31;
	int threads = 40;

	for(unsigned int i = 2; i < arguments.size(); i++)
	{
		if(arguments.at(i) == "--cortex_height")
		{
			cortex_height = Utilities::StrtoI(arguments.at(i+1));
		}

		if(arguments.at(i) == "--cortex_width")
		{
			cortex_width = Utilities::StrtoI(arguments.at(i+1));
		}

		if(arguments.at(i) == "--kmer")
		{
			kMer_size = Utilities::StrtoI(arguments.at(i+1));
		}
	}

	string vcfFile = arguments.at(2);
	string referenceGenome = arguments.at(3);
	string deBruijnGraph = arguments.at(4);
	string outputDirectory = arguments.at(5);

	validateCompleteVCF<1, 31>(vcfFile, referenceGenome, deBruijnGraph, kMer_size, cortex_height, cortex_width, outputDirectory, threads);

	return 0;
}
