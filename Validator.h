/*
 * Validator.h
 *
 *  Created on: 12 Aug 2015
 *      Author: dilthey
 */

#ifndef VALIDATOR_H_
#define VALIDATOR_H_

#include <vector>
#include <string>
#include <set>

#include <iostream>
#include <fstream>
#include <omp.h>
#include <assert.h>

#include "hash/deBruijn/DeBruijnGraph.h"

#include "Utilities.h"

typedef std::vector< std::vector<std::string> > diploidGenomeString;

diploidGenomeString VCF2GenomeString(std::string chromosomeID, int positionStart, int positionStop, std::string VCFpath, std::string referenceGenomePath, std::vector<std::vector<int> >& ret_graph_referencePositions, bool ignoreVCF = false, bool onlyPASS = false);
diploidGenomeString compressGenomeString(diploidGenomeString gS);
std::vector<double> kMer_PP(size_t observedCoverage, size_t upperBound, double coverage);

template<int m, int k, int colours>
void evaluate_dGS(diploidGenomeString& gS, diploidGenomeString& gS_unresolved, const std::set<std::string>& kMers_reference, DeBruijnGraph<m, k, colours>* graph, std::set<std::string>* kMers_in_dGS, std::set<std::string>* kMers_in_dGS_in_sample, std::map<std::string, double>* kMers_in_dGS_optimality, std::string nameForSummary, ofstream& summaryFileStream, std::string pathForSpatialSummary, std::vector<std::vector<int> > chromotypes_referencePositions);

template<int m, int k>
void validateCompleteVCF(std::string VCFfile, std::string referenceGenome, std::string deBruijnGraph, int cortex_height, int cortex_width, std::string outputDirectory, int threads, bool onlyPASS);

template<int m, int k, int colours>
std::pair<diploidGenomeString, diploidGenomeString> greedilyResolveDiploidKMerString(diploidGenomeString& original_gS, DeBruijnGraph<m, k, colours>* graph);


template<int m, int k>
void validateCompleteVCF(std::string VCFfile, std::string referenceGenome, std::string deBruijnGraph, int cortex_height, int cortex_width, std::string outputDirectory, int threads, bool onlyPASS)
{
	omp_set_num_threads(threads);
	assert(omp_get_num_threads() == threads);
	int kMer_size = k;

	std::string configFileOutputPath= outputDirectory + "/validationDetails.txt";

	ofstream configOutputStream;

	configOutputStream.open(configFileOutputPath.c_str());
	if(! configOutputStream.is_open())
	{
		throw std::runtime_error("validateCompleteVCF(..): Want top open config summary file for writing, but can't! Path:\n"+configFileOutputPath);
	}
	configOutputStream << "k = " << kMer_size << "\n";
	configOutputStream << "deBruijnGraph" << " = " << deBruijnGraph << "\n";
	configOutputStream << "VCF" << " = " << VCFfile << "\n";
	configOutputStream << "referenceGenome" << " = " << referenceGenome << "\n";
	configOutputStream << "onlyPASS" << " = " << onlyPASS << "\n";
	configOutputStream << "threads" << " = " << threads << "\n";

	configOutputStream.close();

	std::cout << Utilities::timestamp() << "Allocate Cortex graph object with height = " << cortex_height << ", width = " << cortex_width << " ...\n" << std::flush;
	DeBruijnGraph<m, k, 1> myGraph(cortex_height, cortex_width);
	std::cout << Utilities::timestamp() << "Cortex graph object allocated, loading binary...\n" << std::flush;

	myGraph.loadMultiColourBinary(deBruijnGraph);
	std::cout << Utilities::timestamp() << "\tdone\n" << std::flush;
	std::cout << "\tTotal coverage: " << myGraph.totalCoverage() << "\n";

	std::vector<std::string> chromosomes;
	for(int i = 1; i <= 22; i++)
	{
		chromosomes.push_back(Utilities::ItoStr(i));
	}
	chromosomes.push_back("X");
	chromosomes.push_back("Y");


	// todo remove
	// chromosomes.clear();
	// chromosomes.push_back("22");

	// chromosomes.clear();
	// chromosomes.push_back("6");

	std::string summaryFilePath = outputDirectory + "/summaryPerChromosome.txt";
	ofstream summaryFileStream;
	summaryFileStream.open(summaryFilePath.c_str());
	assert(summaryFileStream.is_open());

	std::vector<std::string> headerFields;
	headerFields.push_back("Method");
	headerFields.push_back("Total characters");
	headerFields.push_back("Total non-gap characters");
	headerFields.push_back("# kMers");
	headerFields.push_back("# kMers invalid");
	headerFields.push_back("# kMers present");
	headerFields.push_back("Unweighted optimality");
	headerFields.push_back("Coverage-weighted optimality");

	summaryFileStream << Utilities::join(headerFields, "\t") << "\n";

	for(unsigned int chromosomeI = 0; chromosomeI < chromosomes.size(); chromosomeI++)
	{
		std::string chromosome = chromosomes.at(chromosomeI);

		std::cout << Utilities::timestamp() << " " << "Chromosome " << chromosome << ", convert VCF to diploidGS...\n" << std::flush;
		std::vector<std::vector<int> > VCF_chromotypes_referencePositions;
		diploidGenomeString VCF_chromotypes = VCF2GenomeString(chromosome, -1, -1, VCFfile, referenceGenome, VCF_chromotypes_referencePositions, false, onlyPASS);
		std::cout << Utilities::timestamp() << " " <<  "\n\tdone\n" << std::flush;

		std::cout << Utilities::timestamp() << " " <<  "Chromosome " << chromosome << ", compress VCF diploidGS...\n" << std::flush;
		diploidGenomeString VCF_chromotypes_compressed = compressGenomeString(VCF_chromotypes);
		std::cout << Utilities::timestamp() << " " <<  "\n\tdone\n" << std::flush;

		std::cout << Utilities::timestamp() << " " <<  "Chromosome " << chromosome << ", resolve chromotypes from VCF..\n" << std::flush;
		diploidGenomeString VCF_chromotypes_resolved = greedilyResolveDiploidKMerString<m, k, 1>(VCF_chromotypes_compressed, &myGraph).second;
		std::cout << Utilities::timestamp() << " " <<  "\n\tdone\n" << std::flush;

		std::set<std::string> set_kMers_reference;
		std::set<std::string> set_kMers_VCF;
		std::set<std::string> set_kMers_VCF_present;
		std::map<std::string, double> set_kMers_VCF_optimalities;

		std::cout << Utilities::timestamp() << " " << "Chromosome " << chromosome << ", evaluate VCF chromotypes\n" << std::flush;
		evaluate_dGS(VCF_chromotypes_resolved, VCF_chromotypes_compressed, set_kMers_reference, &myGraph, 0, 0, 0, "VCF"+chromosome, summaryFileStream, outputDirectory, VCF_chromotypes_referencePositions);

		std::cout << Utilities::timestamp() << "\n\tdone" <<  "\n\n" << std::flush;

		std::cout << Utilities::timestamp() << " " << "Chromosome " << chromosome << ", convert VCF (reference-only) to diploidGS...\n" << std::flush;

		std::vector<std::vector<int> > reference_chromotypes_referencePositions;
		std::cout << "Get non-modified diploidGS for reference genome...\n" << std::flush;
		diploidGenomeString gS_ref = VCF2GenomeString(chromosome, -1, -1, VCFfile, referenceGenome, reference_chromotypes_referencePositions, true);
		std::cout << "\tdone\n" << std::flush;

		std::cout << Utilities::timestamp() << " " <<  "Chromosome " << chromosome << ", compress VCF (reference-only) diploidGS...\n" << std::flush;
		diploidGenomeString gS_ref_compressed = compressGenomeString(gS_ref);
		std::cout << Utilities::timestamp() << " " <<  "\n\tdone\n" << std::flush;

		set_kMers_reference.clear();
		set_kMers_VCF.clear();
		set_kMers_VCF_present.clear();
		set_kMers_VCF_optimalities.clear();

		std::cout << Utilities::timestamp() << " " << "Chromosome " << chromosome << ", evaluate VCF (reference-only) chromotypes\n" << std::flush;
		evaluate_dGS(gS_ref_compressed, gS_ref_compressed, set_kMers_reference, &myGraph, 0, 0, 0, "toReference"+chromosome, summaryFileStream, outputDirectory, reference_chromotypes_referencePositions);

		std::cout << Utilities::timestamp() << "\n\tdone" <<  "\n\n" << std::flush;
	}
}

template<int m, int k, int colours>
void evaluate_dGS(diploidGenomeString& gS, diploidGenomeString& gS_unresolved, const std::set<std::string>& kMers_reference, DeBruijnGraph<m, k, colours>* graph, std::set<std::string>* kMers_in_dGS, std::set<std::string>* kMers_in_dGS_in_sample, std::map<std::string, double>* kMers_in_dGS_optimality, std::string nameForSummary, ofstream& summaryFileStream, std::string pathForSpatialSummary, std::vector<std::vector<int> > chromotypes_referencePositions)
{
	// We only want to deal with completely resolved genomestrings up to level k
	// That is, we want genomestrings which
	// - have all subsequent homozygous stretches connected
	// - have at least k homozygous characters between any two heterozygous positions
	// - so that if we connect them together, all k-mers are totally determined

	size_t totalLevels = 0;

	// std::cout << Utilities::timestamp()  << "evaluate_dGS " << "A" << "\n" << std::flush; // todo remove

	for(unsigned int l = 0; l < gS.size(); l++)
	{
		if(gS.at(l).size() == 1)
		{
			if(!((l == 0) || (gS.at(l-1).size() != 1)))
			{
				std::cerr << "!((l == 0) || (gS.at(l-1).size() != 1)), l = " << l << ", gS.at(l-1).size() = " << gS.at(l-1).size() << ", gS.size()= " << gS.size() << "\n" << std::flush;
			}
			if(!((l == (gS.size() - 1)) || (gS.at(l+1).size() != 1)))
			{
				std::cerr << "!((l == (gS.size() - 1)) || (gS.at(l+1).size() != 1)), l = " << l << ", gS.at(l+1).size() = " << gS.at(l+1).size() << ", gS.size()= " << gS.size() << "\n" << std::flush;
			}
			assert((l == 0) || (gS.at(l-1).size() != 1));
			assert((l == (gS.size() - 1)) || (gS.at(l+1).size() != 1));
		}


		unsigned int requiredSpace = gS.at(l).at(0).size();

		totalLevels += requiredSpace;
	}

	// std::cout << Utilities::timestamp()  << "evaluate_dGS " << "B" << "\n" << std::flush; // todo remove

	std::vector<std::vector<std::string> > uncompressed_chromotypes;
	std::vector<int> uncompressed_chromotypes_nGaps;

	uncompressed_chromotypes.reserve(totalLevels);
	uncompressed_chromotypes_nGaps.reserve(totalLevels);

	for(unsigned int outerI = 0; outerI < gS.size(); outerI++)
	{
		std::vector<std::string>& segment = gS.at(outerI);
		assert((segment.size() == 1) || (segment.size() == 2));
		if(segment.size() == 2)
		{
			assert(segment.at(0).size() == segment.at(1).size());
		}

		unsigned int requiredSpace = segment.at(0).size();

		for(unsigned int posInSegment = 0; posInSegment < requiredSpace; posInSegment++)
		{
			std::vector<std::string> posVector;
			std::vector<double> posCoverageVector;
			std::vector<double> posSupportVector;
			std::vector<double> posInsertionsVector;

			if(segment.size() == 2)
			{
				posVector.push_back(segment.at(0).substr(posInSegment, 1));
				posVector.push_back(segment.at(1).substr(posInSegment, 1));
				posCoverageVector.push_back(0);
				posCoverageVector.push_back(0);
				posSupportVector.push_back(0);
				posSupportVector.push_back(0);
				posInsertionsVector.push_back(0);
				posInsertionsVector.push_back(0);
			}
			else
			{
				posVector.push_back(segment.at(0).substr(posInSegment, 1));
				posCoverageVector.push_back(0);
				posSupportVector.push_back(0);
				posInsertionsVector.push_back(0);
			}

			uncompressed_chromotypes.push_back(posVector);

			int numGaps = 0;
			if(segment.size() == 2)
			{
				if(segment.at(0).substr(posInSegment, 1) == "_")
				{
					numGaps++;
				}
				if(segment.at(1).substr(posInSegment, 1) == "_")
				{
					numGaps++;
				}
			}
			else
			{
				if(segment.at(0).substr(posInSegment, 1) == "_")
				{
					numGaps = 2;
				}
			}
			uncompressed_chromotypes_nGaps.push_back(numGaps);
		}
	}
	assert(uncompressed_chromotypes.size() == uncompressed_chromotypes_nGaps.size());

	// std::cout << Utilities::timestamp()  << "evaluate_dGS " << "C" << "\n" << std::flush; // todo remove

	std::vector<int> uncompressed_chromotypes_referencePositions;
	uncompressed_chromotypes_referencePositions.reserve(totalLevels);

	for(unsigned int outerI = 0; outerI < chromotypes_referencePositions.size(); outerI++)
	{
		for(unsigned int posInSegment = 0; posInSegment < chromotypes_referencePositions.at(outerI).size(); posInSegment++)
		{
			uncompressed_chromotypes_referencePositions.push_back(chromotypes_referencePositions.at(outerI).at(posInSegment));
		}
	}

	if(!(uncompressed_chromotypes.size() == uncompressed_chromotypes_referencePositions.size()))
	{
		std::cerr << "! (uncompressed_chromotypes.size() == uncompressed_chromotypes_referencePositions.size())" << "\n";
		std::cerr << "uncompressed_chromotypes.size()" << ": " << uncompressed_chromotypes.size() << "\n";
		std::cerr << "uncompressed_chromotypes_referencePositions.size()" << ": " << uncompressed_chromotypes_referencePositions.size() << "\n\n" << std::flush;
		std::cerr << "gS.size()" << ": " << gS.size() << "\n";
		std::cerr << "chromotypes_referencePositions.size()" << ": " << chromotypes_referencePositions.size() << "\n\n";

		std::cerr << std::flush;
	}

	assert(uncompressed_chromotypes.size() == uncompressed_chromotypes_referencePositions.size());

	std::cout << Utilities::timestamp()  << "evaluate_dGS " << "D" << "\n" << std::flush; // todo remove

	std::vector<int> uncompressed_chromotypes_diploid;
	std::vector<int> uncompressed_chromotypes_losePhasing;
	uncompressed_chromotypes_diploid.reserve(totalLevels);
	uncompressed_chromotypes_losePhasing.reserve(totalLevels);

	for(unsigned int outerI = 0; outerI < gS_unresolved.size(); outerI++)
	{
		std::vector<std::string>& segment = gS_unresolved.at(outerI);
		assert((segment.size() == 1) || (segment.size() == 2));
		if(segment.size() == 2)
		{
			assert(segment.at(0).size() == segment.at(1).size());
		}

		for(unsigned int posInString = 0; posInString < segment.at(0).length(); posInString++)
		{
			uncompressed_chromotypes_diploid.push_back(segment.size() - 1);
			if((posInString == 0) && (segment.size() == 2))
			{
				uncompressed_chromotypes_losePhasing.push_back(1);
			}
			else
			{
				uncompressed_chromotypes_losePhasing.push_back(0);
			}
		}
	}



	if(!(uncompressed_chromotypes.size() == uncompressed_chromotypes_diploid.size()))
	{
		std::cerr << "! uncompressed_chromotypes.size() == uncompressed_chromotypes_diploid.size()" << "\n";
		std::cerr << "uncompressed_chromotypes.size()" << ": " << uncompressed_chromotypes.size() << "\n";
		std::cerr << "uncompressed_chromotypes_diploid.size()" << ": " << uncompressed_chromotypes_diploid.size() << "\n\n" << std::flush;
		std::cerr << "gS.size()" << ": " << gS.size() << "\n";
		std::cerr << "uncompressed_chromotypes_referencePositions.size()" << ": " << uncompressed_chromotypes_referencePositions.size() << "\n\n";

		std::cerr << std::flush;
	}

	// std::cout << Utilities::timestamp()  << "evaluate_dGS " << "E" << "\n" << std::flush; // todo remove

	assert(uncompressed_chromotypes.size() == uncompressed_chromotypes_diploid.size());
	assert(uncompressed_chromotypes.size() == uncompressed_chromotypes_losePhasing.size());

	std::string s1; std::string s2;
	s1.reserve(totalLevels);
	s2.reserve(totalLevels);
	for(unsigned int l = 0; l < gS.size(); l++)
	{
		if(gS.at(l).size() == 1)
		{
			s1.append(gS.at(l).at(0));
			s2.append(gS.at(l).at(0));

			unsigned int nonGap_characters = 0;
			for(unsigned int i = 0; i < gS.at(l).at(0).size(); i++)
			{
				if(gS.at(l).at(0).at(i) != '_')
				{
					nonGap_characters++;
				}
			}

			if((l != 0) && (l != (gS.size() - 1)))
			{
				assert(nonGap_characters >= k);
			}
		}
		else
		{
			assert(gS.at(l).size() == 2);
			assert((l == 0) || (gS.at(l-1).size() == 1));
			assert((l == (gS.size() - 1)) || (gS.at(l+1).size() == 1));

			s1.append(gS.at(l).at(0));
			s2.append(gS.at(l).at(1));
		}
	}

	std::cout << Utilities::timestamp()  << "evaluate_dGS " << "F" << "\n" << std::flush; // todo remove


	// now count k-Mers
	size_t total_characters = 0;
	size_t total_nonGap_characters = 0;

	std::map<std::string, int> kMer_counts;
	size_t kMersCannotEvaluate = 0;

	std::vector<int> inducedkMers;
	std::vector<int> inducedkMers_OK;
	std::vector<int> inducedkMers_invalid;
	std::vector<int> inducedkMers_inReference;

	std::string filenameForSpatialSummary = pathForSpatialSummary + "/spatialSummary_" + nameForSummary + ".txt";
	ofstream spatialSummaryStream;
	spatialSummaryStream.open(filenameForSpatialSummary.c_str());
	assert(spatialSummaryStream.is_open());

	spatialSummaryStream <<
			"Level" << "\t" <<
			"ReferenceCoordinate" << "\t" <<
			"kMersInduced" << "\t" <<
			"kMersInvalid" << "\t" <<
			"kMersOK" << "\t" <<
			"InducedkMersInReference" << "\t" <<
			"DiploidChromotype" << "\t" <<
			"ChromotypeLostPhase" << "\t" <<
			"GapsAtLevel" << "\n";

	auto partitionStringIntokMers_overGaps = [](std::string s, int kMerSize) -> std::vector<std::string> {
		std::vector<std::string> forReturn;

		forReturn.reserve((int)s.length() - kMerSize);

		for(int i = 0; i <= ((int)s.length() - kMerSize); i++)
		{
			std::string C = s.substr(i, 1);
			if(C == "_")
			{
				forReturn.push_back(std::string(""));
			}
			else
			{
				int foundMoreNonGaps = 0;
				int currentStop = i;
				while((currentStop < ((int)s.length()-1)) && (foundMoreNonGaps < (kMerSize-1)))
				{
					currentStop++;
					std::string justAddedC = s.substr(currentStop,1);
					if(justAddedC != "_")
					{
						foundMoreNonGaps++;
					}
				}

				if(foundMoreNonGaps == (kMerSize-1))
				{
					std::string kMer = C + s.substr(i+1, (currentStop-(i+1)+1));
					kMer.erase(std::remove_if(kMer.begin(),kMer.end(), [&](char c){return ((c == '_') ? true : false);}), kMer.end());
					assert((int)kMer.length() == kMerSize);

					forReturn.push_back(kMer);
				}
				else
				{
					forReturn.push_back(std::string(""));
				}

			}
		}

		assert(forReturn.size() == (s.length() - kMerSize + 1));

		return forReturn;
	};

	std::cout << Utilities::timestamp()  << "evaluate_dGS " << "G" << "\n" << std::flush; // todo remove

	std::string testString = "AHCCHDHRHJSKSJSKKSKLSSOMCSNCSDKDKDKDKDKCMCMCKWL";
	std::vector<std::string> kMers_1 = partitionStringIntokMers(testString, 5);
	std::vector<std::string> kMers_2 = partitionStringIntokMers_overGaps(testString, 5);
	assert(kMers_1.size() == kMers_2.size());
	for(unsigned int i = 0; i < kMers_1.size(); i++)
	{
		assert(kMers_1.at(i) == kMers_2.at(i));
	}

	auto processString = [&](std::string s) {
		for(size_t i = 0; i < s.length(); i++)
		{
			total_characters++;
			if(s.at(i) != '_')
			{
				total_nonGap_characters++;
			}
		}
		//s.erase(std::remove_if(s.begin(),s.end(), [&](char c){return ((c == '_') ? true : false);}), s.end());

		std::vector<std::string> kMers = partitionStringIntokMers_overGaps(s, k);

		for(size_t i = 0; i < kMers.size(); i++)
		{
			std::string kMer = kMers.at(i);

			bool containsGap = (kMer.length() == 0);
			bool kMerOK = true;
			for(unsigned int kI = 0; kI < kMer.length(); kI++)
			{
				char kMerC = kMer.at(kI);
				assert(kMerC != '_');
				if(!((kMerC == 'A') || (kMerC == 'C') || (kMerC == 'G') || (kMerC == 'T')))
				{
					kMerOK = false;
					break;
				}
			}

			if(containsGap)
			{
				// nothing
			}
			else
			{
				inducedkMers.at(i)++;

				if(kMerOK)
				{
					if(kMer_counts.count(kMer) == 0)
						kMer_counts[kMer] = 0;
					kMer_counts[kMer]++;

					if(graph->kMerinGraph(kMer))
					{
						inducedkMers_OK.at(i)++;
					}
				}
				else
				{
					kMersCannotEvaluate++;
					inducedkMers_invalid.at(i)++;
				}

				if(kMers_reference.size() == 0)
				{

				}
				else
				{
					if(kMers_reference.count(kMer))
					{
						inducedkMers_inReference.at(i)++;
					}
				}
			}
		}
	};

	assert(s1.length() > k);
	inducedkMers.resize(s1.length()-k+1,0);
	inducedkMers_invalid.resize(s1.length()-k+1,0);
	inducedkMers_OK.resize(s1.length()-k+1,0);
	inducedkMers_inReference.resize(s1.length()-k+1,0);
	assert(s1.length() == s2.length());
	processString(s1);
	processString(s2);

	for(unsigned int level = 0; level < inducedkMers.size(); level++)
	{
		spatialSummaryStream <<
				level << "\t" <<
				uncompressed_chromotypes_referencePositions.at(level) << "\t" <<
				inducedkMers.at(level) << "\t" <<
				inducedkMers_invalid.at(level) << "\t" <<
				inducedkMers_OK.at(level) << "\t" <<
				inducedkMers_inReference.at(level) << "\t" <<
				uncompressed_chromotypes_diploid.at(level) << "\t" <<
				uncompressed_chromotypes_losePhasing.at(level) << "\t" <<
				uncompressed_chromotypes_nGaps.at(level) << "\n";
	}

	spatialSummaryStream.close();



	size_t total_kMers_present = 0;
	size_t maximumkMerCount = 0;
	size_t maximumAllowedkMerCount = 100;

	std::cout << Utilities::timestamp()  << "evaluate_dGS " << "H" << "\n" << std::flush; // todo remove

	for(std::map<std::string, int>::iterator kMerIt = kMer_counts.begin(); kMerIt != kMer_counts.end(); kMerIt++)
	{
		std::string kMer = kMerIt->first;
		int kMerCount = kMerIt->second;

		if(kMers_in_dGS != 0)
		{
			kMers_in_dGS->insert(kMer);
		}

		if(graph->kMerinGraph(kMer))
		{
			total_kMers_present++;
			if(kMers_in_dGS_in_sample != 0)
			{
				kMers_in_dGS_in_sample->insert(kMer);
			}
		}

		if(kMerCount > maximumkMerCount)
			maximumkMerCount = kMerCount;
	}

	double haploidGenomeSize = 3.2e9;
	size_t graph_totalCoverage = graph->totalCoverage();
	double averageCoverage = double(graph_totalCoverage)/(2.0 * haploidGenomeSize);
	std::cout << "Estimate coverage posterior probabilties for max count of " << maximumkMerCount << ", assuming WG average coverage of " << averageCoverage << " [total coverage on graph: " <<  graph_totalCoverage << "]\n";

	if(maximumkMerCount > maximumAllowedkMerCount)
	{
		maximumkMerCount = maximumAllowedkMerCount;
		std::cout << "Restrict maximum kMer count (underlying) to " << maximumAllowedkMerCount << " - all kMers with higher underlying numbers will be ignored!" << "\n" << std::flush;
	}

	if(averageCoverage == 0)
	{
		averageCoverage = 1;
	}
	assert(averageCoverage != 0);
	std::map<size_t, std::vector<double> > _pp_cache;
	auto getPPs = [&](size_t observedCoverage) -> std::vector<double> {
		if(_pp_cache.count(observedCoverage) > 0)
		{
			return _pp_cache.at(observedCoverage);
		}
		else
		{
			std::vector<double> forReturn = kMer_PP(observedCoverage, maximumkMerCount, averageCoverage);
			if(observedCoverage < (2 * averageCoverage))
			{
				_pp_cache[observedCoverage] = forReturn;
			}
			return forReturn;
		}
	};

	std::cout << Utilities::timestamp()  << "evaluate_dGS " << "I" << "\n" << std::flush; // todo remove

	double expected_kMer_presence = 0;
	double expected_kMer_presence_sum = 0;
	for(std::map<std::string, int>::iterator kMerIt = kMer_counts.begin(); kMerIt != kMer_counts.end(); kMerIt++)
	{
		std::string kMer = kMerIt->first;
		int kMerCount = kMerIt->second;

		int kMerCoverage = graph->kMer_getCoverage(kMer);
		std::vector<double> pps;
		if(kMerCount <= maximumAllowedkMerCount)
		{
			pps = getPPs(kMerCoverage);

			double expected_missing = 0;
			for(int i = 0; i < kMerCount; i++)
			{
				int e_missing = kMerCount - i;
				double e_p = pps.at(i);
				expected_missing += e_missing * e_p;
			}

			assert(expected_missing >= 0);
			assert(expected_missing <= kMerCount);

			double this_kMer_presence = (kMerCount - expected_missing);
			expected_kMer_presence += this_kMer_presence;
			expected_kMer_presence_sum += kMerCount;

			if(kMers_in_dGS_optimality != 0)
			{
				(*kMers_in_dGS_optimality)[kMer] = (double)this_kMer_presence / (double)kMerCount;
			}
		}
	}

	assert(expected_kMer_presence >= 0);
	assert(expected_kMer_presence <= expected_kMer_presence_sum);
	assert(expected_kMer_presence_sum > 0);

	double unweighted_optimality = (kMer_counts.size() != 0) ? (double(total_kMers_present)/double(kMer_counts.size())) : 1;
	double weighted_probabilistic_optimality = expected_kMer_presence / expected_kMer_presence_sum;

	std::cout << "evaluate_dGS(..):\n";
	std::cout << "\tLength of gS: " << gS.size() << "\n";
	std::cout << "\tTotal characters: " << total_characters << ", of which non-gap: " << total_nonGap_characters << "\n";
	std::cout << "\tTotal kMers: " << kMer_counts.size() << "\n";
	std::cout << "\tTotal invalid kMers (driven by non-A/C/G/T): " << kMersCannotEvaluate << "\n";

	std::cout << "\tTotal kMers present: " << total_kMers_present << "\n";
	std::cout << "\tUnweighted optimality: " << unweighted_optimality << "\n";
	std::cout << "\tWeighted probabilistic optimality: " << weighted_probabilistic_optimality << "\n";

	std::vector<std::string> fieldsForPrinting;

	fieldsForPrinting.push_back(nameForSummary);
	fieldsForPrinting.push_back(Utilities::ItoStr(total_characters));
	fieldsForPrinting.push_back(Utilities::ItoStr(total_nonGap_characters));
	fieldsForPrinting.push_back(Utilities::ItoStr(kMer_counts.size()));
	fieldsForPrinting.push_back(Utilities::ItoStr(kMersCannotEvaluate));
	fieldsForPrinting.push_back(Utilities::ItoStr(total_kMers_present));
	fieldsForPrinting.push_back(Utilities::DtoStr(unweighted_optimality));
	fieldsForPrinting.push_back(Utilities::DtoStr(weighted_probabilistic_optimality));

	summaryFileStream << Utilities::join(fieldsForPrinting, "\t") << "\n";
}


template<int m, int k, int colours>
std::pair<diploidGenomeString, diploidGenomeString> greedilyResolveDiploidKMerString(diploidGenomeString& original_gS, DeBruijnGraph<m, k, colours>* graph)
{
	int threads = 40;
	omp_set_num_threads(threads);
	unsigned int _compartmentI;

	size_t gS_length_before = original_gS.size();

	std::vector<std::pair<std::string, std::string> > openPairs;

	std::pair<std::string, std::string> resolved;
	std::pair<std::string, std::string> resolved_noGaps;

	diploidGenomeString resolved_gS;

	int addedHomozygousSeq = 0;

	auto evaluateSequence = [&](std::string seq) -> double {
		seq.erase(std::remove_if(seq.begin(),seq.end(), [&](char c){return ((c == '_') ? true : false);}), seq.end());
		std::vector<std::string> kMers = partitionStringIntokMers(seq, k);
		std::map<std::string, int> kMers_table;
		for(unsigned int i = 0; i < kMers.size(); i++)
		{
			std::string kMer = kMers.at(i);
			if(kMers_table.count(kMer) == 0)
			{
				kMers_table[kMer] = 0;
			}
			kMers_table[kMer]++;
		}

		double optimality_sum = 0;
		double optimality_score = 0;

		for(std::map<std::string, int>::iterator kMerIt = kMers_table.begin(); kMerIt != kMers_table.end(); kMerIt++)
		{
			std::string kMerSeq = kMerIt->first;

			bool kMerOK = true;
			for(unsigned int kI = 0; kI < kMerSeq.length(); kI++)
			{
				char kMerC = kMerSeq.at(kI);
				if(!((kMerC == 'A') || (kMerC == 'C') || (kMerC == 'G') || (kMerC == 'T')))
				{
					kMerOK = false;
					break;
				}
			}

			if(kMerOK)
			{
				int count = kMerIt->second;

				bool present = graph->kMerinGraph(kMerSeq);

				if(present)
					optimality_score += count;

				optimality_sum += count;
			}
		}

		if(optimality_sum == 0)
		{
			return 1;
		}
		else
		{
			return optimality_score/optimality_sum;
		}
	};

	auto evaluateSequence2 = [&](std::string seq, int& optimality_sum, int& optimality_score) -> void {
		seq.erase(std::remove_if(seq.begin(),seq.end(), [&](char c){return ((c == '_') ? true : false);}), seq.end());
		std::vector<std::string> kMers = partitionStringIntokMers(seq, k);
		std::map<std::string, int> kMers_table;
		for(unsigned int i = 0; i < kMers.size(); i++)
		{
			std::string kMer = kMers.at(i);
			if(kMers_table.count(kMer) == 0)
			{
				kMers_table[kMer] = 0;
			}
			kMers_table[kMer]++;
		}

		optimality_sum = 0;
		optimality_score = 0;

		for(std::map<std::string, int>::iterator kMerIt = kMers_table.begin(); kMerIt != kMers_table.end(); kMerIt++)
		{
			std::string kMerSeq = kMerIt->first;

			bool kMerOK = true;
			for(unsigned int kI = 0; kI < kMerSeq.length(); kI++)
			{
				char kMerC = kMerSeq.at(kI);
				if(!((kMerC == 'A') || (kMerC == 'C') || (kMerC == 'G') || (kMerC == 'T')))
				{
					kMerOK = false;
					break;
				}
			}

			if(kMerOK)
			{
				// int count = kMerIt->second;

				bool present = graph->kMerinGraph(kMerSeq);

				if(present)
					optimality_score += 1;

				optimality_sum += 1;
			}
		}
	};


	auto evaluateAndReduce = [&]() {
		assert(resolved.first.length() == resolved.second.length());

		size_t size_before = openPairs.size();
		if(openPairs.size() > 0)
		{
			std::vector<double> optimality;
			optimality.resize(openPairs.size());

			long long max_i = optimality.size() - 1;
			long long chunk_size = max_i / threads;

			#pragma omp parallel
			{
				assert(omp_get_num_threads() == threads);
				long long thisThread = omp_get_thread_num();
				long long firstPair = thisThread * chunk_size;
				long long lastPair = (thisThread+1) * chunk_size - 1;
				if((thisThread == (threads-1)) && (lastPair < max_i))
				{
					lastPair = max_i;
				}

				for(long long i = firstPair; i <= lastPair; i++)
				{
					int copyCharactersFromPrevious = k;
					if((copyCharactersFromPrevious > resolved_noGaps.first.length()) || (copyCharactersFromPrevious > resolved_noGaps.second.length()))
					{
						copyCharactersFromPrevious = (resolved_noGaps.first.length() < resolved_noGaps.second.length()) ? resolved_noGaps.first.length() : resolved_noGaps.second.length() ;
					}

					std::string previous_first = resolved_noGaps.first.substr(resolved_noGaps.first.length() - copyCharactersFromPrevious, copyCharactersFromPrevious);
					std::string previous_second = resolved_noGaps.second.substr(resolved_noGaps.second.length() - copyCharactersFromPrevious, copyCharactersFromPrevious);

					assert(previous_first.length() == copyCharactersFromPrevious);
					assert(previous_second.length() == copyCharactersFromPrevious);
					assert(previous_first == previous_second);

					std::string seq_evaluate_1 = previous_first + openPairs.at(i).first;
					std::string seq_evaluate_2 = previous_second + openPairs.at(i).second;

					int optimality_sum_seq1 = 0;
					int optimality_score_seq1 = 0;
					int optimality_sum_seq2 = 0;
					int optimality_score_seq2 = 0;
					evaluateSequence2(seq_evaluate_1, optimality_sum_seq1, optimality_score_seq1);
					evaluateSequence2(seq_evaluate_2, optimality_sum_seq2, optimality_score_seq2);

					double O = 1;
					if((optimality_sum_seq1 + optimality_sum_seq2) > 0)
					{
						O = (double)(optimality_score_seq1 + optimality_score_seq2)/(double)(optimality_sum_seq1 + optimality_sum_seq2);
					}

					optimality.at(i) = O;
				}
			}

			double maxOptimality; size_t whereOptimality;
			for(size_t i = 0; i < optimality.size(); i++)
			{
				if((i == 0) || (maxOptimality < optimality.at(i)))
				{
					maxOptimality = optimality.at(i);
					whereOptimality = i;
				}
			}

			if(maxOptimality == 0)
			{
				std::vector<std::pair<std::string, std::string> > new_openPairs;
				new_openPairs.push_back(openPairs.at(0));
				openPairs = new_openPairs;
			}
			else
			{
				std::vector<std::pair<std::string, std::string> > new_openPairs;
				if(optimality.size() <= 1000)
				{
					for(size_t i = 0; i < optimality.size(); i++)
					{
						if((optimality.at(i)/maxOptimality) > 0.9)
						{
							new_openPairs.push_back(openPairs.at(i));
						}
					}
				}
				else
				{
					std::vector<unsigned int> optimality_sort_indices;
					optimality_sort_indices.resize(optimality.size());
					for(unsigned int i = 0; i < optimality_sort_indices.size(); i++)
					{
						optimality_sort_indices.at(i) = i;
					}
					std::sort(optimality_sort_indices.begin(), optimality_sort_indices.end(), [&](unsigned int a, unsigned int b) -> bool {
							return (optimality.at(b) < optimality.at(a));
						}
					);

					assert(optimality_sort_indices.size() == optimality.size());
					unsigned int maxIndex = 1000;
					for(unsigned int optimalityI = 0; optimalityI <= maxIndex; optimalityI++)
					{
						unsigned int realOptimalityI = optimality_sort_indices.at(optimalityI);
						double realOptimality = optimality.at(realOptimalityI);
						if(optimalityI > 0)
						{
							unsigned int previous_OptimalityI = optimality_sort_indices.at(optimalityI - 1);
							double previous_Optimality = optimality.at(previous_OptimalityI);
							assert(previous_Optimality >= realOptimality);
						}

						new_openPairs.push_back(openPairs.at(realOptimalityI));
					}
				}
				openPairs = new_openPairs;
			}

			size_t size_after = openPairs.size();
			if(((double)size_after/(double)size_before) > 0.5)
			{
				std::cerr << "\nevaluateAndReduce(..): Not very effective. Reduction from " << size_before << " to " << size_after << ", maxOptimality = " << maxOptimality << " -- but does this really matter??\n" << std::flush;
//				std::cerr << "\t" << "resolved.first.length(): " << resolved.first.length() << "\n";
//				std::cerr << "\t" << "openPairs.at(0).first.length(): " << openPairs.at(0).first.length() << " // " << openPairs.at(0).first << "\n" << std::flush;
//				assert(size_after > 100);
//				std::vector<std::pair<std::string, std::string> > new_openPairs;
//				for(size_t i = 0; i < openPairs.size(); i++)
//				{
//					if(Utilities::randomDouble() < 0.1)
//					{
//						new_openPairs.push_back(openPairs.at(i));
//					}
//				}
//				openPairs = new_openPairs;
//				std::cerr << "Random reduction to " << openPairs.size() << " elements. " << "\n" << std::flush;
			}
		}
	};

	auto evaluateAndResolve = [&]() {
		assert(resolved.first.length() == resolved.second.length());

		if(openPairs.size() > 0)
		{
			std::vector<double> optimality;
			optimality.resize(openPairs.size());

			long long max_i = optimality.size() - 1;
			long long chunk_size = max_i / threads;

			#pragma omp parallel
			{
				assert(omp_get_num_threads() == threads);
				long long thisThread = omp_get_thread_num();
				long long firstPair = thisThread * chunk_size;
				long long lastPair = (thisThread+1) * chunk_size - 1;
				if((thisThread == (threads-1)) && (lastPair < max_i))
				{
					lastPair = max_i;
				}

				for(long long i = firstPair; i <= lastPair; i++)
				{
					int copyCharactersFromPrevious = k;
					if((copyCharactersFromPrevious > resolved_noGaps.first.length()) || (copyCharactersFromPrevious > resolved_noGaps.second.length()))
					{
						copyCharactersFromPrevious = (resolved_noGaps.first.length() < resolved_noGaps.second.length()) ? resolved_noGaps.first.length() : resolved_noGaps.second.length() ;
					}

					std::string previous_first = resolved_noGaps.first.substr(resolved_noGaps.first.length() - copyCharactersFromPrevious, copyCharactersFromPrevious);
					std::string previous_second = resolved_noGaps.second.substr(resolved_noGaps.second.length() - copyCharactersFromPrevious, copyCharactersFromPrevious);

					assert(previous_first.length() == copyCharactersFromPrevious);
					assert(previous_second.length() == copyCharactersFromPrevious);
					assert(previous_first == previous_second);

					std::string seq_evaluate_1 = previous_first + openPairs.at(i).first;
					std::string seq_evaluate_2 = previous_second + openPairs.at(i).second;

					double O = evaluateSequence(seq_evaluate_1) + evaluateSequence(seq_evaluate_2);

					optimality.at(i) = O;
				}
			}

			double maxOptimality; size_t whereOptimality;
			for(size_t i = 0; i < optimality.size(); i++)
			{
				if((i == 0) || (maxOptimality < optimality.at(i)))
				{
					maxOptimality = optimality.at(i);
					whereOptimality = i;
				}
			}

			resolved.first.append(openPairs.at(whereOptimality).first);
			resolved.second.append(openPairs.at(whereOptimality).second);


			std::string chosen_s1 = openPairs.at(whereOptimality).first;
			std::string chosen_s2 = openPairs.at(whereOptimality).second;

			std::string chosen_s1_noGaps = chosen_s1;
			std::string chosen_s2_noGaps = chosen_s2;
			chosen_s1_noGaps.erase(std::remove_if(chosen_s1_noGaps.begin(),chosen_s1_noGaps.end(), [&](char c){return ((c == '_') ? true : false);}), chosen_s1_noGaps.end());
			chosen_s2_noGaps.erase(std::remove_if(chosen_s2_noGaps.begin(),chosen_s2_noGaps.end(), [&](char c){return ((c == '_') ? true : false);}), chosen_s2_noGaps.end());
			resolved_noGaps.first.append(chosen_s1_noGaps);
			resolved_noGaps.second.append(chosen_s2_noGaps);


			std::string chosen_hom;

			size_t last_character_s1 = chosen_s1.size();
			size_t last_character_s2 = chosen_s2.size();
			size_t have_hom_s1 = 0;
			size_t have_hom_s2 = 0;

			if(addedHomozygousSeq > 0)
			{
				do{
					last_character_s1--;
					if(chosen_s1.at(last_character_s1) != '_')
					{
						have_hom_s1++;
					}
				} while (have_hom_s1 != addedHomozygousSeq);

				do{
					last_character_s2--;
					if(chosen_s2.at(last_character_s2) != '_')
					{
						have_hom_s2++;
					}
				} while (have_hom_s2 != addedHomozygousSeq);

				std::string suffix_s1 = chosen_s1.substr(last_character_s1);
				std::string suffix_s2 = chosen_s2.substr(last_character_s2);

				assert(suffix_s1 == suffix_s2);
				std::string original_chosen_1 = chosen_s1;
				chosen_hom = suffix_s1;

				chosen_s1 = chosen_s1.substr(0, last_character_s1);
				chosen_s2 = chosen_s2.substr(0, last_character_s2);
				assert((chosen_s1 + chosen_hom) == original_chosen_1);
				assert(chosen_hom.length() >= addedHomozygousSeq);
			}

			std::vector<std::string> newCompartment;
			newCompartment.push_back(chosen_s1);
			newCompartment.push_back(chosen_s2);
			resolved_gS.push_back(newCompartment);

			if(chosen_hom.size() > 0)
			{
				std::vector<std::string> newHomCompartment;
				newHomCompartment.push_back(chosen_hom);
				resolved_gS.push_back(newHomCompartment);
			}

			openPairs.clear();
		}
		addedHomozygousSeq = 0;
	};

	auto addOpenPair = [&](std::pair<std::string, std::string> p) {
		assert(p.first.length() == p.second.length());

		if(openPairs.size() == 0)
		{
			openPairs.push_back(p);
		}
		else
		{
			std::vector<std::pair<std::string, std::string> > new_openPairs;

			for(unsigned int existingI = 0; existingI < openPairs.size(); existingI++)
			{
				std::pair<std::string, std::string> newpair_1 = openPairs.at(existingI);
				std::pair<std::string, std::string> newpair_2 = openPairs.at(existingI);

				newpair_1.first.append(p.first);
				newpair_1.second.append(p.second);

				newpair_2.first.append(p.second);
				newpair_2.second.append(p.first);

				new_openPairs.push_back(newpair_1);
				new_openPairs.push_back(newpair_2);
			}

			openPairs = new_openPairs;

			if(openPairs.size() > 100000)
			{
				evaluateAndReduce();
			}

			assert(openPairs.size() < 10000000);
		}

		addedHomozygousSeq = 0;
	};

	auto addOpenSeq = [&](std::string seq) {
		for(unsigned int existingI = 0; existingI < openPairs.size(); existingI++)
		{
			openPairs.at(existingI).first.append(seq);
			openPairs.at(existingI).second.append(seq);
		}

		for(int i = 0; i < seq.length(); i++)
		{
			if(seq.at(i) != '_')
				addedHomozygousSeq++;
		}

		if(addedHomozygousSeq >= k)
		{
			evaluateAndResolve();
		}
	};

	size_t expectedResolvedLength = 0;
	std::cout << "\n";
	for(unsigned int compartmentI = 0; compartmentI < original_gS.size(); compartmentI++)
	{
		_compartmentI = compartmentI;

		if((compartmentI % 100) == 0)
		{
			std::cout << "\r" << "compartmentI " << compartmentI << " / " << original_gS.size() << "   openPairs.size(): " << openPairs.size() << "     " << std::flush;
		}

		assert(original_gS.at(compartmentI).size() >= 1);
		assert(original_gS.at(compartmentI).size() <= 2);

		if(original_gS.at(compartmentI).size() == 2)
		{
			assert(original_gS.at(compartmentI).at(0).size() == original_gS.at(compartmentI).at(1).size());
		}

		unsigned int compartmentLength = original_gS.at(compartmentI).at(0).size();
		bool diploid = (original_gS.at(compartmentI).size() == 2);

		if(diploid)
		{
			std::pair<std::string, std::string> oP;
			oP.first = original_gS.at(compartmentI).at(0);
			oP.second = original_gS.at(compartmentI).at(1);

			addOpenPair(oP);
		}
		else
		{
			std::string homozygousSeq = original_gS.at(compartmentI).at(0);
			if(openPairs.size() > 0)
			{
				addOpenSeq(homozygousSeq);
			}
			else
			{
				if(resolved_gS.size() == 0)
				{
					std::vector<std::string> firstCompartment;
					firstCompartment.push_back(homozygousSeq);
					resolved_gS.push_back(firstCompartment);
				}
				else
				{
					if(resolved_gS.at(resolved_gS.size() - 1).size() == 1)
					{
						resolved_gS.at(resolved_gS.size() - 1).at(0).append(homozygousSeq);
					}
					else
					{
						std::vector<std::string> newCompartment;
						newCompartment.push_back(homozygousSeq);
						resolved_gS.push_back(newCompartment);
					}
				}

				resolved.first.append(homozygousSeq);
				resolved.second.append(homozygousSeq);

				std::string homozygousSeq_noGaps = homozygousSeq;
				homozygousSeq_noGaps.erase(std::remove_if(homozygousSeq_noGaps.begin(),homozygousSeq_noGaps.end(), [&](char c){return ((c == '_') ? true : false);}), homozygousSeq_noGaps.end());

				resolved_noGaps.first.append(homozygousSeq_noGaps);
				resolved_noGaps.second.append(homozygousSeq_noGaps);
			}
		}

		expectedResolvedLength += compartmentLength;

		// std::cout << "greedilyResolveDiploidKMerString(..): Compartment " << compartmentI << " of " << original_gS.size() << "\n";
		// std::cout << "\topenPairs.size(): " << openPairs.size() << "\n";
		// std::cout << "\taddedHomozygousSeq: " << addedHomozygousSeq << "\n";
		// std::cout << "\texpectedResolvedLength: " << expectedResolvedLength << "\n";
		// std::cout << std::flush;
	}

	std::cout << "\n";

	evaluateAndResolve();
	assert(resolved.first.length() == expectedResolvedLength);
	assert(resolved.second.length() == expectedResolvedLength);

	std::pair<diploidGenomeString, diploidGenomeString> forReturn;
	std::vector<string> forReturn_oneCompartment;
	forReturn_oneCompartment.push_back(resolved.first);
	forReturn_oneCompartment.push_back(resolved.second);
	forReturn.first.push_back(forReturn_oneCompartment);
	forReturn.second = resolved_gS;

	std::cout << "\tgreedilyResolveDiploidKMerString(..): length of gS before: " << gS_length_before << ", and after: " << resolved_gS.size() << "\n" << std::flush;


	return forReturn;
}




#endif /* VALIDATOR_H_ */
