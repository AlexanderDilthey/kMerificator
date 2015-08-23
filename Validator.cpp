#include "Validator.h"

#include <vector>
#include <map>
#include <set>
#include "Utilities.h"

#include "hash/deBruijn/DeBruijnGraph.h"

#include <boost/math/distributions/poisson.hpp>

using namespace boost::math::policies;
using namespace boost::math;

typedef boost::math::poisson_distribution< double, policy < discrete_quantile < integer_round_inwards> > > poisson_up;

std::vector<double> kMer_PP(size_t observedCoverage, size_t upperBound, double coverage)
{
	std::vector<double> forReturn;
	forReturn.resize(upperBound + 1, 0);

	assert(upperBound < 1000);
	assert(coverage > 0);

	double sum_to_upper = 0;
	if(observedCoverage >= (upperBound * coverage))
	{
		forReturn.clear();
		forReturn.resize(upperBound+1, 0);
		forReturn.at(upperBound) = 1;
		sum_to_upper = 1;
	}
	else
	{
		for(unsigned int underlyingCount = 0; underlyingCount <= (upperBound+1); underlyingCount++)
		{
			double lambda;
			if(underlyingCount == 0)
			{
				lambda = 0.01;
			}
			else
			{
				lambda = underlyingCount * coverage;
			}

			poisson_up poisson(lambda);
			double likelihood = pdf(poisson, observedCoverage);

			if(underlyingCount <= upperBound)
			{
				forReturn.at(underlyingCount) = likelihood;
				sum_to_upper += likelihood;
			}
			else
			{
				if(likelihood > forReturn.at(upperBound))
				{
					forReturn.clear();
					forReturn.resize(upperBound+1, 0);
					forReturn.at(upperBound) = 1;
					sum_to_upper = 1;
				}
			}
		}
	}

	if(sum_to_upper == 0)
	{
		std::cerr << "sum_to_upper = 0!!!\n";
		std::cerr << "observedCoverage: " << observedCoverage << "\n";
		std::cerr << "upperBound: " << upperBound << "\n";
		std::cerr << "coverage: " << coverage << "\n";

		for(unsigned int underlyingCount = 0; underlyingCount <= (upperBound+1); underlyingCount++)
		{
			double lambda;
			if(underlyingCount == 0)
			{
				lambda = 0.01;
			}
			else
			{
				lambda = underlyingCount * coverage;
			}

			poisson_up poisson(lambda);
			double likelihood = pdf(poisson, observedCoverage);

			std::cerr << "\t underlyingCount = " << underlyingCount << " => likelihood " << likelihood << "\n";
		}
	}

	assert(sum_to_upper != 0);


	for(unsigned int underlyingCount = 0; underlyingCount <= upperBound; underlyingCount++)
	{
		forReturn.at(underlyingCount) = forReturn.at(underlyingCount) / sum_to_upper;
	}

	double p_sum = 0;
	for(unsigned int underlyingCount = 0; underlyingCount <= upperBound; underlyingCount++)
	{
		double thisElement = forReturn.at(underlyingCount);
		assert(thisElement >= 0);
		assert(thisElement <= 1);
		p_sum += thisElement;
	}
	assert(abs(p_sum - 1) < 1e-5);

	return forReturn;
}

diploidGenomeString compressGenomeString(diploidGenomeString gS)
{
	size_t gS_length_before = gS.size();

	diploidGenomeString forReturn;
	for(size_t l = 0; l < gS.size(); l++)
	{
		if(gS.at(l).size() == 1)
		{
			if((forReturn.size() == 0) || (forReturn.at(forReturn.size() - 1).size() != 1))
			{
				forReturn.push_back(gS.at(l));
			}
			else
			{
				forReturn.at(forReturn.size() - 1).at(0).append(gS.at(l).at(0));
			}
		}
		else
		{
			forReturn.push_back(gS.at(l));
		}
	}

	std::cout << "\tcompressGenomeString(..): length of gS before: " << gS_length_before << ", and after: " << forReturn.size() << "\n" << std::flush;

	return forReturn;
}

diploidGenomeString VCF2GenomeString(std::string chromosomeID, int positionStart, int positionStop, std::string VCFpath, std::vector<std::vector<int> >& ret_graph_referencePositions, bool ignoreVCF, bool onlyPASS, const map<string, string>& referenceGenome)
{

	if(! referenceGenome.count(chromosomeID))
	{
		std::cerr << "Error: reference genome does not seem to contain chromosome " << chromosomeID << "\n";
		std::cerr << "Available chromosomes:\n";
		for(map<string, string>::const_iterator chrIt = referenceGenome.begin(); chrIt != referenceGenome.end(); chrIt++)
		{
			std::cerr << " - " << chrIt->first << "\n";
		}
		std::cerr << std::flush;
	}
	assert(referenceGenome.count(chromosomeID));

	ret_graph_referencePositions.clear();

	std::string referenceChromosome = referenceGenome.at(chromosomeID);

	std::ifstream VCFstream;

	VCFstream.open(VCFpath.c_str());

	if(! VCFstream.is_open())
	{
		throw std::runtime_error("Cannot open VCF file: "+ VCFpath);
	}

	bool ignoreBoundaries = ((positionStart == -1) && (positionStop == -1));

	if(!ignoreBoundaries)
		assert(positionStop > positionStart);

	if(ignoreBoundaries)
	{
		positionStart = 1;
		positionStop = referenceGenome.at(chromosomeID).length();
	}

	std::vector< std::vector<std::string> > chars_2_graph;

	std::string line;
	size_t lineCounter = 0;

	std::vector<std::string> header_fields;
	std::string sampleID;
	size_t fieldIndex_sample_genotypes;
	int lastExtractedPosition = 0;

	if(! ignoreVCF)
	{
		while(VCFstream.good())
		{
			std::getline(VCFstream, line);
			lineCounter++;

			Utilities::eraseNL(line);

			if(line.length() == 0)
				continue;

			if((line.length() > 1) && (line.substr(0, 2) == "##"))
				continue;

			if(line.substr(0, 1) == "#")
			{
				header_fields = Utilities::split(line, "\t");
				if(!(header_fields.size() >= 5))
				{
					std::cerr << "Not enough fields in header line of file " << VCFpath << "; expect at least 5, have " << header_fields.size() << "\n" << std::flush;
				}
				assert(header_fields.size() >= 5);

				for(unsigned int fI = 0; fI < header_fields.size(); fI++)
				{
					if(header_fields.at(fI) == "FORMAT")
					{
						assert(fI == (header_fields.size() - 2));
						fieldIndex_sample_genotypes = fI + 1;
						sampleID = header_fields.at(fieldIndex_sample_genotypes);
					}
				}
				assert(header_fields.at(3) == "REF");
				assert(header_fields.at(4) == "ALT");
				assert(header_fields.at(6) == "FILTER");

				continue;
			}

			assert(sampleID != "");

			std::vector<std::string> line_fields = Utilities::split(line, "\t");
			assert(line_fields.size() == header_fields.size());

			std::string chromosome = line_fields.at(0);
			int position = Utilities::StrtoI(line_fields.at(1));
			assert(position > 0);

			if((position % 1000) == 0)
				std::cerr << "\r" << "Chromosome " << chromosomeID << " position " << position << std::flush;

			if(chromosome != chromosomeID)
				continue;

			if(position < positionStart)
				continue;

			if(position <= lastExtractedPosition)
				continue;

			if(onlyPASS)
			{
				if(line_fields.at(6) != "PASS")
				{
					continue;
				}
			}
			if(position > positionStop)
			{
				break;
			}

			if(lastExtractedPosition != (position - 1))
			{
				if(lastExtractedPosition == 0)
				{
					for(int pos = positionStart; pos < position; pos++)
					{
						std::vector<std::string> thisPosVec;
						thisPosVec.push_back(referenceChromosome.substr(pos - 1, 1));
						chars_2_graph.push_back(thisPosVec);

						std::vector<int> thisPosVec_referencePositions;
						thisPosVec_referencePositions.push_back(pos);
						ret_graph_referencePositions.push_back(thisPosVec_referencePositions);
					}
				}
				else
				{
					for(int pos = lastExtractedPosition + 1; pos < position; pos++)
					{
						std::vector<std::string> thisPosVec;
						thisPosVec.push_back(referenceChromosome.substr(pos - 1, 1));
						chars_2_graph.push_back(thisPosVec);

						std::vector<int> thisPosVec_referencePositions;
						thisPosVec_referencePositions.push_back(pos);
						ret_graph_referencePositions.push_back(thisPosVec_referencePositions);
					}
				}
			}

			std::string reference_allele = line_fields.at(3);

			std::vector<int> thisPosVec_referencePositions;

			for(size_t pos = 0; pos < reference_allele.size(); pos++)
			{
				if(referenceChromosome.at(position + pos - 1) != reference_allele.at(pos))
				{
					std::cerr << "Reference allele error at position " << position << ", file " << VCFpath << "[in-allele position " << pos << "]\n";
					std::cerr << "\tVCF says it is " << reference_allele.at(pos) << ", but reference genome says " << referenceChromosome.at(position + pos - 1) << "\n" << std::flush;
				}
				assert(referenceChromosome.at(position + pos - 1) == reference_allele.at(pos));
			}

			std::vector<std::string> alleles;
			alleles.push_back(reference_allele);

			std::vector<std::string> alternative_alleles = Utilities::split(line_fields.at(4), ",");
			alleles.insert(alleles.end(), alternative_alleles.begin(), alternative_alleles.end());

			std::string formatString = line_fields.at(fieldIndex_sample_genotypes - 1);
			std::vector<std::string> formatString_elements = Utilities::split(formatString, ":");
			assert(formatString_elements.at(0) == "GT");

			std::string dataString = line_fields.at(fieldIndex_sample_genotypes);
			std::vector<std::string> dataString_elements = Utilities::split(dataString, ":");

			std::string genotypesString = dataString_elements.at(0);
			std::vector<std::string> genotypesString_elements = Utilities::split(genotypesString, "/");
			assert(genotypesString_elements.size() == 2);

			std::vector<std::string> called_alleles;
			unsigned int alleles_maxLength = 0;
			for(unsigned int eI = 0; eI < genotypesString_elements.size(); eI++)
			{
				int genotypeIndex = Utilities::StrtoI(genotypesString_elements.at(eI));
				std::string allele = alleles.at(genotypeIndex);
				if((eI > 0) && (allele == called_alleles.at(0)))
				{
					continue;
				}
				called_alleles.push_back(allele);
				alleles_maxLength = (allele.length() > alleles_maxLength) ? allele.length() : alleles_maxLength;
			}
			assert(alleles_maxLength != 0);

			for(unsigned int aI = 0; aI < called_alleles.size(); aI++)
			{
				int missingGaps = alleles_maxLength - called_alleles.at(aI).length();
				if(missingGaps > 0)
				{
					for(int mI = 0; mI < missingGaps; mI++)
					{
						called_alleles.at(aI).push_back('_');
					}
				}
			}

			for(unsigned int aI = 0; aI < called_alleles.size(); aI++)
			{
				assert(called_alleles.at(aI).size() == called_alleles.at(0).size());
			}

			for(unsigned int pI = 0; pI < alleles_maxLength; pI++)
			{
				if(pI < reference_allele.size())
				{
					thisPosVec_referencePositions.push_back(position + pI);
				}
				else
				{
					thisPosVec_referencePositions.push_back(-1);
				}
			}
			assert(thisPosVec_referencePositions.size() == alleles_maxLength);

			chars_2_graph.push_back(called_alleles);

			lastExtractedPosition = position + reference_allele.length() - 1;

			ret_graph_referencePositions.push_back(thisPosVec_referencePositions);
		}
	}

	if(chars_2_graph.size() > 0)
	{
		for(int pos = lastExtractedPosition + 1; pos <= positionStop; pos++)
		{
			if(pos >= (int)referenceChromosome.length())
				break;

			std::vector<std::string> thisPosVec;
			thisPosVec.push_back(referenceChromosome.substr(pos - 1, 1));
			chars_2_graph.push_back(thisPosVec);

			std::vector<int> thisPosVec_referencePositions;
			thisPosVec_referencePositions.push_back(pos);
			ret_graph_referencePositions.push_back(thisPosVec_referencePositions);
		}
	}
	else
	{
		for(int pos = positionStart; pos <= positionStop; pos++)
		{
			if(pos >= (int)referenceChromosome.length())
				break;

			std::vector<std::string> thisPosVec;
			thisPosVec.push_back(referenceChromosome.substr(pos - 1, 1));
			chars_2_graph.push_back(thisPosVec);

			std::vector<int> thisPosVec_referencePositions;
			thisPosVec_referencePositions.push_back(pos);
			ret_graph_referencePositions.push_back(thisPosVec_referencePositions);
		}
	}

	VCFstream.close();

	assert(ret_graph_referencePositions.size() == chars_2_graph.size());
	for(unsigned int segmentI = 0; segmentI < chars_2_graph.size(); segmentI++)
	{
		for(unsigned int haploI = 0; haploI < chars_2_graph.at(segmentI).size(); haploI++)
		{
			assert(chars_2_graph.at(segmentI).at(haploI).length() == ret_graph_referencePositions.at(segmentI).size());
		}
	}
	return chars_2_graph;
}


