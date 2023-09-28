#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <sstream>
#include <zlib.h>
#include <cmath>
#include <set>
#include <tuple>

struct Variant
{
    std::string chromosome;
    long position;
    std::string reference;
    std::string alternative;
};

Variant parseVariant(const std::string &variant_str)
{
    Variant variant;
    std::stringstream ss(variant_str);
    std::string token;
    std::getline(ss, variant.chromosome, ':');
    std::getline(ss, token, ':');
    variant.position = std::stol(token);
    std::getline(ss, variant.reference, ':');
    std::getline(ss, variant.alternative, ':');
    return variant;
}

void printUsage(const char *path)
{
    std::cerr << "\nProgram: variantPairFinder v1.0.0\n"
              << std::endl;
    std::cerr << "Usage: " << path << " --geno <Genotype File> [--max-bp-dist <Maximum Base Pair Distance>] \n"
              << std::endl;

    std::cerr << "\nDescription:" << std::endl;
    std::cerr << "  This program identifies pairs of genetic variants within a specified distance from each other in a gzipped file.\n"
              << std::endl;

    std::cerr << "\nInput:" << std::endl;
    std::cerr << "  --geno/-g <Genotype File>    : Required. A gzipped file containing genotype data." << std::endl;
    std::cerr << "                                The file should have sample id, variant id, and genotype separated by whitespace." << std::endl;
    std::cerr << "                                Each line might look like the following: 'Sample1 chr21:12314:C:T 1|0'." << std::endl;
    std::cerr << "  --max-bp-dist/-d <Distance>  : Optional. Maximum base pair distance between variant pairs. Default is 500." << std::endl;

    std::cerr << "\nOutput Format:" << std::endl;
    std::cerr << "  The output will be printed to the console with the following columns:" << std::endl;
    std::cerr << "  Sample, Distance, Position1, Position2, Variant1, Variant2" << std::endl;
}

int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        printUsage(argv[0]);
        return 1;
    }

    std::string genoFile;
    long maxBpDist = 500;
    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if (arg == "--geno" && i + 1 < argc)
        {
            genoFile = argv[++i];
        }
        else if (arg == "--max-bp-dist" && i + 1 < argc)
        {
            maxBpDist = std::stol(argv[++i]);
        }
        else
        {
            printUsage(argv[0]);
            return 1;
        }
    }

    if (genoFile.empty())
    {
        printUsage(argv[0]);
        return 1;
    }

    gzFile file = gzopen(genoFile.c_str(), "rb");
    if (!file)
    {
        std::cerr << "Error: Cannot open gzipped file for reading: " << genoFile << std::endl;
        return 1;
    }

    std::map<std::string, std::vector<Variant>> sample_variants;
    char buf[1024];
    while (gzgets(file, buf, sizeof(buf)))
    {
        std::stringstream ss(buf);
        std::string sample, variant_str, genotype;
        ss >> sample >> variant_str >> genotype;
        if (genotype == "1|0" || genotype == "0|1")
        {
            Variant variant = parseVariant(variant_str);
            sample_variants[sample].push_back(variant);
        }
    }
    gzclose(file);

    std::set<std::tuple<std::string, std::string, std::string>> printedPairs;
    int variantSet = 0; // Initialize variantSet counter

    // Iterate through each sample and validate ordering, then find variant pairs
    for (const auto &sample_pair : sample_variants)
    {
        const std::vector<Variant> &variants = sample_pair.second;

        // Validate that variants are ordered by position
        for (size_t i = 1; i < variants.size(); ++i)
        {
            if (variants[i - 1].position > variants[i].position)
            {
                std::cerr << "Error: Variants for sample " << sample_pair.first << " are not ordered by position. Please sort the input.\n";
                return 1;
            }
        }

        // Find variant pairs within the specified distance
        for (size_t i = 0; i < variants.size(); ++i)
        {
            bool newSetStarted = false;
            // Iterate only in one direction since we track the processed pairs
            for (size_t k = i + 1; k < variants.size() && variants[k].position - variants[i].position <= maxBpDist; ++k)
            {
                if (variants[i].chromosome == variants[k].chromosome)
                {
                    long distance = std::abs(variants[i].position - variants[k].position);
                    // Construct a pair of variant identifiers
                    std::string variantId1 = variants[i].chromosome + ':' + std::to_string(variants[i].position) + ':' + variants[i].reference + ':' + variants[i].alternative;
                    std::string variantId2 = variants[k].chromosome + ':' + std::to_string(variants[k].position) + ':' + variants[k].reference + ':' + variants[k].alternative;
                    std::string sampleId = sample_pair.first; // or however you get the sample ID
                    if (printedPairs.find({sampleId, variantId1, variantId2}) == printedPairs.end() &&
                        printedPairs.find({sampleId, variantId2, variantId1}) == printedPairs.end())
                    {

                        if (!newSetStarted)
                        {
                            variantSet++; // Start a new variant set
                            newSetStarted = true;
                        }

                        // Print the pair and add it to the set of printed pairs
                        std::cout << sampleId << '\t'
                                  << distance << '\t'
                                  << variantId1 << '\t'
                                  << variantId2 << '\t'
                                  << variantSet << std::endl;
                        printedPairs.insert({sampleId, variantId1, variantId2});
                    }
                }
            }
        }
    }

    return 0;
}


