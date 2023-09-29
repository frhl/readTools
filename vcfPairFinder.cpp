#include <htslib/vcf.h>
#include <map>
#include <list>
#include <string>
#include <iostream>


void printUsage(const char *path)
{
    std::cerr << "\nProgram: vcfPairFinder v1.0.0\n" << std::endl;
    std::cerr << "Usage: " << path << " --vcf <VCF File> [--max-bp-dist <Maximum Base Pair Distance>] \n" << std::endl;

    std::cerr << "\nDescription:" << std::endl;
    std::cerr << "  Identifies pairs of genetic heterozygous variants within a specified distance from each other in a VCF file.\n" << std::endl;

    std::cerr << "\nInput:" << std::endl;
    std::cerr << "  --vcf/-v <VCF File>          : Required. Input VCF file." << std::endl;
    std::cerr << "  --max-bp-dist/-d <Distance>  : Optional. Maximum base pair distance between variant pairs. Default is 500." << std::endl;

    std::cerr << "\nOutput Format:" << std::endl;
    std::cerr << "  Sample, Distance, Variant1, Variant2" << std::endl;
}


int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        printUsage(argv[0]);
        return 1;
    }

    std::string vcfFile;
    long maxBpDist = 500;
    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if ((arg == "--vcf" || arg == "-v") && i + 1 < argc) // Ensure parentheses here
        {
            vcfFile = argv[++i];
        }
        else if ((arg == "--max-bp-dist" || arg == "-d") && i + 1 < argc) // And here
        {
            maxBpDist = std::stol(argv[++i]);
        }
        else
        {
            printUsage(argv[0]);
            return 1;
        }
    }

    if (vcfFile.empty())
    {
        std::cerr << "Error: VCF file path is required." << std::endl;
        printUsage(argv[0]);
        return 1;
    }

    htsFile *inf = hts_open(vcfFile.c_str(), "r");
    bcf_hdr_t *hdr = bcf_hdr_read(inf);
    bcf1_t *rec = bcf_init();


    // Map with sample name as key and a list of variant IDs as value
    std::map<std::string, std::list<std::string>> window_variants;

    while (bcf_read(inf, hdr, rec) == 0)
    {
        bcf_unpack(rec, BCF_UN_STR);
        std::string ref = rec->d.allele[0];
        std::string alt = rec->d.allele[1];
        int variant_position = rec->pos + 1; // 1-based position
	    std::string variant_id = "chr" + std::to_string(rec->rid + 1) + ":" +
                         std::to_string(variant_position) + ":" + ref + ":" + alt;


        int ngt = 0;
        int *genotype = NULL;
        ngt = bcf_get_genotypes(hdr, rec, &genotype, &ngt);
        int nalleles = ngt / bcf_hdr_nsamples(hdr);

        for (int i = 0; i < bcf_hdr_nsamples(hdr); i++)
        {
            int *sample_genotype = genotype + i * nalleles;
            bool is_heterozygous = nalleles == 2 && sample_genotype[0] != sample_genotype[1];
            std::string sample_name = std::string(hdr->samples[i]);

            // Skip if not heterozygous
            if (!is_heterozygous)
                continue;

            auto it = window_variants[sample_name].begin();
            while (it != window_variants[sample_name].end())
            {
                int other_variant_position = std::stoi(it->substr(it->find(':') + 1));
                if (abs(variant_position - other_variant_position) > maxBpDist)
                {
                    it = window_variants[sample_name].erase(it);
                }
                else
                {
                    for (const auto &pos : window_variants[sample_name])
                    {
                        int pos_variant_position = std::stoi(pos.substr(pos.find(':') + 1));
                        if (pos_variant_position < other_variant_position)
                        {
                            int distance = other_variant_position - pos_variant_position;
                            std::cout << sample_name << "\t" << distance << "\t" << pos << "\t" << *it << std::endl;
                        }
                    }
                    ++it;
                }
            }
            window_variants[sample_name].push_back(variant_id);
        }

        free(genotype);
    }

    // Release resources
    bcf_hdr_destroy(hdr);
    bcf_destroy(rec);
    hts_close(inf);

    return 0;
}

