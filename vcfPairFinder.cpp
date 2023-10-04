#include <htslib/vcf.h>
#include <map>
#include <list>
#include <string>
#include <iostream>
#include <set>

void printUsage(const char *path)
{
    std::cerr << "\nProgram: vcfPairFinder v1.0.0\n"
              << std::endl;
    std::cerr << "Usage: " << path << " --vcf <VCF File> [--max-bp-dist <Maximum Base Pair Distance>] \n"
              << std::endl;

    std::cerr << "\nDescription:" << std::endl;
    std::cerr << "  Identifies pairs of genetic heterozygous variants within a specified distance from each other in a VCF file.\n"
              << std::endl;

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


     // Map with sample name as key and a set of variant IDs as value
    std::map<std::string, std::set<std::pair<int, std::string>>> heterozygous_variants;

    while (bcf_read(inf, hdr, rec) == 0) {
        bcf_unpack(rec, BCF_UN_STR);
        int variant_position = rec->pos + 1;
        std::string ref = rec->d.allele[0];
        std::string alt = rec->d.allele[1];
        std::string variant_id = "chr" + std::to_string(rec->rid + 1) + ":" +
                                std::to_string(variant_position) + ":" + ref + ":" + alt;

        int ngt = 0;
        int *genotype = NULL;
        ngt = bcf_get_genotypes(hdr, rec, &genotype, &ngt);

        for (int i = 0; i < bcf_hdr_nsamples(hdr); i++) {
            int *sample_genotype = genotype + i * 2; // Assuming diploid organisms
            int allele1 = bcf_gt_allele(sample_genotype[0]);
            int allele2 = bcf_gt_allele(sample_genotype[1]);

            bool is_heterozygous = allele1 != allele2 && allele1 >= 0 && allele2 >= 0;
            std::string sample_name = std::string(hdr->samples[i]);

            if (is_heterozygous) {
                heterozygous_variants[sample_name].emplace(variant_position, variant_id);
            }
        }

        free(genotype);
    }

    // print header
    std::cout << "sample" << "\t" << "dist" << "\t"
        << "rsid1" << "\t" << "rsid2" << std::endl;

    // Now iterate over samples, then over their heterozygous variants
    for (const auto &sample_entry : heterozygous_variants) {
        const auto &sample_name = sample_entry.first;
        const auto &variants = sample_entry.second;

        for (auto it = variants.begin(); it != variants.end(); ++it) {
            const auto &current_variant = *it;

            // Iterate over the remaining variants in the set to find close variants
            for (auto jt = std::next(it); jt != variants.end(); ++jt) {
                const auto &other_variant = *jt;

                // If distance is within the limit, process the variant pair
                if (abs(other_variant.first - current_variant.first) <= maxBpDist) {
                    int distance = other_variant.first - current_variant.first;
                    std::cout << sample_name << "\t" << distance << "\t"
                            << current_variant.second << "\t" << other_variant.second << std::endl;
                    // Since the set is sorted, and we have output this pair, 
                    // there is no need to check them again in the reverse order.
                    break; 
                }
                // Since the set is sorted, we can break early if we passed the possible close variants
                else if (other_variant.first - current_variant.first > maxBpDist) {
                    break;
                }
            }
        }
    }


    // Release resources
    bcf_hdr_destroy(hdr);
    bcf_destroy(rec);
    hts_close(inf);

    return 0;
}

