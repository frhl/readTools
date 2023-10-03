#include <htslib/sam.h>
#include <htslib/hts.h>
#include <sstream>
#include <iostream>
#include <zlib.h>
#include <vector>
#include <iterator>
#include <chrono>

// Function to split a string by a delimiter character
std::vector<std::string> split(const std::string &s, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

void printUsage(const char *path)
{
    std::cerr << "\nProgram: readOverlap\n"
              << std::endl;
    std::cerr << "Usage: " << path << " --pairs <Pairs File> --cram <CRAM File> \n"
              << std::endl;

    std::cerr << "\nDescription:" << std::endl;
    std::cerr << "  Counts reads supporting pairs of variants.\n"
              << std::endl;

    std::cerr << "\nInput:" << std::endl;
    std::cerr << "  --pairs/-p <Pairs File>     : Required. Input gzipped file of variant pairs." << std::endl;
    std::cerr << "  --cram/-c <CRAM File>       : Required. Input CRAM file." << std::endl;
    std::cerr << "  --ref/-r <Reference File>   : Required. Input reference genome file." << std::endl;
}

int main(int argc, char *argv[]) {

    if (argc < 2)
    {
        printUsage(argv[0]);
        return 1;
    }

    std::string pairsFile;
    std::string cramFile;
    std::string refFile;

    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if ((arg == "--pairs" || arg == "-p") && i + 1 < argc) 
        {
            pairsFile = argv[++i];
        }
        else if ((arg == "--cram" || arg == "-c") && i + 1 < argc)
        {
            cramFile = argv[++i];
        }
        else if ((arg == "--ref" || arg == "-r") && i + 1 < argc)
        {
            refFile = argv[++i];
        }
        else
        {
            printUsage(argv[0]);
            return 1;
        }
    }

    if (pairsFile.empty() || cramFile.empty() || refFile.empty())
    {
        std::cerr << "Error: pairs fil, CRAM and reference file paths are required." << std::endl;
        printUsage(argv[0]);
        return 1;
    }

    // Open the gzipped input file
    gzFile file = gzopen(pairsFile.c_str(), "r");
    if (!file) {
        std::cerr << "Could not open input file" << std::endl;
        return 1;
    }

    // Open the CRAM file
    std::cerr << "Opening .cram file.." << std::endl;
    bam_hdr_t *header = NULL;
    htsFile *in = hts_open(cramFile.c_str(), "r");
    
    // set reference 
    if (hts_set_fai_filename(in, refFile.c_str()) == -1) {
        std::cerr << "Could not set reference genome file. "
                << "Ensure the file exists, is readable, and indexed (with .fai file)." << std::endl;
        return 1;
    }

    header = sam_hdr_read(in);

    std::cerr << "Loading .crai file.." << std::endl;
    hts_idx_t *idx = sam_index_load(in, cramFile.c_str());
    if (idx == NULL) {
        std::cerr << "Couldn't load index" << std::endl;
        return 1;
    }

    std::cerr << "Reading through gzfile" << std::endl;
    char buf[1024];
    while (gzgets(file, buf, sizeof(buf))) {
        std::stringstream ss(buf);
        std::string sample, rsid1, rsid2;
        int dist;
        ss >> sample >> dist >> rsid1 >> rsid2;

        // Parse rsid1 and rsid2 to construct region
        std::vector<std::string> tokens1 = split(rsid1, ':');
        std::vector<std::string> tokens2 = split(rsid2, ':');
        if (tokens1.size() < 2 || tokens2.size() < 2) {
            std::cerr << "Invalid rsid format" << std::endl;
            continue; // skip this line or handle error
        }

        std::string chr = tokens1[0];
        std::string start_pos = tokens1[1];
        std::string end_pos = tokens2[1];
        std::string region = chr + ":" + start_pos + "-" + end_pos;

	    hts_itr_t *iter = sam_itr_querys(idx, header, region.c_str());

        bam1_t *b = bam_init1();

        int count_both = 0;
        int count_rsid1 = 0;
        int count_rsid2 = 0;

        std::cerr << "Opening region " << region << std::endl;
        while (sam_itr_next(in, iter, b) >= 0) {
            uint8_t *seq = bam_get_seq(b);
            uint32_t q2 = b->core.qual ; //mapping quality

            std::cerr << "Mapping quality: " << q2 << std::endl;

            int pos1 = stoi(tokens1[1]) - 1; // Convert to 0-based index
            int pos2 = stoi(tokens2[1]) - 1; // Convert to 0-based index

            std::cerr << "converted index...." << std::endl;

            if (seq == NULL) {
                std::cerr << "seq is NULL. Skipping iteration." << std::endl;
                continue;
            }

            std::cerr << "Checking variables before accessing sequence..." << std::endl;
            if (b == NULL || seq == NULL) {
                std::cerr << "Null pointer detected. Skipping iteration." << std::endl;
                continue; // Skip to the next iteration if null pointer detected.
            }

            std::cerr << "Validating indices..." << std::endl;
            // You might need a function to get the sequence length, adjust accordingly.
            if (pos1 < 0 || pos2 < 0 /* || pos1 >= seq_length || pos2 >= seq_length */) {
                std::cerr << "Invalid index detected. Skipping iteration." << std::endl;
                continue; // Skip to the next iteration if index is invalid.
            }

            int index1 = bam_seqi(seq, pos1);
            int index2 = bam_seqi(seq, pos2);

            std::cerr << "Attempting to access sequence data..." << std::endl;

            // You may also need to consider the read orientation and other factors
            char nt1 = seq_nt16_str[bam_seqi(seq, pos1)];
            char nt2 = seq_nt16_str[bam_seqi(seq, pos2)];

            std::cerr << "Sequence data accessed successfully." << std::endl;

            // Comparing nucleotides; you might need more sophisticated comparison for real data
            if (nt1 == tokens1[2][0] && nt2 == tokens2[2][0]) {
                count_both++;
            } else if (nt1 == tokens1[2][0]) {
                count_rsid1++;
            } else if (nt2 == tokens2[2][0]) {
                count_rsid2++;
            }
        }

        // After the loop, you might want to output or store the counts
        std::cout << rsid1 << "\t" << rsid2 << "\t" << count_both << "\t" << count_rsid1 << "\t" << count_rsid2 << std::endl;


        // Clean up for this iteration
        bam_destroy1(b);
        hts_itr_destroy(iter);
    }

    // Global clean up
    gzclose(file);
    hts_idx_destroy(idx);
    sam_hdr_destroy(header);
    hts_close(in);

    return 0;
}

