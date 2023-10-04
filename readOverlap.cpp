
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <sstream>
#include <iostream>
#include <zlib.h>
#include <vector>
#include <iterator>
#include <chrono>

// Function to split a string by a delimiter character
std::vector<std::string> split(const std::string &s, char delimiter)
{
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter))
    {
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

int main(int argc, char *argv[])
{

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
    if (!file)
    {
        std::cerr << "Could not open input file" << std::endl;
        return 1;
    }

    // Open the CRAM file
    bam_hdr_t *header = NULL;
    htsFile *in = hts_open(cramFile.c_str(), "r");

    // set reference
    if (hts_set_fai_filename(in, refFile.c_str()) == -1)
    {
        std::cerr << "Could not set reference genome file. "
                  << "Ensure the file exists, is readable, and indexed (with .fai file)." << std::endl;
        return 1;
    }

    header = sam_hdr_read(in);

    hts_idx_t *idx = sam_index_load(in, cramFile.c_str());
    if (idx == NULL)
    {
        std::cerr << "Couldn't load index" << std::endl;
        return 1;
    }
    char buf[1024];
    while (gzgets(file, buf, sizeof(buf)))
    {
        std::stringstream ss(buf);
        std::string sample, rsid1, rsid2;
        int dist;
        ss >> sample >> dist >> rsid1 >> rsid2;

        // Parse rsid1 and rsid2 to construct region
        std::vector<std::string> tokens1 = split(rsid1, ':');
        std::vector<std::string> tokens2 = split(rsid2, ':');
        if (tokens1.size() < 2 || tokens2.size() < 2)
        {
            std::cerr << "Invalid rsid format" << std::endl;
            continue;
        }

        // setup region
        std::string region = tokens1[0] + ":" + tokens1[1] + "-" + tokens2[1];
        hts_itr_t *iter = sam_itr_querys(idx, header, region.c_str());
        bam1_t *b = bam_init1();

        int count_cis = 0;
        int count_trans = 0;

        int pos1 = stoi(tokens1[1]) - 1; // Convert to 0-based index
        int pos2 = stoi(tokens2[1]) - 1; // Convert to 0-based index

        int total_reads_covering_pos1 = 0;
        int total_reads_covering_pos2 = 0;

        while (sam_itr_next(in, iter, b) >= 0)
        {
            uint8_t *seq = bam_get_seq(b);
            uint32_t q2 = b->core.qual; // mapping quality

            // Adjust indices to be relative to start of read
            int read_start = b->core.pos; // start position of read
            int adjusted_pos1 = pos1 - read_start;
            int adjusted_pos2 = pos2 - read_start;

            // Check if the read covers pos1 
            if (adjusted_pos1 >= 0 && adjusted_pos1 < b->core.l_qseq)
            {
                total_reads_covering_pos1++;
            }

            // check if the read covers pos2
            if (adjusted_pos2 >= 0 && adjusted_pos2 < b->core.l_qseq)
            {
                total_reads_covering_pos2++;
            }

            // We are interested in the reads covering two variants at once
            if (adjusted_pos1 >= 0 && adjusted_pos1 < b->core.l_qseq && adjusted_pos2 >= 0 && adjusted_pos2 < b->core.l_qseq)
            {
                int index1 = bam_seqi(seq, adjusted_pos1);
                int index2 = bam_seqi(seq, adjusted_pos2);

                // get nucletide for pos1 and pos2
                char nt1 = seq_nt16_str[index1];
                char nt2 = seq_nt16_str[index2];

                // std::cerr << nt1 << "\t" << region << std::endl;
                // std::cerr << nt2 << "\t" << region << std::endl;

                // Comparing nucleotides; you might need more sophisticated comparison for real data
                if (nt1 == tokens1[2][0] && nt2 == tokens2[2][0])
                {
                    count_cis++;
                }
                else
                {
                    count_trans++;
                }
            }
        }

        // After the loop, you might want to output or store the counts
        std::cout << rsid1 << "\t" << rsid2 << "\t" << count_cis << "\t" << count_trans
            << "\t" << total_reads_covering_pos1 << "\t" << total_reads_covering_pos2 << std::endl;

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
