#include "deBruijn/DeBruijnGraph.h"
#include "deBruijn/DeBruijnGraphAlt.h"

#include "lib/bioio.hpp"
#include "lib/cxxopts.h"

#include "spdlog/spdlog.h"
#include "spdlog/stopwatch.h"

#include <cassert>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sys/resource.h>
#include <thread>
#include <vector>

int main( int argc, char **argv ) {

    // TODO add program option
    // set log level to trace
    spdlog::set_level( spdlog::level::trace );

    // TODO add program option
    unsigned int thread_count = 4;

    // TODO add program option
    std::string textFilePath = "ecoli1.fna";

    // TODO add program option
    int kmerL = 30;

    //------------------------------------------------------------------------------------------------------------------

    // load text from file
    std::string text;
    {
        spdlog::stopwatch sw;
        std::ifstream myfile( textFilePath );
        auto record = bioio::read_fasta( myfile );
        for ( const auto &it : record ) {
            text.append( it.sequence );
        }

        spdlog::info( "Record({}) loaded in {} seconds", record.size(), sw );
        spdlog::info( "Text loaded:" );
        spdlog::info( "Bytes    : {}", text.size() );
        spdlog::info( "KiloBytes: {:.4f}", text.size() / 1024.0 );
        spdlog::info( "MegaBytes: {:.4f}", text.size() / std::pow(1024, 2) );
        spdlog::info( "GigaBytes: {:.4f}", text.size() / std::pow(1024, 3) );
    }

    std::string fail = "AGGCCCTGAAGC";
    std::string fail2 = "TAAGCTGATGTT"; // 4 good, 3bad
    std::string fail3 = "ATGCTGTAGCTAGATATCGTAGCTATGCTAGCTAATAGCTATTTCGATGCGGTAGCTAGTGCTAGCATGCGTATGCATGCGTACGGCTAGCTAG"
                        "TAGAGCTCGACTACGACGACGAGAGGGCATCGACGATTAGAGACTAGCGACTACGAGCTAGCGACT";

    //------------------------------------------------------------------------------------------------------------------

    {
        spdlog::stopwatch sw;

        auto graph = DeBruijnGraphAlt::create( std::move( text ), kmerL, thread_count );
        // auto graph = DeBruijnGraphAlt::create( std::move( text ), kmerL, 1 );

        spdlog::info( "Text buidling took: {} seconds", sw );
    }

    //------------------------------------------------------------------------------------------------------------------

    // auto a = DeBruijnGraphAlt( fail3, 4 );
    // std::cout << "Graph build" << a.kmerToNode.size() << std::endl;
    // auto tour = a.hasEulerianWalkdOrCycle();
    // LOG(INFO) << "HEAD:     " + a.head->kmer ;
    // a.toDot();
    // TODO add tour to_dot
    // TODO find out if g is multimap or just 1 to many
}
