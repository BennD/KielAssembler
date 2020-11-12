#include "DeBruijnGraphAlt.h"

#include <cmath>
#include <fstream>
#include <functional>
#include <future>
#include <sstream>
#include <thread>

#include "spdlog/fmt/ostr.h" // overload << for std::thread::id
#include "spdlog/spdlog.h"

DeBruijnGraphAlt::DeBruijnGraphAlt( const std::string_view sequenceToAssemble, size_t kmerLength )
    : m_kmer_length( kmerLength ) {

    auto thread_id = std::this_thread::get_id();
    spdlog::trace( "Thread {} Graph-build: started", thread_id );

    for ( size_t i = 0; i < sequenceToAssemble.size() - ( m_kmer_length + 1 ); i++ ) {
        auto kmerL = sequenceToAssemble.substr( i, m_kmer_length );
        auto kmerR = sequenceToAssemble.substr( i + 1, m_kmer_length );
        auto iNodeL = find_or_create_node( kmerL );
        auto iNodeR = find_or_create_node( kmerR );
        m_edgesOut[iNodeL].push_back( iNodeR );
        m_edgesIn[iNodeR].push_back( iNodeL );
    }

    spdlog::trace( "Thread {} Graph-build: finished", thread_id );
}

void DeBruijnGraphAlt::calculate_graph_properties() {
    // check if graph has an eulerian walk/cycle
    // find head/tail
    {
        // size_t balanced = 0;
        size_t semiBalanced = 0;
        size_t neither = 0;

        for ( size_t i = 0; i < m_kmer.size(); i++ ) {
            if ( is_node_balanced( i ) ) {
                // balanced++;
            } else if ( is_node_semi_balanced( i ) ) {
                semiBalanced++;

                if ( get_node_out_degree( i ) > get_node_in_degree( i ) ) {
                    m_head = i;
                } else {
                    m_tail = i;
                }
            } else {
                neither++;
            }
        }
        m_hasEulerianWalk = ( neither == 0 && semiBalanced == 2 );
        m_hasEulerianCycle = ( neither == 0 && semiBalanced == 0 );
    }
}

void DeBruijnGraphAlt::set_sequence( std::string &&sequence ) {
    m_sequence = sequence;
}

size_t DeBruijnGraphAlt::find_or_create_node( std::string_view kmer, std::optional<Hash> hashOpt ) {
    size_t index;
    auto hash = hashOpt.value_or( std::hash<std::string_view>{}( kmer ) );
    const auto it = m_kmerMap.find( hash );

    if ( it == m_kmerMap.end() ) {
        index = create_node( kmer, hash );
    } else {
        index = it->second;
    }

    return index;
}

size_t DeBruijnGraphAlt::create_node( std::string_view kmer, Hash hash ) {
    auto index = m_kmer.size();

    m_kmerMap[hash] = index;
    m_kmer.emplace_back( kmer );
    m_kmerHashes.emplace_back( hash );
    m_edgesIn.emplace_back();
    m_edgesOut.emplace_back();
    m_mergedWith.emplace_back( NotMerged );
    m_isActive.emplace_back( true );

    return index;
}

bool DeBruijnGraphAlt::is_node_balanced( size_t index ) const {
    long inDegree = get_node_in_degree( index );
    long outDegree = get_node_out_degree( index );

    return inDegree == outDegree;
}

bool DeBruijnGraphAlt::is_node_semi_balanced( size_t index ) const {
    long inDegree = get_node_in_degree( index );
    long outDegree = get_node_out_degree( index );

    return std::abs( inDegree - outDegree ) == 1;
}

size_t DeBruijnGraphAlt::get_node_degree( size_t index ) const {
    return get_node_in_degree( index ) + get_node_out_degree( index );
}

size_t DeBruijnGraphAlt::get_node_in_degree( size_t index ) const {
    return m_edgesIn[index].size();
}

size_t DeBruijnGraphAlt::get_node_out_degree( size_t index ) const {
    return m_edgesOut[index].size();
}

DeBruijnGraphAlt DeBruijnGraphAlt::create( std::string &&sequence, size_t kmerLength, size_t thread_count ) {
    if ( thread_count == 0 ) {
        thread_count = 1;
    }

    // no multithreading if one thread is selected
    if ( thread_count == 1 ) {
        auto graph = DeBruijnGraphAlt( std::string_view( sequence ), kmerLength );
        graph.set_sequence( std::move( sequence ) );
        return graph;
    }

    // create sub views
    std::vector<std::string_view> subViews;
    {
        auto sequenceView = std::string_view( sequence );
        auto viewLength = sequenceView.size() / thread_count;

        // create first view
        subViews.push_back( sequenceView.substr( 0, viewLength ) );

        for ( size_t i = 1; i < ( thread_count - 1 ); i++ ) {
            // subsequent views overlap the previous views by one character
            subViews.push_back( sequenceView.substr( ( i * viewLength ) - kmerLength, viewLength ) );
        }

        // last view is longer if (sequence.size() % thread_count != 0)
        // .substring() caps the view to the length of the sequence
        auto lastIndex = thread_count - 1;
        subViews.push_back(
            sequenceView.substr( ( lastIndex * viewLength ) - kmerLength, std::numeric_limits<size_t>::max() ) );
    }

    assert( subViews.size() == thread_count );

    // generate sub graphs
    std::vector<DeBruijnGraphAlt> subGraphs;
    {
        std::vector<std::future<DeBruijnGraphAlt>> futures;
        auto lambda = [=]( std::string_view subSequence ) { return DeBruijnGraphAlt( subSequence, kmerLength ); };

        for ( auto view : subViews ) {
            futures.emplace_back( std::async( std::launch::async, lambda, view ) );
        }

        for ( auto &future : futures ) {
            subGraphs.emplace_back( future.get() );
        }
    }

    auto graph = std::move( subGraphs.front() );

    for( size_t i = 1; i < subGraphs.size(); i++ ) {
        spdlog::info("What is this?");
        graph.merge_graph( subGraphs[i] );
    }

    return graph;
}

bool DeBruijnGraphAlt::is_eulerian() const {
    return has_eulerian_cycle() || has_eulerian_walk();
}

bool DeBruijnGraphAlt::has_eulerian_walk() const {
    return m_hasEulerianWalk;
}

bool DeBruijnGraphAlt::has_eulerian_cycle() const {
    return m_hasEulerianCycle;
}

std::vector<size_t> DeBruijnGraphAlt::get_euler_path() const {
    // stack St;
    std::vector<size_t> stack;
    std::vector<size_t> euler_path;
    auto edges = m_edgesOut;

    // put start vertex in St;
    stack.push_back( 0 );

    // until St is empty
    while ( !stack.empty() ) {

        // let V be the value at the top of St;
        const auto v = stack.back();

        auto &v_out_edges = edges[v];

        // if degree(V) = 0, then % probably meant to be outDegree
        if ( v_out_edges.empty() ) {
            // add V to the answer;
            euler_path.push_back( v );
            // remove V from the top of St;
            stack.pop_back();
        }
        // otherwise
        else {
            // find any edge coming out of V;
            auto tgt = v_out_edges.back();
            // remove it from the graph;
            v_out_edges.pop_back();
            // put the second end of this edge in St;
            stack.push_back( tgt );
        }
    }

    return euler_path;
}

void DeBruijnGraphAlt::toDot( const std::string &filename ) const {
    std::ofstream file;
    file.open( filename, std::ios::out );
    file << "digraph {\n";

    for ( size_t i = 0; i < m_edgesOut.size(); i++ ) {
        for ( const auto &j : m_edgesOut[i] ) {
            file << m_kmer[i] << "->" << m_kmer[j] << "\n";
        }
    }

    file << "}";
    file.close();
}

void DeBruijnGraphAlt::merge_graph( const DeBruijnGraphAlt &otherGraph ) {
    for ( const auto &altKmerNode : otherGraph.m_kmerMap ) {
        auto const &hash = altKmerNode.first;
        auto const &altKmerId = altKmerNode.second;
        auto const& kmer = otherGraph.m_kmer[altKmerId];

        // add node from other graph to this graph
        auto const kmerId = find_or_create_node( kmer, hash );

        // add all outgoing edges
        for ( const auto &altTargetKmerId : otherGraph.m_edgesOut[altKmerId] ) {
            auto const &targetKmer = otherGraph.m_kmer[altTargetKmerId];
            auto const &targetHash = otherGraph.m_kmerHashes[altTargetKmerId];
            auto const targetKmerId = find_or_create_node( targetKmer, targetHash );
            m_edgesOut[kmerId].push_back( targetKmerId );
            m_edgesIn[targetKmerId].push_back( kmerId );
        }
    }
}
