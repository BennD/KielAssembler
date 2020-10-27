//
// Created by Benno Doerr on 10/21/2020.
//

#ifndef KIELASSEMBLER_GRAPH_H
#define KIELASSEMBLER_GRAPH_H

#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include <functional>
#include <future>
#include <sstream>
#include <thread>

class DeBruijnGraphAlt {
  public:
    using Hash = size_t;
    using Index = size_t;

    static constexpr size_t NONE = std::numeric_limits<size_t>::max();

    std::string m_sequence;
    // Maps hashes to ids
    std::unordered_map<Hash, Index> m_kmerMap;
    std::vector<std::string_view> m_kmer;
    std::vector<std::vector<Index>> m_edgesIn;
    std::vector<std::vector<Index>> m_edgesOut;
    std::vector<Index> m_mergedWith;
    std::vector<bool> m_isActive;

    size_t m_head;
    size_t m_tail;

    size_t m_kmer_length;

    bool m_hasEulerianWalk;
    bool m_hasEulerianCycle;

  private:
    DeBruijnGraphAlt( const std::string_view sequenceToAssemble, size_t kmerPairLength )
        : m_kmer_length( kmerPairLength - 1 ) {

        for ( size_t i = 0; i < sequenceToAssemble.size() - m_kmer_length; i++ ) {
            auto kmerL = sequenceToAssemble.substr( i, m_kmer_length );
            auto kmerR = sequenceToAssemble.substr( i + 1, m_kmer_length );
            auto iNodeL = find_or_create_node( kmerL );
            auto iNodeR = find_or_create_node( kmerR );
            m_edgesOut[iNodeL].push_back( iNodeR );
            m_edgesIn[iNodeR].push_back( iNodeL );
        }

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

    void set_sequence( std::string &&sequence ) {
        m_sequence = sequence;
    }

    size_t find_or_create_node( std::string_view kmer ) {
        size_t index;
        const auto it = m_kmerMap.find( std::hash<std::string_view>{}( kmer ) );

        if ( it == m_kmerMap.end() ) {
            index = create_node( kmer );
        } else {
            index = it->second;
        }

        return index;
    }

    size_t create_node( std::string_view kmer ) {
        auto index = m_kmer.size();

        m_kmerMap[std::hash<std::string_view>{}( kmer )] = index;

        m_kmer.emplace_back( kmer );
        m_edgesIn.emplace_back();
        m_edgesOut.emplace_back();
        m_mergedWith.emplace_back( NONE );
        m_isActive.emplace_back( true );

        return index;
    }

    bool is_node_balanced( size_t index ) const {
        return m_edgesIn[index].size() == m_edgesOut[index].size();
    }

    bool is_node_semi_balanced( size_t index ) const {
        return std::abs( static_cast<long>( m_edgesIn[index].size() ) -
                         static_cast<long>( m_edgesOut[index].size() ) ) == 1;
    }

    size_t get_node_degree( size_t index ) {
        return get_node_in_degree( index ) + get_node_out_degree( index );
    }

    size_t get_node_in_degree( size_t index ) {
        return m_edgesIn[index].size();
    }

    size_t get_node_out_degree( size_t index ) {
        return m_edgesOut[index].size();
    }

  public:
    DeBruijnGraphAlt( DeBruijnGraphAlt &&graph ) {
        m_sequence = graph.m_sequence;
        m_kmer = graph.m_kmer;
        m_kmerMap = graph.m_kmerMap;
        m_edgesIn = graph.m_edgesIn;
        m_edgesOut = graph.m_edgesOut;
        m_isActive = graph.m_isActive;
        m_mergedWith = graph.m_mergedWith;
        m_head = graph.m_head;
        m_tail = graph.m_tail;
        m_kmer_length = graph.m_kmer_length;
    }

    static DeBruijnGraphAlt create( std::string &&sequence, size_t kmerPairLength, size_t thread_count = 0 ) {
        auto graph = DeBruijnGraphAlt( std::string_view( sequence ), kmerPairLength );
        graph.set_sequence( std::move( sequence ) );
        return graph;
    }

    bool is_eulerian() {
        return has_eulerian_cycle() || has_eulerian_walk();
    }

    bool has_eulerian_walk() {
        return m_hasEulerianWalk;
    }

    bool has_eulerian_cycle() {
        return m_hasEulerianCycle;
    }

    /*
    std::future<std::unique_ptr<DeBruijnGraphAlt>> merge( std::unique_ptr<DeBruijnGraphAlt> graph_to_merge ) {

        // Highest index in "merged onto" graph
        auto index = m_kmer.size();
        int counter = 0;
        for ( const auto &currentNode : graph_to_merge->m_kmerMap ) {

            auto kmerToAdd = m_kmerMap.find( currentNode.first );
            if ( kmerToAdd != m_kmerMap.end() ) {
                // Add edges
                auto inEdgesToAdd = graph_to_merge->m_edgesIn[currentNode.second];
                auto outEdgesToAdd = graph_to_merge->m_edgesOut[currentNode.second];
                m_edgesOut[kmerToAdd->second].insert( outEdgesToAdd.begin(), outEdgesToAdd.end(),
                                                      outEdgesToAdd.begin() );
                m_edgesIn[kmerToAdd->second].insert( inEdgesToAdd.begin(), inEdgesToAdd.end(), inEdgesToAdd.begin() );
                // copy edges onto us
            } else {
                // If we dont find entry with the same hash, we add it
                // and add the max Index tour the right index, so that its stays unique
                // TODO But I think it will fuck up the vector
                m_kmerMap[kmerToAdd->first] = index + counter;
                // Add edges
                auto inEdgesToAdd = graph_to_merge->m_edgesIn[currentNode.second];
                auto outEdgesToAdd = graph_to_merge->m_edgesOut[currentNode.second];
                m_edgesOut[kmerToAdd->second].insert( outEdgesToAdd.begin(), outEdgesToAdd.end(),
                                                      outEdgesToAdd.begin() );
                m_edgesIn[kmerToAdd->second].insert( inEdgesToAdd.begin(), inEdgesToAdd.end(), inEdgesToAdd.begin() );
                // add the string view
                m_kmer.push_back( graph_to_merge->m_kmer[kmerToAdd->second] );

                ++counter;
            }
        }
    }
    */

    /*
    std::vector<size_t> get_euler_path() const {
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
    */

    void toDot( const std::string &filename ) {
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
};

#endif // KIELASSEMBLER_GRAPH_H
