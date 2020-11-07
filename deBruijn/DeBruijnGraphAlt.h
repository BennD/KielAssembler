//
// Created by Benno Doerr on 10/21/2020.
//

#ifndef KIELASSEMBLER_GRAPH_H
#define KIELASSEMBLER_GRAPH_H

#include <limits>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

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

  public:
    DeBruijnGraphAlt() = delete;
    DeBruijnGraphAlt( const DeBruijnGraphAlt & ) = delete;
    DeBruijnGraphAlt( DeBruijnGraphAlt &&graph ) = default;

    DeBruijnGraphAlt& operator=( const DeBruijnGraphAlt & ) = delete;
    DeBruijnGraphAlt& operator=( DeBruijnGraphAlt && ) = default;

    /**
     * Generate graph.
     *
     * Graph is generated without creating additional threads if thread_count <= 1
     *
     * @sequence view of sequence used to generate graph on
     * @kmerLength initial length of kmers
     * @thread_count number of threads used for graph generation
     */
    static DeBruijnGraphAlt create( std::string &&sequence, size_t kmerLength, size_t thread_count = 1 );

    /**
     * True if graph is eularian.
     */
    bool is_eulerian() const;

    /**
     * True if graph has eulerian walk.
     */
    bool has_eulerian_walk() const;

    /**
     * True if graph has eulerian cycle.
     */
    bool has_eulerian_cycle() const;

    /**
     * Generate euler path.
     */
    std::vector<size_t> get_euler_path() const;

    /**
     * Save dot file of graph to file.
     */
    void toDot( const std::string &filename ) const;

  private:
    DeBruijnGraphAlt( const std::string_view sequenceToAssemble, size_t kmerLength );

    void calculate_graph_properties();
    void set_sequence( std::string &&sequence );
    size_t find_or_create_node( std::string_view kmer );
    size_t create_node( std::string_view kmer );
    bool is_node_balanced( size_t index ) const;
    bool is_node_semi_balanced( size_t index ) const;
    size_t get_node_degree( size_t index ) const;
    size_t get_node_in_degree( size_t index ) const;
    size_t get_node_out_degree( size_t index ) const;
};

#endif // KIELASSEMBLER_GRAPH_H
