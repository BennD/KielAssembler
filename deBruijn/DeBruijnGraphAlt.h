//
// Created by Benno Doerr on 10/21/2020.
//

#ifndef KIELASSEMBLER_GRAPH_H
#define KIELASSEMBLER_GRAPH_H

#include <limits>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

class DeBruijnGraphAlt {
  public:
    using Hash = size_t;
    using KmerId = size_t;

    static constexpr size_t NotMerged = std::numeric_limits<size_t>::max();

    /**
     * Basis for all kmer string_views.
     */
    std::string m_sequence;

    /**
     * Maps hashes to kmer ids.
     */
    std::unordered_map<Hash, KmerId> m_kmerMap;

    /**
     * Maps kmer ids to kmers.
     */
    std::vector<std::string_view> m_kmer;

    /**
     * Maps kmer ids to kmer hashes.
     */
    std::vector<Hash> m_kmerHashes;

    /**
     * List all ingoing egdges.
     * edges[i][j] -> edges[i]
     */
    std::vector<std::vector<KmerId>> m_edgesIn;

    /**
     * List all outgoing egdges.
     * edges[i] -> edges[i][j]
     */
    std::vector<std::vector<KmerId>> m_edgesOut;

    /**
     * Kmer nodes can be chained with other kmer nodes.
     * DeBruijngraphalt::NotMerged if node has not been chained to another node.
     */
    std::vector<KmerId> m_mergedWith;

    /**
     * Kmer nodes are dead once they have been absorbed by another kmer node.
     */
    std::vector<bool> m_isActive;

    /**
     * Head of the graph. Needed for eulerian cycle.
     */
    size_t m_head;

    /**
     * Tail of the graph. Needed for X.
     * TODO
     */
    size_t m_tail;

    /**
     * Kmer length.
     * TODO does this have any meaning?
     */
    size_t m_kmer_length;

    /**
     * True if graph has a eulerian walk.
     */
    bool m_hasEulerianWalk;

    /**
     * True if graph has a eulerian cycle.
     */
    bool m_hasEulerianCycle;

  public:
    DeBruijnGraphAlt() = delete;
    DeBruijnGraphAlt( const DeBruijnGraphAlt & ) = delete;
    DeBruijnGraphAlt( DeBruijnGraphAlt &&graph ) = default;

    DeBruijnGraphAlt &operator=( const DeBruijnGraphAlt & ) = delete;
    DeBruijnGraphAlt &operator=( DeBruijnGraphAlt && ) = default;

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

    /**
     * Needs to be run once after graph has been assembled/merged. Sets
     * - head
     * - tail
     * - has eulerian cycle
     * - has eulerian walk
     */
    void calculate_graph_properties();
    void set_sequence( std::string &&sequence );
    size_t find_or_create_node( std::string_view kmer, std::optional<Hash> hashOpt = std::nullopt );
    size_t create_node( std::string_view kmer, Hash hash );
    bool is_node_balanced( size_t index ) const;
    bool is_node_semi_balanced( size_t index ) const;
    size_t get_node_degree( size_t index ) const;
    size_t get_node_in_degree( size_t index ) const;
    size_t get_node_out_degree( size_t index ) const;

    void merge_graph( const DeBruijnGraphAlt &otherGraph );
};

#endif // KIELASSEMBLER_GRAPH_H
