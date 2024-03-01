# Genome assembly algorithm

Assembly algorithm destined to work on single-end reads originating from the same strand of a single chromosome written in C++.

Made for Bioinformatics algorithms for genomic data analysis classes at MIMUW in the winter semester of 2023/2024.

## Assumptions:
Typical parameters of input dataset:
- number of reads: 1000,
- read length: 80bp,
- average percentage of mismatches: $\leq$ 5%,
- average coverage: $\geq$ 5x,
- alignments of length $<$ 300bp shoild be excluded.

## Short description:
### Correcting Reads
The first step involves correcting the reads. I base it on the `correct1mm` function presented in the lecture. The parameters chosen for the `correct_reads` function are:
- `k = 19`: the length of kmers whose frequency in reads we compare.
- `threshold = 2`: the frequency value of a kmer, exceeding which it will be considered correct.

```
std::vector<std::string> correct_reads(std::vector<std::string>&& reads) {
    auto timer = Timer("Correcting reads");
    std::vector<std::string> corrected;
    auto histogram = count_kmers(reads);
    for (auto& read : reads) {
        for (size_t i = 0; i <= read.size() - k; ++i) {
            std::string kmer = read.substr(i, k);
            if (map_get(histogram, kmer, (uint32_t) 0) > threshold)
                continue;
            auto neighbours = get_neighbours(kmer);
            for (const auto& neighbour : neighbours) {
                if (map_get(histogram, neighbour, (uint32_t) 0) > threshold) {
                    read.replace(i, k, neighbour);
                    break;
                }
            }

        }
        corrected.push_back(read);         
    }
    return corrected;
}
```

### Building Overlap Graph

The next step is to build an overlap graph from the corrected reads. The graph is directed. An edge comes out from the vertex of a read whose suffix overlaps with the prefix of the connecting read. To create edges, I use the KMP algorithm to find occurrences of a pattern in a text. The parameter indicating the minimum overlap between reads taken into account is **11**.

An optimization is to remove redundant edges, those directly connecting two vertices that could be traversed by passing through another vertex. A mention of such reduction appears in the article [1].

```
    void removeRedundantEdges() {
        for (auto& [_, node] : nodes) {
            auto candidates = node.getSortedSuccessors();
            for (auto it1 = candidates.begin(); it1 != candidates.end(); ++it1) {
                for (auto it2 = candidates.rbegin(); it2 != candidates.rend(); ++it2) {
                    if (*it1 == *it2) continue;
                    if (node.outsContain(*it2) && nodes.at(*it2).outsContain(*it1)) {
                        node.removeOut(*it1);
                        nodes.at(*it1).removeIn(node.id);
                    }
                }
            }
        }
    }
```

[1] Rizzi, R., Beretta, S., Patterson, M., Pirola, Y., Previtali, M., Della Vedova, G. and Bonizzoni, P. (2019), Overlap graphs and de Bruijn graphs: data structures for de novo genome assembly in the big data era. Quantitative Biology, 7: 278-292.

### Creating Contigs

For creating contigs, I use a simple greedy algorithm. To establish contig assembly, I start with the read that has the smallest ID among those with the fewest incoming edges.

```
    while (!graph.empty()) {
        auto id = graph.getSmallestSinkId();
        auto node = graph.removeNode(id);
        contigs.push_back(node.sequence);
        while (node.hasOutgoingEdge()) {
            auto [id, w] = node.getBiggesOut();
            node = graph.removeNode(id);
            contigs.back() += node.sequence.substr(w);
        }
    }
```

Since only contigs equal to or longer than 300 bp are considered for evaluation, I save only such contigs to the file.

## Usage:
```
g++ -O3 -Wall -Wextra -Wpedantic -std=c++17 assembly.cpp -o assembly
```
```
./assembly reads_path contigs_path
```