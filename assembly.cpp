#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>
#include <unordered_map>
#include <algorithm>


constexpr uint32_t overlap = 11;
constexpr uint32_t k = 19;
constexpr uint32_t threshold = 2;
constexpr uint32_t min_superstring_length = 300;
constexpr char alphabet[] = {'A', 'C', 'G', 'T'};

class Timer {
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    std::string msg;
public:
    Timer(std::string&& _msg) : start(std::chrono::high_resolution_clock::now()), msg(_msg) {}
    ~Timer() {
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        std::cout << msg << " took " << duration << " ms" << std::endl;
    }
};

std::pair<std::string, std::string> read_args(int argc, const char** argv) {
    if (argc != 3) {
        std::cout << "Usage: " << argv[0] << " <reads file> <output file>" << std::endl;
        exit(1);
    }
    return std::make_pair(argv[1], argv[2]);
}

std::vector<std::string> load_reads(const std::string& path) {
    auto timer = Timer("Loading reads");
    std::ifstream reads_file(path);
    std::vector<std::string> reads;

    std::string line;
    while (std::getline(reads_file, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') continue;
        reads.push_back(line);
    }

    return reads;
}

void write_contigs(const std::string& path, const std::vector<std::string>& contigs) {
    auto timer = Timer("Writing contigs");
    std::ofstream contigs_file(path);
    for (size_t i = 0; i < contigs.size(); ++i) {
        contigs_file << ">contig_" << i << std::endl;
        contigs_file << contigs[i] << std::endl;
    }
}

class GraphOfOverlaps;

class Node {
public:
    uint32_t id;
    std::string sequence;
private:
    std::unordered_map<uint32_t, uint32_t> out;
    std::unordered_map<uint32_t, uint32_t> ins;
public:
    Node(uint32_t _id, std::string&& _sequence) : id(_id), sequence(_sequence) {}

    void addEdge(Node& other, uint32_t weight) {
        out[other.id] = weight;
        other.ins[id] = weight;
    }

    bool hasOutgoingEdge() const {
        return !out.empty();
    }

    std::vector<uint32_t> getSortedSuccessors() const {
        std::vector<std::pair<uint32_t, uint32_t>> tmp(out.begin(), out.end());
        std::sort(tmp.begin(), tmp.end(), [](const auto& lhs, const auto& rhs) {
            return lhs.second < rhs.second;
        });
        std::vector<uint32_t> successors;
        for (const auto& [id, _] : tmp) {
            successors.push_back(id);
        }
        return successors;
    }

    void removeIn(uint32_t id) {
        ins.erase(id);
    }

    void removeOut(uint32_t id) {
        out.erase(id);
    }

    bool outsContain(uint32_t id) const {
        return out.find(id) != out.end();
    }

    std::pair<uint32_t, uint32_t> getBiggesOut() const {
        std::pair<uint32_t, uint32_t> biggest(0, 0);
        for (const auto& [id, weight] : out) {
            if (weight > biggest.second) {
                biggest = std::make_pair(id, weight);
            }
        }
        return biggest;
    }

    size_t getInsCount() const {
        return ins.size();
    }

    friend bool operator==(const Node& lhs, const Node& rhs) {
        return lhs.id == rhs.id || lhs.sequence == rhs.sequence;
    }

    friend class GraphOfOverlaps;
};


class GraphOfOverlaps {
public:
    std::unordered_map<uint32_t, Node> nodes;
public:
    GraphOfOverlaps(std::vector<std::string>&& reads) {
        auto timer = Timer("Building graph of overlaps");
        for (uint32_t i = 0; i < reads.size(); ++i) {
            nodes.insert({i, Node(i, std::move(reads[i]))});
        }
    }

    std::unordered_map<uint32_t, Node>& getNodes() {
        return nodes;
    }

    bool empty() const {
        return nodes.empty();
    }

    void removeRedundantEdges() {
        for (auto& [_, node] : nodes) {
            auto candidates = node.getSortedSuccessors();
            for (auto it1 = candidates.begin(); it1 != candidates.end(); ++it1) {
                for (auto it2 = candidates.begin(); it2 != candidates.end(); ++it2) {
                    if (*it1 == *it2) continue;
                    if (node.outsContain(*it2) && nodes.at(*it2).outsContain(*it1)) {
                        node.removeOut(*it1);
                        nodes.at(*it1).removeIn(node.id);
                    }
                }
            }
        }
    }
        
    uint32_t getSmallestSinkId() {
        uint32_t smallest = UINT32_MAX;
        size_t lowest_ins = SIZE_MAX;
        for (auto& [_, node] : nodes) {
            if (node.getInsCount() <= lowest_ins) {
                if (node.getInsCount() == lowest_ins){
                    if (node.id < smallest) {
                        smallest = node.id;
                }
                } else {
                    lowest_ins = node.getInsCount();
                    smallest = node.id;
                }
            } 
        }
        return smallest;
    }

    Node removeNode(uint32_t id) {
        auto node = nodes.at(id);
        for (const auto& [id, _] : node.ins) {
            nodes.at(id).removeOut(node.id);
        }
        for (const auto& [id, _] : node.out) {
            nodes.at(id).removeIn(node.id);
        }
        nodes.erase(id);
        return node;
    }
};

std::vector<size_t> lps_table(const std::string& pattern) {
    std::vector<size_t> lps(pattern.size());
    size_t len = 0;
    lps[0] = 0;
    size_t i = 1;
    while (i < pattern.size()) {
        if (pattern[i] == pattern[len]) {
            ++len;
            lps[i] = len;
            ++i;
        } else {
            if (len != 0) {
                len = lps[len - 1];
            } else {
                lps[i] = 0;
                ++i;
            }
        }
    }
    return lps;
}

std::vector<size_t> find_occurences(
    const std::string& text, 
    const std::string& pattern,
    const std::vector<size_t>& lps
) {
    std::vector<size_t> occurences;
    size_t i = 0;
    size_t j = 0;
    while (i < text.size()) {
        if (pattern[j] == text[i]) {
            ++i;
            ++j;
        }
        if (j == pattern.size()) {
            occurences.push_back(i - j);
            j = lps[j - 1];
        } else if (i < text.size() && pattern[j] != text[i]) {
            if (j != 0) {
                j = lps[j - 1];
            } else {
                ++i;
            }
        }
    }
    return occurences;
}

bool ends_with(const std::string& str, const std::string& suffix) {
    if (str.size() < suffix.size()) return false;
    return 0 == str.compare(str.size() - suffix.size(), suffix.size(), suffix);
}

std::vector<std::string> assemble(std::vector<std::string>&& reads) {
    auto timer = Timer("Assembling");
    std::vector<std::string> contigs;
    GraphOfOverlaps graph(std::move(reads));
    
    for (auto& [_, node] : graph.getNodes()) {
        const std::string suffix = node.sequence.substr(node.sequence.size() - overlap);
        const auto lps = lps_table(suffix);
        for (auto& [_, other] : graph.getNodes()) {
            if (node == other) continue;
            const auto occurences = find_occurences(other.sequence, suffix, lps);
            for (const auto& occurence : occurences) {
                if (ends_with(node.sequence, other.sequence.substr(0, occurence + overlap))) {
                    node.addEdge(other, occurence + overlap);
                    break;
                }
            }
        }
    }
    
    graph.removeRedundantEdges();

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

    std::vector<std::string> filtered_contigs;
    for (const auto& contig : contigs) {
        if (contig.size() >= min_superstring_length) {
            filtered_contigs.push_back(contig);
        }
    }
    return filtered_contigs;
}

std::unordered_map<std::string, uint32_t> count_kmers(std::vector<std::string>& reads) {
    std::unordered_map<std::string, uint32_t> kmers;
    for (const auto& read : reads) {
        for (size_t i = 0; i <= read.size() - k; ++i) {
            std::string kmer = read.substr(i, k);
            if (kmers.find(kmer) == kmers.end()) {
                kmers[kmer] = 1;
            } else {
                ++kmers[kmer];
            }
        }
    }
    return kmers;
}

std::vector<std::string> get_neighbours(const std::string& kmer) {
    std::vector<std::string> neighbours;
    for (size_t i = 0; i < kmer.size() ; ++i) {
        for (const auto& c : alphabet) {
            if (kmer[i] == c) continue;
            std::string neighbour = kmer;
            neighbour[i] = c;
            neighbours.push_back(neighbour);
        }
    }
    return neighbours;
}

uint32_t map_get(const std::unordered_map<std::string, uint32_t>& map, const std::string& key, const uint32_t& default_value) {
    auto it = map.find(key);
    if (it == map.end()) {
        return default_value;
    }
    return it->second;
}

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

int main(int argc, const char** argv) {
    auto timer = Timer("Total");
    auto args = read_args(argc, argv);
    auto reads = load_reads(args.first);
    auto corrected = correct_reads(std::move(reads));
    auto contigs = assemble(std::move(corrected));
    write_contigs(args.second, contigs);
}