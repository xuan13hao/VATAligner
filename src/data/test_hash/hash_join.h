#include <iostream>
#include <vector>
#include <unordered_map>

// Tuple structure representing a record in a relation
struct Tuple {
    int key;
    std::vector<int> value;

    Tuple(int k, const std::vector<int>& v) : key(k), value(v) {}
};

// Hash join algorithm
std::vector<std::pair<Tuple, Tuple>> hashJoin(const std::vector<Tuple>& relation1, const std::vector<Tuple>& relation2) {
    // Build hash table for relation1
    std::unordered_map<int, std::vector<Tuple>> hashTable;
    for (const auto& tuple : relation1) {
        hashTable[tuple.key].push_back(tuple);
    }

    // Perform the hash join
    std::vector<std::pair<Tuple, Tuple>> result;
    for (const auto& tuple : relation2) {
        int key = tuple.key;
        auto it = hashTable.find(key);

        // Check if there are matching tuples in relation1
        if (it != hashTable.end()) {
            for (const auto& matchingTuple : it->second) {
                result.push_back({matchingTuple, tuple});
            }
        }
    }

    return result;
}

int main() {
    // Example relations
    std::vector<Tuple> relation1 = { {1, {10, 20, 30}}, {2, {40, 50}}, {3, {60, 70, 80}} };
    std::vector<Tuple> relation2 = { {2, {90, 100}}, {3, {110, 120, 130}}, {4, {140, 150}} };

    // Perform hash join
    std::vector<std::pair<Tuple, Tuple>> result = hashJoin(relation1, relation2);

    // Display the result
    for (const auto& pair : result) {
        std::cout << "(" << pair.first.key << ", [";
        for (int val : pair.first.value) {
            std::cout << val << ", ";
        }
        std::cout << "]) JOIN (";
        std::cout << pair.second.key << ", [";
        for (int val : pair.second.value) {
            std::cout << val << ", ";
        }
        std::cout << "])" << std::endl;
    }

    return 0;
}
