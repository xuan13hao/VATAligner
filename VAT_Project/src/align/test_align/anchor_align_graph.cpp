#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>

using namespace std;
/*


Anchor-Align is a popular algorithm used in computational genomics to detect the order and orientation of synteny blocks, which are regions of conserved genomic structure between different species or individuals.

The Anchor-Align algorithm takes as input a set of anchor points, which are highly conserved regions in the genomes being compared, and a set of synteny blocks, which are subsequences of the genomes that share common structure. 
The algorithm constructs a graph where the anchor points are nodes and the edges represent the adjacency of the anchor points within the synteny blocks. 
The edges are weighted by the number of times that a pair of anchor points are adjacent in the synteny blocks.
The algorithm then orders the anchor points by their degree in the graph, which is the number of edges that each anchor point has. 
The intuition behind this ordering is that anchor points that are adjacent to many other anchor points are likely to be in the correct order and orientation, while those with few connections are more likely to be in the wrong order or orientation.
Finally, the algorithm outputs the ordered set of anchor points, which can be used to infer the order and orientation of the synteny blocks.

Anchor-Align has been shown to be effective in many applications of comparative genomics, such as genome assembly, gene prediction, and phylogenetics.

Create a hash table to store the indices of the anchor points.
Create an anchor-based graph and initialize its edges to 0.
Iterate over the synteny blocks and update the graph edges based on the block's anchor points.
Initialize the order of the anchor points to their original indices.
Sort the anchor points by their degree in the graph.
*/
vector<int> anchor_align(vector<int>& anchors, vector<vector<int>>& blocks) {
    // Create a hash table to store the indices of the anchor points
    unordered_map<int, int> anchor_indices;
    for (int i = 0; i < anchors.size(); i++) {
        anchor_indices[anchors[i]] = i;
    }

    // Create an anchor-based graph and initialize its edges to 0
    int num_anchors = anchors.size();
    vector<vector<int>> graph(num_anchors, vector<int>(num_anchors, 0));

    // Iterate over the synteny blocks and update the graph edges
    for (int i = 0; i < blocks.size(); i++) {
        // Sort the anchor points in the block by their index
        sort(blocks[i].begin(), blocks[i].end(), [&](int a, int b) {
            return anchor_indices[a] < anchor_indices[b];
        });

        // Update the graph edges based on the block's anchor points
        for (int j = 0; j < blocks[i].size() - 1; j++) {
            int u = anchor_indices[blocks[i][j]];
            int v = anchor_indices[blocks[i][j+1]];
            graph[u][v]++;
            graph[v][u]++;
        }
    }

    // Initialize the order of the anchor points to their original indices
    vector<int> order(num_anchors);
    for (int i = 0; i < num_anchors; i++) {
        order[i] = i;
    }

    // Sort the anchor points by their degree in the graph
    sort(order.begin(), order.end(), [&](int a, int b) {
        int degree_a = 0;
        int degree_b = 0;
        for (int i = 0; i < num_anchors; i++) {
            degree_a += graph[a][i];
            degree_b += graph[b][i];
        }
        return degree_a > degree_b;
    });

    return order;
}

int main() {
    // Example usage
    vector<int> anchors = {1, 3, 5, 7};
    vector<vector<int>> blocks = {{1, 2, 3}, {3, 4, 5}, {5, 6, 7}};
    vector<int> order = anchor_align(anchors, blocks);
    for (int i = 0; i < order.size(); i++) {
        cout << anchors[order[i]] << " ";
    }
    cout << endl;
    return 0;
}
