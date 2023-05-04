#include <iostream>
#include <vector>
#include <unordered_map>

using namespace std;
/*
This implementation takes as input the vector of anchor points anchors and the vector of synteny blocks blocks, where each synteny block is represented as a vector of anchor points. It returns a vector order that contains the indices of the anchor points in the order determined by the Anchor-Align algorithm.

The implementation uses an unordered map anchor_index to map each anchor point to its index in the anchors vector. It then initializes a vector visited to keep track of which anchor points have already been visited, and a variable pos to keep track of the current position in the order.

The algorithm then iterates over each anchor point, starting from the first, and performs a breadth-first search to determine the order of the remaining anchor points. At each step, it finds the unvisited anchor point that is closest to the last visited anchor point in the current synteny block and adds it to the order.

Finally, the implementation returns the vector order containing the indices of the anchor points in the determined order.

Note that this is a simple implementation and may not be optimal for large-scale genome comparisons. More efficient implementations may use data structures such as priority queues or dynamic programming to improve performance.
*/
vector<int> anchor_align(vector<int>& anchors, vector<vector<int>>& blocks) {
    int n = anchors.size();
    unordered_map<int, int> anchor_index;
    for (int i = 0; i < n; i++) {
        anchor_index[anchors[i]] = i;
    }
    vector<int> order(n);
    vector<bool> visited(n);
    int pos = 0;
    for (int i = 0; i < n; i++) {
        if (!visited[i]) {
            visited[i] = true;
            order[pos++] = i;
            int last = anchors[i];
            while (pos < n) {
                int max_len = -1;
                int next = -1;
                for (int j = 0; j < n; j++) {
                    if (!visited[j]) {
                        int len = -1;
                        for (int k = 0; k < blocks[last].size(); k++) {
                            int index = anchor_index[blocks[last][k]];
                            if (index == j) {
                                len = k;
                                break;
                            }
                        }
                        if (len > max_len) {
                            max_len = len;
                            next = j;
                        }
                    }
                }
                visited[next] = true;
                order[pos++] = next;
                last = anchors[next];
            }
        }
    }
    return order;
}

vector<int> anchor_align(vector<int>& anchors, vector<vector<int>>& blocks);

void test_anchor_align() {
    vector<int> anchors = {1, 3, 5, 7};
    vector<vector<int>> blocks = {{1, 2, 3}, {3, 4, 5}, {5, 6, 7}};
    vector<int> expected_order = {0, 1, 2, 3};
    vector<int> actual_order = anchor_align(anchors, blocks);
    if (actual_order == expected_order) {
        cout << "Test case passed." << endl;
    } else {
        cout << "Test case failed." << endl;
        cout << "Expected order: ";
        for (int i = 0; i < expected_order.size(); i++) {
            cout << expected_order[i] << " ";
        }
        cout << endl;
        cout << "Actual order: ";
        for (int i = 0; i < actual_order.size(); i++) {
            cout << actual_order[i] << " ";
        }
        cout << endl;
    }
}

int main() {
    test_anchor_align();
    return 0;
}
