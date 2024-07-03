#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

struct Anchor {
    int query_pos;
    int subject_pos;
    int length;
    double score;
};

struct SyntenyBlock {
    int query_begin;
    int subject_begin;
    int length;
    double score;
};

vector<Anchor> find_anchors(const string& query, const string& subject, int anchor_length, int anchor_spacing) {
    vector<Anchor> anchors;
    for (int i = 0; i + anchor_length <= query.size(); i += anchor_spacing) {
        string anchor_query = query.substr(i, anchor_length);
        int anchor_subject_pos = subject.find(anchor_query);
        if (anchor_subject_pos != string::npos) {
            anchors.push_back({i, anchor_subject_pos, anchor_length, 0.0});
        }
    }
    return anchors;
}

double compute_anchor_score(const Anchor& anchor, const string& query, const string& subject, int gap_open_penalty, int gap_extend_penalty, int mismatch_penalty) {
    double score = 0.0;
    int query_pos = anchor.query_pos;
    int subject_pos = anchor.subject_pos;
    int anchor_length = anchor.length;
    for (int i = 0; i < anchor_length; ++i) {
        if (query[query_pos] == subject[subject_pos]) {
            score += 1.0;
        } else {
            score += mismatch_penalty;
        }
        ++query_pos;
        ++subject_pos;
    }
    return score;
}

vector<SyntenyBlock> find_synteny_blocks(const string& query, const string& subject, const vector<Anchor>& anchors, int gap_open_penalty, int gap_extend_penalty, int mismatch_penalty) {
    vector<SyntenyBlock> blocks;
    for (int i = 0; i < anchors.size(); ++i) {
        for (int j = i + 1; j < anchors.size(); ++j) {
            const Anchor& anchor1 = anchors[i];
            const Anchor& anchor2 = anchors[j];
            int query_anchor_distance = anchor2.query_pos - anchor1.query_pos - anchor1.length;
            int subject_anchor_distance = anchor2.subject_pos - anchor1.subject_pos - anchor1.length;
            if (query_anchor_distance == subject_anchor_distance) {
                double score = compute_anchor_score(anchor1, query, subject, gap_open_penalty, gap_extend_penalty, mismatch_penalty)
                             + compute_anchor_score(anchor2, query, subject, gap_open_penalty, gap_extend_penalty, mismatch_penalty);
                int length = query_anchor_distance + 2 * anchor1.length;
                blocks.push_back({anchor1.query_pos, anchor1.subject_pos, length, score});
            }
        }
    }
    return blocks;
}

bool compare_synteny_blocks(const SyntenyBlock& a, const SyntenyBlock& b) {
    return a.subject_begin < b.subject_begin;
}

vector<SyntenyBlock> filter_synteny_blocks(const vector<SyntenyBlock>& blocks, int min_length, double min_score) 
{
    vector<SyntenyBlock> filtered_blocks;
    for (const SyntenyBlock& block : blocks) {
        if (block.length >= min_length && block.score >= min_score * block.length) {
            filtered_blocks.push_back(block);
        }
    }
    return filtered_blocks;
}

vector<pair<int, int>> compute_matches(const vector<SyntenyBlock>& blocks) 
{
    vector<pair<int, int>> matches;
    for (const SyntenyBlock& block : blocks)
    {
        matches.push_back({block.query_begin, block.subject_begin});
    }
    sort(matches.begin(), matches.end());
    return matches;
}

int main() {
string query = "ACGTACGTAGTACGTAGTACGTAGTACGTAGTACGT";
string subject = "ACGTACGTAGTACGTAGTACGTAGTACGTAGTACGT";
int anchor_length = 5;
int anchor_spacing = 5;
int gap_open_penalty = 2;
int gap_extend_penalty = 1;
int mismatch_penalty = 2;
int min_length = 10;
double min_score = 0.5;vector<Anchor> anchors = find_anchors(query, subject, anchor_length, anchor_spacing);
vector<SyntenyBlock> blocks = find_synteny_blocks(query, subject, anchors, gap_open_penalty, gap_extend_penalty, mismatch_penalty);
vector<SyntenyBlock> filtered_blocks = filter_synteny_blocks(blocks, min_length, min_score);
vector<pair<int, int>> matches = compute_matches(filtered_blocks);

for (const auto& match : matches) {
    cout << "Match: " << match.first << "-" << match.first + anchor_length << " (" << query.substr(match.first, anchor_length) << ")"
         << " vs " << match.second << "-" << match.second + anchor_length << " (" << subject.substr(match.second, anchor_length) << ")" << endl;
}

return 0;
}