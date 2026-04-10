#include <iostream>
#include <vector>
#include <fstream>
#include <queue>
#include <list>

using namespace std;
std::ofstream order("order.txt",std::ios::out);

// Structure to represent an edge
struct Edge {
    int to;
    int weight;
};

// Structure for priority queue
typedef pair<int, int> pii; // {weight, node}

void preorder(int u, int t, vector<vector<pair<int, int>>>& adj, vector<bool>& visited) {
    visited[u] = true;

    if(u!=t) order << u << " ";
//    std::cout << " u: " << u << " t: " << t;
    for (auto& edge : adj[u]) {
        int v = edge.first;
        if (!visited[v]) {
//	    std::cout << "v: " << v;
            preorder(v, t, adj, visited);
        }
    }
}

int mst_preorder(int s, int t) {
    int nodes, edges;
    ifstream file("spanner.txt");
    if (!file.is_open()) return 1;

    file >> nodes >> edges;
  //  std::cout << nodes << " " << edges << std::endl;
    vector<vector<pair<int, int>>> adj(nodes + 1); // Original Graph
    vector<vector<pair<int, int>>> mst_adj(nodes + 1); // MST Graph

    for (int i = 0; i < edges; ++i) {
        int u, v, w;
        file >> u >> v >> w;
//	std::cout << u << " " << v << " " << w << std::endl;
        adj[u].push_back({v, w});
        adj[v].push_back({u, w});
    }
    
    std::cout << "spanner read complete." << std::endl;
    // Prim's Algorithm
    priority_queue<pii, vector<pii>, greater<pii>> pq;
    vector<int> key(nodes + 1, 1e9);
    vector<int> parent(nodes + 1, -1);
    vector<bool> inMST(nodes + 1, false);

    pq.push({0, s}); // Start from node 1
    key[1] = 0;

    int mst_weight = 0;
    while (!pq.empty()) {
        int u = pq.top().second;
        pq.pop();

        if (inMST[u]) continue;
        inMST[u] = true;
        mst_weight += key[u];

        // Store MST edge
        if (parent[u] != -1) {
            mst_adj[parent[u]].push_back({u, key[u]});
            mst_adj[u].push_back({parent[u], key[u]});
        }

        for (auto& edge : adj[u]) {
            int v = edge.first;
            int weight = edge.second;
            if (!inMST[v] && key[v] > weight) {
                key[v] = weight;
                parent[v] = u;
                pq.push({key[v], v});
            }
        }
    }

    cout << "Total MST Weight: " << mst_weight << endl;
    cout << "Preorder Traversal: ";
    vector<bool> visited(nodes + 1, false);
    preorder(s, t, mst_adj, visited);
    order << t;
    cout << std::endl << "Traversal Finished." << std::endl;
    order.close();
    return 0;
}
