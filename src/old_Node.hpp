// class Node;
#include <iostream>
using std::ostream;

class Node {
public:
  // Node data
  int sample_ID;
  int site;
  bool genotype;

  // Node pointers
  Node* above;
  Node* below;
  Node* left;
  Node* right;
  Node* next_0;
  Node* next_1;
  Node* prev_0;
  Node* prev_1;

  // Constructors
  Node(int sample_ID, int site, bool genotype);

  // Lookers-uppers
  Node* first_uncle(bool g);
  Node* first_nephew(bool g);

  // Node movers-arounders
  void insert_above(Node* node);
  void insert_below(Node* node);
  void insert_right(Node* node);
  void insert_left(Node* node);

  // Recursive 0/1-updaters
  void update_prev_1_recursive(Node* old_node, Node* new_node);
  void update_next_1_recursive(Node* old_node, Node* new_node);
  void update_prev_0_recursive(Node* old_node, Node* new_node);
  void update_next_0_recursive(Node* old_node, Node* new_node);

  // Comparators
  bool equivalent_to(Node* other);
  friend bool operator==(const Node& node1, const Node& node2) {
    return (node1.sample_ID == node2.sample_ID && node1.site == node2.site) && node1.genotype == node2.genotype;
  }

  // Output
  friend ostream& operator<<(ostream& os, const Node& node);
};