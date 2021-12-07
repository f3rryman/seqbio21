import argparse


class Node:
    """
    Node object with two parameters:
    pos = position in the text
    char2edge = dictionary of outgoing edges
    """
    def __init__(self, pos: int = -1, char2edge: dict = None):
        self.pos = pos
        self.char2edge = char2edge or {}

    def is_leaf(self):
        """
        Checks if the Node is a leaf, by its position pos.
        :return: True if pos != -1 (default)
        """
        return self.pos != -1


class Edge:
    """
    Edge object with two parameters:
    child = a Node to which the Edge points
    sequence = the corresponding sequence
    """
    def __init__(self, child: Node, sequence: str):
        self.child = child
        self.sequence = sequence

    def get_seq_len(self):
        """
        Calcs sequence length of the Edge.
        :return: len(sequences)
        """
        return len(self.sequence)

    def get_seq(self):
        """
        Get Sequence
        :return: self.sequence
        """
        return self.sequence

    def match_first_char(self, c: chr):
        """
        Checks if a provided char c matches with the first char in self.sequence.
        :param c: chr to check if first chr in sequence
        :return: True if c matches seq[0], False otherwise
        """
        if c == self.sequence[0]:
            return True
        else:
            return False

    def has_child(self):
        """
        Checks if Edge has a child.
        :return: True if Edge has child, False otherwise
        """
        if not self.child:
            return False
        else:
            return True


class SuffixTree:
    """
    Full suffix tree based on its root.
    """
    def __init__(self):
        self.root = None

    def compute_tree(self, text):
        """
        Computes a suffix tree from a text.
        :param text: String from which the suff tree is computed
        :return: Suffixtree based on its root
        """
        if text[-1:] != "$":
            text += "$"

        self.root = Node()
        root = self.root
        for i in range(len(text)):
            self.__add_suffix(root, i, i, text)
        #print(root.char2edge)

    def compute_common_prefix_length(self, s1, s2):
        """
        compute the shortest common prefix length of s1 and s2
        :param s1: string 1 to get common prefix from
        :param s2: string 2 to get common prefix from
        :return: int common prefix length
        """
        n = len(s1)
        m = len(s2)
        i = 0
        j = 0
        count = 0
        while i < n and j < m:
            if s1[i] == s2[j]:
                count += 1
            else:
                return count
            i += 1
            j += 1
        return count

    def __add_suffix(self, v: Node, start: int, pos: int, text: str):
        """
        Hidden method to add a suffix to the suffix tree
        :param v: actual Node in tree
        :param start: int, start position in the text
        :param pos: int, actual position in the text
        :param text: string, from which the tree is computed
        :return: void
        """
        current_char = text[pos]
        if not v.char2edge:
            edge = None
        else:
            edge = v.char2edge.get(current_char)

        if not edge:
            leaf = Node()
            leaf.pos = start
            edge = Edge(leaf, text[pos:])
            v.char2edge.update({current_char:edge})
        else:
            length = self.compute_common_prefix_length(text[pos:], edge.sequence)
            if length == edge.get_seq_len():
                self.__add_suffix(edge.child, start, pos + edge.get_seq_len(), text)
            else:
                topString = edge.sequence[0:length]
                bottomString = edge.sequence[length:]

                branch = Node()
                topEdge = Edge(branch, topString)
                v.char2edge.update({current_char: topEdge})
                bottomEdge = Edge(edge.child, bottomString)
                branch.char2edge.update({bottomEdge.sequence[0]: bottomEdge})
                self.__add_suffix(branch, start, pos + len(topString), text)

    def contains(self, query):
        """
        Checks if query is in the suffix tree.
        :param query: string, query to search for in the suffix tree.
        :return: False if first char in query can't be found from root, otherwise returns contains_rec
        """
        root = self.root
        if query[0] not in root.char2edge:
            return False
        else:
            return self.contains_rec(root, query)

    def contains_rec(self, v: Node, query):
        """
        Searches for the query further down of the root within the suffix tree.
        :param v: actual Node
        :param query: string to search for
        :return: True if query is matched, False if query is not in tree, Calls itself if its necessary to check
        leafs further down
        """
        if len(query) == 0:
            return True
        else:
            edge = v.char2edge[query[0]]
            if not edge:
                return False
            else:
                length = self.compute_common_prefix_length(edge.sequence, query)
                if length == len(query):
                    return True
                elif length == len(edge.sequence):
                    return self.contains_rec(edge.child, query[length:])
                else:
                    return False

    def find_all_occurrences(self, query):
        """
        Finds all occurrences of the query by recursively check all edges of the tree and returns its positions in the tree as list.
        :param query: string to search all occurrences in tree
        :return: list of int with positions of query occurrences in tree
        """
        root = self.root
        occur = []
        if root:
            self.find_all_occurrences_rec(root, query, occur)
        return occur

    def find_all_occurrences_rec(self, v: Node, query: str, list_of_occurrences: list):
        """
        recursion for find_all_occurrences
        :param v: actual Node
        :param query: query to look for
        :param list_of_occurrences: list of occurrences
        :return: void, fills list of occurrences
        """
        if len(query) > 0:
            edge = v.char2edge[query[0]]
            if edge:
                length = self.compute_common_prefix_length(edge.sequence, query)
                if length == len(query):
                    self.collect_all_suff_below_rec(edge.child, list_of_occurrences)
                elif length == len(edge.sequence):
                    self.find_all_occurrences_rec(edge.child, query[length:], list_of_occurrences)

    def collect_all_suff_below_rec(self, v: Node, list_of_occurrences: list):
        """
        Collects all suffixes below the actual node v and appends it to the list of occurences if Node v is a leaf.
        else checks childs of the Node recursively
        :param v: current Node
        :param list_of_occurences: list of occurences of a query were all leafes are appended
        :return: void, fills list of occurences
        """
        if v.is_leaf():
            list_of_occurrences.append(v.pos)
        else:
            for b in v.char2edge.values():
                below = b.child
                self.collect_all_suff_below_rec(below, list_of_occurrences)

    def print_suffix_tree(self):
        """
        Prints the simple suffix tree to console in a very basic manner. For debugging purposes.
        :return: void, prints tree to console
        """
        if self.root is None:
            print("Tree is empty.")
            return

        def show():
            for a, b in self.root.char2edge.items():
                print("id:", a, "pos:", b.child.pos, b.sequence)
                for c, d in b.child.char2edge.items():
                    print("\t","id:",  c, "pos:", d.child.pos, d.sequence)
                    for e, f in d.child.char2edge.items():
                        print("\t\t", "id:", e, "pos:", f.child.pos, f.sequence)
                        for g, h in f.child.char2edge.items():
                            print("\t\t\t", "id:", g, "pos:", h.child.pos, h.sequence)
        show()


def main():
    parser = argparse.ArgumentParser(description="Computes a suffix tree from provided text and checks for the "
                                                 "provided queries")
    parser.add_argument("text",
                        type=str,
                        help="Input text to create suffix tree from.")
    parser.add_argument("queries",
                        nargs="+",
                        type=str,
                        help="Multiple queries to search for in the suffix tree.")
    parser.add_argument("-T",
                        action="store_true",
                        help="If set prints the suffix tree to console in a very basic manner.",
                        default=False)
    args = parser.parse_args()
    text = args.text
    queries = args.queries
    print_tree = args.T

    tree = SuffixTree()
    tree.compute_tree(text)
    print("text: " + text)
    if print_tree:
        tree.print_suffix_tree()
    for q in queries:
        c = tree.contains(q)
        occ = tree.find_all_occurrences(q)
        print("Is query \"{}\" in tree? {}\n {}".format(q, c, occ))


if __name__ == "__main__":
    main()
