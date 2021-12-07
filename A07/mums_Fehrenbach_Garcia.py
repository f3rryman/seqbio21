import argparse

from suffixtree_Fehrenbach_Garcia import SuffixTree

class suffixTreeMums(SuffixTree):
    """
    Extension of the class suffixTree to build a tree of two strings and computing all MUMs
    """
    def find_all_mums(self, text1: str, text2: str, min_length: int = 1):
        """
        Builds a suffix tree of the two input strings text1%text2$
        :param text1: string, input to compute suffix tree from
        :param text2: string, input to compute suffix tree from
        :param min_length: int, minimal length L of the MUMs
        :return: out = list of ints, positions of the MUMs within the texts,
        text1[out[0]:out[1]+1], text2[out[2]:out[3]+1]
        """
        concat_text = text1 + "%" + text2 + "$"
        mid = len(text1)
        self.compute_tree(concat_text)
        #suff_tree = SuffixTree.compute_tree(self, concat_text)
        #print(suff_tree)
        out = []
        self.find_all_mums_rec(self.root, concat_text, mid, 0, min_length, out)
        return out

    def find_all_mums_rec(self, v, concat_text: str, middle: int, depth: int, min_length: int, out):
        """
        recursion to find all MUMs within the suffix tree
        :param v: current Node
        :param concat_text: string, concatenate text1%text2$
        :param middle: int, middle of text1
        :param depth: int, depth within the suffix tree
        :param min_length: int, minimal length L of the MUMs
        :param out: list of ints, positions of the MUMs within the texts,
        text1[out[0]:out[1]+1], text2[out[2]:out[3]+1]
        :return: void, fills the output list out
        """
        count_leaf_childs = 0
        count_internal_childs = 0
        for val in v.char2edge.values():
            child = val.child
            if child.is_leaf():
                count_leaf_childs += 1
            else:
                self.find_all_mums_rec(child, concat_text, middle, depth + len(val.sequence), min_length, out)
                count_internal_childs += 1

        if depth >= min_length and count_leaf_childs == 2 and count_internal_childs == 0:
            children = list(v.char2edge.values())
            left_child = children[0].child
            right_child = children[1].child
            if left_child.pos < middle <= right_child.pos:
                if (left_child.pos == 0 or
                        right_child.pos == 0 or
                        concat_text[left_child.pos - 1] != concat_text[right_child.pos - 1]):
                    start1 = min(left_child.pos, right_child.pos)
                    start2 = max(left_child.pos, right_child.pos) - middle - 1
                    out.append([start1, start1 + depth - 1, start2, start2 + depth - 1])


def main():
    parser = argparse.ArgumentParser(description="Takes two strings as input, "
                                                 "computes the suffix tree of both and calculates all MUMs. "
                                                 "This works only for two texts!")
    parser.add_argument("text1",
                        type=str,
                        help="First text from which the suffix tree is calculated.")
    parser.add_argument("text2",
                        type=str,
                        help="Second text from which the suffix tree is calculated.")
    parser.add_argument("-L",
                        type=int,
                        help="Defines the minimal length of the MUMs. Default == 1.",
                        default=1)
    args = parser.parse_args()
    text1 = args.text1
    text2 = args.text2
    min_length = args.L

    mums = suffixTreeMums().find_all_mums(text1, text2, min_length)

    print("text 1: " + text1)
    print("text 2: " + text2)
    if len(mums) == 0:
        print("No MUMs of length >= L = {}".format(min_length))
    else:
        print("MUMs:")
        for i in range(len(mums)):
            s1 = mums[i][0]
            e1 = mums[i][1]
            mum1 = text1[s1:e1+1]
            s2 = mums[i][2]
            e2 = mums[i][3]
            mum2 = text2[s2:e2+1]
            print("{}-{} ({}) - {}-{} ({})".format(s1, e1, mum1, s2, e2, mum2))
    #print(mums)


if __name__ == "__main__":
    main()
