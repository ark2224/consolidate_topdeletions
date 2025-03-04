import numpy as np
import heapq

# from sequence_manipulation import ATGC_energy


def test_1():
    seq1 = 'CAAAAAGGAGCTCTCTC'
    seq2 = 'AAATTGGGGAATCTCTGC'
    # print(find_path(seq1,seq2))
    # print(find_graph(seq1,seq2))
    p, e, s = find_graph(seq1,seq2)
    print(s)
    process_output_path(p)
    print(e)
    print()
    # ans: AAAGGGCTCTCTC
    return False


def test_2():
    seq1 = 'Iwonderifthiswillwork'
    seq2 = 'maybeitwillmaybeitwont'
    # print(find_path(seq1,seq2))
    # print(find_graph(seq1,seq2))
    p, e, s = find_graph(seq1,seq2)
    print(s)
    process_output_path(p)
    print(e)
    print()
    # ans: AAAGGCTCTCTC
    return False


def test_3():
    # seq1455_4 deletion 38-153
    # After deletion begining/end
    seq1 = 'AAGAGAGAAAAGAAGAGTAA'
    seq2 = 'ACCCAGAGCCCCAGCAGCCT'
    # print(find_path(seq1,seq2))
    # print(find_graph(seq1,seq2))
    p, e, s = find_graph(seq1,seq2)
    print(s)
    process_output_path(p)
    print(e)
    print()

    # Before deletion beginning/end
    seq1 = 'ACGACACACTATAAGGAAAT'
    seq2 = 'CCAGACCCGCCATCCTGATG'
    # print(find_path(seq1,seq2))
    # print(find_graph(seq1,seq2))
    p, e, s = find_graph(seq1,seq2)
    print(s)
    process_output_path(p)
    print(e)
    print()
    return False


def find_prev_max(mat: np.array, mismatch_penalty=1.):
    '''Input:
        - mat: array to find max closest to bottom right corner (presumably (i,j) in dp)
        - mismatch_penalty: how to penalize jumps when creating string
    Return
        - index into mat (relative to an array (i,j) in size) 
    '''
    i = len(mat)
    j = len(mat[0])
    max = 0
    max_idx = [-1,-1]
    for i_ in range(i-2, -1, -1):
        for j_ in range(j-2, -1, -1):
            curr_score = mat[i_,j_] - mismatch_penalty * (min(i-i_-1 + j-j_-1, 0))
            if max < curr_score:
                max = curr_score
                max_idx = [i_, j_]
    return max_idx


def ATGC_energy(s1, mismatch_penalty=1.5) -> int:
    score = 0
    for c in s1:
        if c == 'A' or c == 'T':
            score += 2
        elif c == 'G' or c == 'C':
            score += 3
        else:
            score -= mismatch_penalty
    return score


# def find_path(s1: str, s2: str):
#     '''Inputs:
#         - s1: sequence outside deletion (corresponds to d1)
#         - s2: sequence inside deletion (corresponds to d2)
#     Returns:
#         - d1: distance from deletion (outside)
#         - d2: distance from deletion (inside)
#         - energy
#     '''
#     m = len(s1)
#     n = len(s2)
#     dp = np.zeros((m,n))
#     paths = {}
#     for i in range(m):
#         for j in range(n):
#             paths[str(i) + ',' + str(j)] = ['', 0, 0]#[path, d1, d2]
#     d1 = 0
#     d2 = 0
#     e = 0
#     seq = ''
#     for i in range(m):
#         for j in range(n):
#             next_step = ''
#             if s1[i] == s2[j]:
#                 dp[i,j] += 2 if s1[i] == 'A' or s1[i] == 'T' else 3
#                 next_step = s1[i]
#             else:
#                 continue

#             i_,j_ = find_prev_max(dp[:i+1,:j+1])
#             if i_ >= 0:
#                 dp[i,j] += dp[i_, j_]
#             else:
#                 continue
#             paths[str(i) + ',' + str(j)][0] = paths[str(i_)+','+str(j_)][0] + next_step
#             if paths[str(i) + ',' + str(j)][0] == next_step:
#                 paths[str(i) + ',' + str(j)][1] = i
#                 paths[str(i) + ',' + str(j)][2] = j
#             else:
#                 paths[str(i) + ',' + str(j)][1:] = paths[str(i_)+','+str(j_)][1:]
#             if e < dp[i,j]:
#                 e = dp[i,j]
#                 seq, d1, d2 = paths[str(i) + ',' + str(j)]
#     # now how the fuck are you getting d1/d2?
#     print(dp)

#     return d1, d2, e, seq

# Finding repeats instead of individual bases:
def find_path(s1: str, s2: str):
    '''Inputs:
        - s1: sequence outside deletion (corresponds to d1)
        - s2: sequence inside deletion (corresponds to d2)
    Returns:
        - d1: distance from deletion (outside)
        - d2: distance from deletion (inside)
        - energy
    '''
    m = len(s1)
    n = len(s2)
    dp = np.zeros((m,n))
    paths = {}
    for i in range(m):
        for j in range(n):
            paths[str(i) + ',' + str(j)] = ['', 0, 0]#[path, d1, d2]
    d1 = 0
    d2 = 0
    e = 0
    seq = ''
    for i in range(1,m):
        for j in range(1,n):
            next_step = ''
            if s1[i] == s2[j] and s1[i-1] == s2[j-1]:
                dp[i,j] += 2 if s1[i] == 'A' or s1[i] == 'T' else 3
                if dp[i-1,j-1] < 2:
                    dp[i,j] += 2 if s1[i-1] == 'A' or s1[i-1] == 'T' else 3
                    next_step = s1[i-1:i+1]
                else:
                    dp[i,j] += dp[i-1,j-1]#continuing a repeat through the diagonal
                    next_step = s1[i]
            else:
                continue

            i_,j_ = find_prev_max(dp[:i,:j])
            if i_ >= 0:
                dp[i,j] += dp[i_, j_]
            else:
                continue
            paths[str(i) + ',' + str(j)][0] = paths[str(i_)+','+str(j_)][0] + next_step
            if paths[str(i) + ',' + str(j)][0] == next_step:
                paths[str(i) + ',' + str(j)][1] = i
                paths[str(i) + ',' + str(j)][2] = j
            else:
                paths[str(i) + ',' + str(j)][1:] = paths[str(i_)+','+str(j_)][1:]
            if e < dp[i,j]:
                e = dp[i,j]
                seq, d1, d2 = paths[str(i) + ',' + str(j)]
    return d1, d2, e, seq


# graph solution
class node(object): 
    def __init__(self, seq='', d1_=0,d2_=0, en=0):
        self.edges = []
        self.edge_weights = []
        self.sequence = seq
        self.d1 = d1_
        self.d2 = d2_
        self.energy = en

    def __iter__(self):
        return {
            self.sequence,
            self.d1,
            self.d2
        }
    
    def __lt__(self, other):
        if isinstance(other, node):
            return self.energy < other.energy
    def __le__(self, other):
        if isinstance(other, node):
            return self.energy <= other.energy
    def __gt__(self, other):
        if isinstance(other, node):
            return self.energy > other.energy
    def __ge__(self, other):
        if isinstance(other, node):
            return self.energy >= other.energy
    
    def add(self, n):
        dist1 = self.d1 + len(self.sequence) - n.d1
        dist2 = self.d2 + len(self.sequence) - n.d2
        if dist1 <= 0 and dist2 <= 0:
            self.edges.append(n)
            self.edge_weights.append((dist1 + dist2) * 0.5)#mismatch_penalty, not an average


def process_output_path(path):
    returnable_path = {}
    prev = node('',0,0)
    for i,n in enumerate(path):
        m = {
            f'length_{i}': len(n.sequence),
            f'energy_{i}': n.energy,
            f'sequence_{i}': n.sequence,
            f'gap{i}_d1': n.d1 - len(prev.sequence) - prev.d1,
            f'gap{i}_d2': n.d2 - len(prev.sequence) - prev.d2,
        }
        returnable_path.update(m)
        prev = n
    # for k,v in returnable_path.items():
    #     print(k,v)
    return returnable_path

def find_graph(s1: str, s2: str):
    m = len(s1)
    n = len(s2)
    # ============MAKING GRAPH============
    # create nodes
    nodes = [node('',0,0)]
    for i in range(m):
        for j in range(n):
            idx = 0
            while idx+i<m and idx+j<n and s1[i+idx] == s2[j+idx]:
                idx += 1
            if idx > 1:
                nodes.append(
                    node(s1[i:i+idx], i, j, ATGC_energy(s1[i:i+idx]))
                )
    nodes = sorted(nodes, key=lambda x: min(x.d1+len(x.sequence), x.d2+len(x.sequence)))#sort by min to prevent overlaps as you move forward
    # populate edges
    for i in range(len(nodes)):
        for j in range(i+1, len(nodes)):
            nodes[i].add(nodes[j])
    
    heap = [[0, nodes[0], list(), '']]
    heapq.heapify(heap)
    visited = set()

    maxPath = []
    pathSequence = ''
    maxEnergy = 0
    while heap:
        curr_score, n, path, seq = heapq.heappop(heap)
        if curr_score < maxEnergy:
            maxEnergy = curr_score
            maxPath = path
            pathSequence = seq
        visited.add(n)

        for i, no in enumerate(n.edges):
            if no not in visited:
                next_score = curr_score - n.edge_weights[i] - no.energy
                next_path = path.copy()
                next_path.append(no)
                next_seq = seq + '|' + no.sequence

                heapq.heappush(heap, [next_score, no, next_path, next_seq])
    return maxPath, -maxEnergy, pathSequence
        



                
# return:
#   match1_length, 
#   match1_energy,
#   match1_sequence,
#   gap1_d1,
#   gap1_d2,

#   match2_length,
#   match2_energy
#   match2_sequence
#   gap2_d1
#   gap2_d2...

# def main():
#     test_1()
#     test_2()
#     test_3()
    # try graphing solution tomorrow
    # return


# if __name__ == '__main__':
#     main()