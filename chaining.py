import math
import logging

'''
一个用于描述序列段（种子）的类
'''


class Anchor:
    def __init__(self, read_offset, chromosome, position, node=None, offset=None, is_reverse=False):
        # 读取偏移量
        self.read_offset = read_offset
        # 所在染色体的编号
        self.chromosome = chromosome
        # 开始位置
        self.position = position
        # 是否逆序
        self.is_reverse = is_reverse
        # 节点值
        self.node = node
        # 偏移量
        self.offset = offset

    '''
    判断另一个序列段有可能匹配到当前序列构成的序列链之后
    :param other_anchor 考察的序列段
    :return False表示不可能匹配，True表示有可能
    '''

    def fits_to_left_of(self, other_anchor):
        # TODO: If both on reverse strand, left means right
        if other_anchor.chromosome != self.chromosome:
            # 不属于同一个染色体，不考虑是否匹配
            return False

        if other_anchor.is_reverse != self.is_reverse:
            # 两个序列必须顺序一致
            return False

        if other_anchor.position < self.position:
            # 被匹配的序列段必须从当前序列的位置之后开始
            return False

        # Check that read offsets kind of match
        if abs((other_anchor.read_offset - self.read_offset) - (other_anchor.position - self.position)) > 64:
            # 两个序列段的长度差不能超过64
            return False

        return True

    '''
    内置方法，实现序列段的小于号比较
    :param other_anchor 待比较的序列段
    :return True表示小于，False表示不小于
    '''

    def __lt__(self, other_anchor):
        if self.chromosome < other_anchor.chromosome:
            # 所在染色体编号的大小决定了序列段的大小
            return True
        elif self.chromosome > other_anchor.chromosome:
            return False
        else:
            if self.position < other_anchor.position:
                # 处于同一染色体时，序列段的开始位置决定了序列段的大小
                return True
            else:
                return False

    '''
        比较两个序列段
        :param other_anchor 待比较的序列段
        :return 处于同一染色体时,返回开始位置的差，否则返回染色体编号的差
    '''

    def __compare__(self, other_anchor):
        # TODO: IF on reverse strand, reverse comparison on position level
        if other_anchor.chromosome != self.chromosome:
            return other_anchor.chromosome - self.chromosome
        else:
            return other_anchor.position - self.position

    '''
        内置方法，实现对象的字符串表达
        :return 打印出对象的各属性的字符串值
    '''

    def __str__(self):
        return "%s:%d/%d, %s:%s (%d)" % (
        self.chromosome, self.position, self.is_reverse, self.node, self.offset, self.read_offset)

    '''
        内置方法，实现对象的查看
        :return 打印出对象的字符串表达
    '''

    def __repr__(self):
        return self.__str__()

    '''
        内置方法，实现序列段的等于号比较
        :param other 待比较的序列段
        :return True表示两个序列段的所在染色体和开始位置一致，否则为False
    '''

    def __eq__(self, other):
        return self.chromosome == other.chromosome and self.position == other.position and \
               self.is_reverse == other.is_reverse


'''
    用于描述序列链的类
'''


class Chain:
    def __init__(self):
        # 存储内部的序列段的数组
        self.anchors = []
        # 得分，初始为-1
        self.chaining_score = -1
        #
        self.mapq = None

    '''
    判定序列链是否为空
    :return True表示内部的序列段数组长度为0，否则为False
    '''

    def is_empty(self):
        return len(self.anchors) == 0

    '''
        在序列链里添加一个序列段
        :param anchor 待添加的序列段
    '''

    def add_anchor(self, anchor):
        self.anchors.append(anchor)

    '''
        判断在序列链里添加某一个序列段是否合适。
        判定依据是必须在同一个染色体且离第一个序列段不能超过给定最大距离
        :param anchor 待添加的序列段
        :param max_distance_allowed 容许的最大距离，缺省为100
        :return True表示可以添加，否则为False
    '''

    def anchor_fits_with_chain(self, anchor, max_distance_allowed=100):
        first_anchor = self.anchors[0]
        if anchor.chromosome != first_anchor.chromosome:
            return False
        elif anchor.position > first_anchor.position + max_distance_allowed:
            return False
        else:
            return True

    '''
        内置方法，实现对象的字符串表达
        :return 打印出各个序列段以---连接，最后打印得分
    '''

    def __str__(self):
        return ' --- '.join([str(a) for a in self.anchors]) + ", score: %.2f" % self.chaining_score

    '''
        内置方法，实现对象的查看
        :return 打印出对象的字符串表达
    '''

    def __repr__(self):
        return self.__str__()

    '''
        内置方法，实现序列段的遍历，在for ... in 时返回每一个序列段
    '''

    def __iter__(self):
        for anchor in self.anchors:
            yield anchor

    '''
        内置方法，实现序列链的大小比较
        :param other 待比较的序列链
        :return 返回各自最后一个序列段的比较结果
    '''

    def __lt__(self, other):
        return self.anchors[-1] < other.anchors[-1]

    '''
        比较方法，实现序列链的比较
        :param other 待比较的序列链
        :return 返回各自最后一个序列段的比较结果
    '''

    def __compare__(self, other):
        return self.anchors[-1].__compare__(other)

    '''
        内置方法，实现序列链的长度输出
        :return 内部序列段数组的长度
    '''

    def __len__(self):
        return len(self.anchors)


'''
用于描述序列链处理器的类
'''


class Chainer:
    def __init__(self, anchors, mean_seed_length=21, w=21, max_anchor_distance=150):
        # 平均序列长度，缺省值为21
        self._mean_seed_length = mean_seed_length
        # 最佳得分下限
        self.w = w
        # 序列链的最大序列段距离
        self._max_anchor_distance = max_anchor_distance
        # 待处理序列段数组
        self.anchors = anchors
        # 分配好的序列链
        self.chains = []

    '''
    获取两个序列段之间的距离
    :param i 第一个序列段的编号
    :param j 第二个序列段的编号
    :param fast_mode 是否启用快速模式。未实现
    :return 两个序列段的读偏移量之差。如超过最佳得分下限，返回该下限值
    '''

    def distance_between_anchors(self, i, j, fast_mode=True):
        distance = min(self.anchors[i].read_offset - self.anchors[j].read_offset, self.w)
        return distance

    '''
    获取两个序列段之间的长度差
    :param i 第一个序列段的编号
    :param j 第二个序列段的编号
    :return 两个序列段的长度差。
    '''

    def lfunc(self, i, j):
        return ((self.anchors[i].position - self.anchors[j].position) -
                (self.anchors[i].read_offset - self.anchors[j].read_offset))

    '''
    获取两个序列段匹配是否可能存在gap的得分值。越小越好
    :param i 第一个序列段的编号
    :param j 第二个序列段的编号
    :param fast_mode 是否启用快速模式。如果启用快速模式，该得分为负值
    :return 得分值。如不属于同一染色体或位置超前、距离超大等情况下，得分超大；
                    否则返回一个正数
    '''

    def anchor_gap_cost(self, i, j, fast_mode=False):
        anchor1 = self.anchors[j]
        anchor2 = self.anchors[i]
        # print(anchor1.position, anchor2.position)
        if anchor1.chromosome != anchor2.chromosome:
            return 10e15

        if anchor1.position >= anchor2.position:
            # print("Anchor1 pos > anchor2 pos")
            return 10000000

        if anchor2.position - anchor1.position > self._max_anchor_distance or \
                anchor2.read_offset - anchor1.read_offset > self._max_anchor_distance:
            return 10000

        distance = abs(self.lfunc(i, j))
        if distance == 0:
            return 0
        elif distance > 150:
            return 1000000
        else:
            if fast_mode:
                if distance > self._mean_seed_length:
                    return -(self._mean_seed_length - 0.01 * distance * self._mean_seed_length)
                else:
                    return -distance

            else:
                return 0.01 * self._mean_seed_length * abs(distance) + 0.5 * math.log2(abs(distance))

    '''
    获取两个序列段匹配可能的得分。越大越好
    :param i 第一个序列段的编号
    :param j 第二个序列段的编号
    :return 得分值。如不属于同一染色体或位置超前、距离超大等情况下，得分超小；
                    否则返回一个负数
    '''

    def anchors_score(self, i, j):
        anchor1 = self.anchors[j]
        anchor2 = self.anchors[i]
        # print(anchor1.position, anchor2.position)
        if anchor1.chromosome != anchor2.chromosome:
            return -10e15

        if anchor1.position >= anchor2.position:
            # print("Anchor1 pos > anchor2 pos")
            return -10000000

        if anchor2.position - anchor1.position > self._max_anchor_distance or \
                anchor2.read_offset - anchor1.read_offset > self._max_anchor_distance:
            return -10000

        distance = self.lfunc(i, j)
        if distance == 0:
            return 0
        return min(distance, 21) - (0.01 * 21 * abs(distance) + 0.5 * math.log2(abs(distance)))

    '''
    构造序列链，存放到chains数组里
    :param minimium_chaining_score 最小得分，缺省值为20
    :param minimum_minimizers_on_chain 第小序列链长度，缺省值为2
    '''

    def get_chains(self, minimium_chaining_score=20, minimum_minimizers_on_chain=2):
        chaining_scores = [self.w]  # Maximal score up to the ith anchor at ith position
        best_predecessors = {0: -1}  # Best predecessor anchors of anchor at key
        # 为每一个序列段寻找其最优后续序列段和得分值
        # Formula 1 and 2 from Heng Li's Minimap2 paper
        for i in range(1, len(self.anchors)):
            # print("Finding score for i=%d" % i)
            scores = {}  # Scores (value) between chain j (key) and i
            for j in range(i - 1, -1, -1):  # All j's lower than i

                score = max(chaining_scores[j] + self.distance_between_anchors(i, j) - self.anchor_gap_cost(i, j),
                            self.w)
                # score = max(chaining_scores[j] + self.anchors_score(i, j), self.w)
                # print("   j=%d. Dist: %.3f. Score: %.2f. Gap cost: %d" % (j, self.distance_between_anchors(i, j), score, self.anchor_gap_cost(i, j)))
                scores[j] = score
                if j < i - 1:
                    # Using heuristic that chaining this far down probably gives a lower score
                    break
            # print("    scores: %s" % str(scores))
            # Find best predecessor as the one with max score
            best_predecessor, best_predecessor_score = max(scores.items(), key=lambda s: s[1] + 0.000001 * s[
                0])  # Important to prioritize highest key if equal score
            best_predecessor_score = scores[best_predecessor]
            if best_predecessor_score == self.w:
                best_predecessor = -1
            # print("    Best predecessor of %d: %d" % (i, best_predecessor))
            best_predecessors[i] = best_predecessor
            chaining_scores.append(best_predecessor_score)
        # print(score)
        # Backtracking
        # print(best_predecessors)
        chains = []
        # print(" == Backtracking ==")
        # 将所有序列段分到不同的序列链里，保证每一个序列链里前后序列段都是最佳后续关系
        # 并且序列链的长度和得分满足给定的要求
        used_anchors = set()
        for i in range(len(self.anchors) - 1, -1, -1):
            # print("Checking anchor %d" % i)
            if i in used_anchors:
                continue
            current_chain = Chain()
            current_anchor = i

            while True:
                used_anchors.add(current_anchor)
                current_chain.anchors.append(self.anchors[current_anchor])
                # Find best predecessor
                best = best_predecessors[current_anchor]
                # print("  Found best %d" % best)
                if best == -1 or best in used_anchors:
                    break
                current_anchor = best
            current_chain.chaining_score = chaining_scores[i]
            if current_chain.chaining_score >= minimium_chaining_score and \
                    len(current_chain.anchors) >= minimum_minimizers_on_chain:
                # Hack for now, we need to reverse order to be compatible with old setup
                current_chain.anchors = current_chain.anchors[::-1]
                chains.append(current_chain)
        # print(chains)

        # if len(chains) > 100:
        #    print("\n".join([str(c) for c in chains]))

        # print("=========")
        # print("\n".join([str(c) for c in chains]))
        self.chains = chains


'''
测试函数
'''


def test_close_anchors():
    anchors = [
        Anchor(20, 1, 10),
        Anchor(30, 1, 20),
        Anchor(32, 1, 22),
        Anchor(42, 1, 32),
    ]

    chainer = Chainer(anchors)
    chainer.get_chains(minimium_chaining_score=0)
    chains = chainer.chains
    print("\n".join(str(c) for c in chains))


if __name__ == "__main__":
    test_close_anchors()
