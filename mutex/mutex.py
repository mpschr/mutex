__author__ = 'mpschr'

import pandas as pd
import numpy as np



class MutExResult(object):
    def __init__(self, coverage, signal, higher_coverage_count, lower_coverage_count, n, mean_sim_coverage,
                 sample_size):
        self.sample_size = sample_size
        self.mean_sim_coverage = mean_sim_coverage
        self.n = n
        self.lower_coverages = lower_coverage_count
        self.higher_coverages = higher_coverage_count
        self.signal = signal
        self.coverage = coverage

        self.mutex_pvalue = higher_coverage_count / n
        self.co_occurence_pvalue = lower_coverage_count / n

    def __str__(self):
        return "MuTexResult\n  Mutual Exclusive p-value:\t{}\n  Co-occurence p-value:\t\t{}\n  permutations:\t\t\t\t{}".format(
            self.mutex_pvalue, self.co_occurence_pvalue, self.n
        )

    def __repr__(self):
        self.__str__()

class MutEx(object):
    def __init__(self, background: pd.DataFrame, n: int=100):
        """
        :param background: A data frame containing all the observations as binary data 1 and 0 or True and False where
                rows represent observations and columns represent samples.
        :param n: how many permutations by default
        :return:
        """
        self.n = n
        self.background = background
        self.sample_weights = background.apply(sum) / background.apply(sum).pipe(sum)
        self.cummulative_sum = np.cumsum(self.sample_weights)
        self.sample_indices = [x for x in range(0,background.shape[1])]


    def calculate(self, indices: list, use_empty_samples = True, n = None) -> MutExResult:
        """
        :param indices: A list of indices for which to test the MutEx. The indices refer the the background-data row-ids.
        :return: MutExResult
        """

        if not all([x in self.background.index for x in indices]):
            raise Exception("Not all indices found in background")

        target = self.background.loc[indices]

        coverage = target.apply(max).pipe(sum)
        observation_signal = target.apply(sum, axis=1)
        signal = sum(observation_signal)

        if n == None:
            n = self.n

        sim_coverages = []
        higher_coverage = []
        lower_coverage = []
        print("running {} permutations".format(n))

        for i in range(0, n):

            sim = self.__simulate(observation_signal)

            sim_cov = sim.apply(max).pipe(sum)
            sim_coverages.append(sim_cov)

            higher_coverage.append(sim_cov >= coverage)
            lower_coverage.append(sim_cov <= coverage)

        return MutExResult(coverage=coverage, signal=signal,
                           higher_coverage_count=np.sum(higher_coverage),
                           lower_coverage_count=np.sum(lower_coverage), n=n,
                           mean_sim_coverage=np.mean(sim_coverages),
                           sample_size=len(self.sample_weights)
                           )



    def __simulate(self, observation_signal):
        simulations = []
        for observations in observation_signal:
            simulations.append(self.__weighted_choice(observations))
        return pd.DataFrame.from_records(simulations).fillna(0)

    def __weighted_choice(self, amount: int):
        return {x:1 for x in  np.random.choice(self.sample_indices, amount, False, self.sample_weights)}




def test():
    """

    :rtype : None
    """

    import scipy.sparse as sparse


    row, col = 100, 100
    np.random.seed(77)
    df = pd.DataFrame(sparse.random(100,100, density=0.2).A).apply(np.ceil)

    df.loc[0] = [1 if x < 20 else 0 for x in range(0,df.shape[1]) ]
    df.loc[1] = [1 if x > 13 and x < 35 else 0 for x in range(0,df.shape[1]) ]
    df.loc[2] = [1 if x > 80 else 0 for x in range(0,df.shape[1]) ]

    m = MutEx(background=df, n=100)

    pd.set_option('display.max_columns', 100)
    print(df.loc[[0,1,2]])
    print(m.calculate([0, 1, 2]))


    print(m.calculate([4, 5, 6]))



if __name__ == "__main__":
    test()