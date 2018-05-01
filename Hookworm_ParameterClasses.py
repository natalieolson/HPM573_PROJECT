from enum import Enum
import numpy as np
import scipy.stats as stat
import Hookworm_InputData as Data
import scr.MarkovClasses as MarkovCls
import scr.RandomVariantGenerators as Random
import scr.FittingProbDist_MM as Est


class HealthStats(Enum):
    """ health states of patients with HIV """
    WELL = 0
    INFECTED = 1
    TREATMENT = 2

class Therapies(Enum):
    """ mono vs. combination therapy """
    ANNUAL = 0
    SEMI = 1

class _Parameters:
    def __init__(self, therapy):
    # selected therapy
        self._therapy = therapy

        # simulation time step
        self._delta_t = Data.DELTA_T

        # calculate the adjusted discount rate
        self._adjDiscountRate = Data.DISCOUNT*Data.DELTA_T

        self._initialHealthState = HealthStats.WELL

        # annual treatment cost
        if self._therapy == Therapies.ANNUAL:
            self._annualTreatmentCost = Data.MDA_COST

        if self._therapy == Therapies.SEMI:
            self._annualTreatmentCost = Data.MDA_COST*2

# transition probability matrix of the selected therapy
        self._rate_matrix = []
        self._prob_matrix = []
        self._treatmentRR = 0

        self._annualStateCosts = []
        self._annualStateUtilities = []


    def get_initial_health_state(self):
        return self._initialHealthState

    def get_delta_t(self):
        return self._delta_t

    def get_adj_discount_rate(self):
        return self._adjDiscountRate

    def get_transition_prob(self, state):
        return self._prob_matrix[state.value]

    def get_annual_state_cost(self, state):
        if state == HealthStats.WELL:
            return 0
        else:
            return self._annualStateCosts[state.value]

    def get_annual_state_utility(self, state):
            return self._annualStateUtilities[state.value]

    def get_annual_treatment_cost(self):
        return self._annualTreatmentCost


class _ParametersFixed(_Parameters):

    def __init__(self, therapy):

                #initialize base class
        _Parameters.__init__(self,therapy)

        if therapy == Therapies.ANNUAL:
            self._rate_matrix = Data.TRANS_MATRIX
            self._prob_matrix[:], p = MarkovCls.continuous_to_discrete(self._rate_matrix, Data.DELTA_T)
        else:
            self._rate_matrix = Data.TRANS_MATRIX_SEMI
            self._prob_matrix[:], p = MarkovCls.continuous_to_discrete(self._rate_matrix, Data.DELTA_T)

        self._annualStateCosts = Data.ANNUAL_STATE_COST
        self._annualStateUtilities = Data.ANNUAL_STATE_UTILITY


class ParametersProbabilistic(_Parameters):
    def __init__(self, seed, therapy):

        #initialize base class
        _Parameters.__init__(self,therapy)

        self._rng = Random.RNG(seed)    # random number generator to sample from parameter distributions
        self._infectionProbMatrixRVG = []  # list of dirichlet distributions for transition probabilities
        self._lnRelativeRiskRVG = None  # random variate generator for the natural log of the treatment relative risk
        self._annualStateCostRVG = []       # list of random variate generators for the annual cost of states
        self._annualStateUtilityRVG = []    # list of random variate generators for the annual utility of states

 #transition probabilities
      #  j = 0
        # for prob in Data.TRANS_MATRIX:
          #  self._infectionProbMatrixRVG.append(Random.Dirichlet(prob[j:]))
          #  j += 1

        # annual state cost
        for cost in Data.ANNUAL_STATE_COST:
            # find shape and scale of the assumed gamma distribution
            estDic = Est.get_gamma_params(mean=cost, st_dev=cost/4)
            # append the distribution
            self._annualStateCostRVG.append(
                Random.Gamma(a=estDic["a"], loc=0, scale=estDic["scale"]))

        # annual state utility
        for utility in Data.ANNUAL_STATE_UTILITY:
            # find alpha and beta of the assumed beta distribution
            estDic = Est.get_beta_params(mean=utility, st_dev=utility/4)
            # append the distribution
            self._annualStateUtilityRVG.append(
                Random.Beta(a=estDic["a"], b=estDic["b"]))

        # resample parameters
        self.__resample()

    def __resample(self):


       # self._prob_matrix = []
       # for s in HealthStats:
           # self._prob_matrix.append([0] * len(HealthStats))

        # for all health states
       # for s in HealthStats:
                # sample from the dirichlet distribution to find the transition probabilities between states
               # sample = self._infectionProbMatrixRVG[s.value].sample(self._rng)
               # for j in range(len(sample)):
                   # self._prob_matrix[s.value][s.value+j] = sample[j]


        # sample from gamma distributions that are assumed for annual state costs
        self._annualStateCosts = []
        for dist in self._annualStateCostRVG:
            self._annualStateCosts.append(dist.sample(self._rng))

        # sample from beta distributions that are assumed for annual state utilities
        self._annualStateUtilities = []
        for dist in self._annualStateUtilityRVG:
            self._annualStateUtilities.append(dist.sample(self._rng))

