import scr.SamplePathClasses as PathCls
import scr.StatisticalClasses as StatCls
import scr.RandomVariantGenerators as rndClasses
import scr.EconEvalClasses as EconCls
import Hookworm_ParameterClasses as P
import Hookworm_InputData as Data

# patient class simulates patient, patient monitor follows patient, cohort simulates a cohort,
#  cohort outcome extracts info from simulation and returns it back


class Patient:
    def __init__(self, id, parameters):
        """ initiates a patient
        :param id: ID of the patient
        :param parameters: parameter object
        """
        self._infectionTime = []
        self._id = id
        # random number generator
        self._rng = None
        # parameters
        self._param = parameters
        # state monitor
        self._stateMonitor = PatientStateMonitor(parameters)
        # simulate time step
        self._delta_t = parameters.get_delta_t() # length of time step

    def simulate(self, sim_length):
        """ simulate the patient over the specified simulation length """
        # random number generator for this patient
        self._rng = rndClasses.RNG(self._id)  # from now on use random number generator from support library

        k = 0  # current time step

        # while the patient is alive and simulation length is not yet reached
        while k*self._delta_t < sim_length:
            # find transition probabilities of future state
            trans_prob = self._param.get_transition_prob(self._stateMonitor.get_current_state())
            # create an empirical distribution
            empirical_dist = rndClasses.Empirical(trans_prob)
            # sample from the empirical distribution to get a new state
            # (return an intger from {0, 1, 2, ...}
            new_state_index = empirical_dist.sample(self._rng) # pass RNG

            # update health state
            self._stateMonitor.update(k, P.HealthStats(new_state_index))

            # increment time step
            k += 1

    def get_infection_duration(self):
        """ returns the patient's infection time"""
        return self._stateMonitor.get_infection_duration()

    def get_number_of_infections(self):
        """ returns the patient's time to the POST_STROKE state """
        return self._stateMonitor.get_number_of_infections()

    def get_number_of_treated(self):
        return self._stateMonitor.get_number_of_treated()

    def get_total_discounted_cost(self):
        return self._stateMonitor.get_total_discounted_cost()

    def get_total_discounted_utility(self):
        return self._stateMonitor.get_total_discounted_utility()


class PatientStateMonitor:
    """ to update patient outcomes (years survived, cost, etc.) throughout the simulation """
    def __init__(self, parameters):
        """
        :param parameters: patient parameters
        """
        # current health state
        self._currentState = parameters.get_initial_health_state()
        self._delta_t = parameters.get_delta_t()
        self._transmissionTime = 0
        self._infectionTime = 0
        self._ifDevelopedInfection = False
        self._infectioncount = 0
        self._numberTreated = 0

        self._costUtilityOutcomes = PatientCostUtilityMonitor(parameters)

    def update(self, k, next_state):
        """
        :param k: current time step
        :param next_state: next state
        """

        # update infection count
        if self._currentState == P.HealthStats.INFECTED:
            self._ifDevelopedInfection= True
            self._infectioncount += 1

        #update treatment count
        if self._currentState ==P.HealthStats.TREATMENT:
            self._numberTreated +=1


        self._costUtilityOutcomes.update(k, self._currentState, next_state)

        self._currentState = next_state


        #update start of infection
        if self._currentState == P.HealthStats.WELL:
            if next_state == P.HealthStats.INFECTED:
                self._transmissionTime = (k+0.5) * self._delta_t

        # update infection time
        if self._currentState == P.HealthStats.INFECTED:
            if next_state in [P.HealthStats.WELL, P.HealthStats.TREATMENT]:
                self._infectionTime = (k+0.5) * self._delta_t - self._transmissionTime
                self._infectionTime.append((k+0.5) * self._delta_t - self._transmissionTime)


    def get_if_infected(self):
        if self._currentState == P.HealthStats.INFECTED:
            result = True
        if self._currentState in [P.HealthStats.WELL, P.HealthStats.TREATMENT]:
            result = False
            return result

    def get_current_state(self):
        return self._currentState

    def get_infection_duration(self):
        """ returns the patient infection time """
        # return infection time only if the patient has died
        return self._infectionTime

    def get_number_of_infections(self):
        return self._infectioncount

    def get_number_of_treated(self):
        return self._numberTreated

    def get_total_discounted_cost(self):
        return self._costUtilityOutcomes.get_total_discounted_cost()

    def get_total_discounted_utility(self):
        return self._costUtilityOutcomes.get_total_discounted_utility()


class PatientCostUtilityMonitor:

    def __init__(self, parameters):
        self._param = parameters
        self._totalDiscountedCost = 0
        self._totalDiscountedUtility = 0

    def update(self, k, current_state, next_state):

        # state cost and utility
        cost = 0.5*(self._param.get_annual_state_cost(current_state)
                    +(self._param.get_annual_state_cost(next_state))) * self._param.get_delta_t()

        utility = 0.5 * (self._param.get_annual_state_utility(current_state) +
                         (self._param.get_annual_state_utility(next_state))) * self._param.get_delta_t()

        # treatment cost (incurred only in post-stroke state)
        if current_state is P.HealthStats.TREATMENT:
            if next_state is P.HealthStats.WELL:
                cost += 0.5*self._param.get_annual_treatment_cost() * self._param.get_delta_t()
            else:
                cost += 1*self._param.get_annual_treatment_cost() * self._param.get_delta_t()

        self._totalDiscountedCost += EconCls.pv(cost, self._param.get_adj_discount_rate()/2, 2*k+1)
        self._totalDiscountedUtility += EconCls.pv(utility, self._param.get_adj_discount_rate()/2, 2*k+1)

    def get_total_discounted_cost(self):
        return self._totalDiscountedCost

    def get_total_discounted_utility(self):
        return self._totalDiscountedUtility


class Cohort:

    def __init__(self, id, therapy):
        """ create a cohort of patients
        :param id: an integer to specify the seed of the random number generator
        """
        self._initial_pop_size = Data.POP_SIZE
        self._patients = []      # list of patients

        # populate the cohort
        for i in range(self._initial_pop_size):
            # create a new patient (use id * pop_size + i as patient id)
            if Data.PSA_ON:
                patient = Patient(id* self._initial_pop_size + i, P.ParametersProbabilistic(i, therapy))
            else:
                patient = Patient(id * self._initial_pop_size + i, P._ParametersFixed(therapy))
            # add the patient to the cohort
            self._patients.append(patient)

    def simulate(self):
        """ simulate the cohort of patients over the specified number of time-steps
        :returns outputs from simulating this cohort
        """

        # simulate all patients
        for patient in self._patients:
            patient.simulate(Data.SIM_LENGTH)

        # return the cohort outputs
        return CohortOutputs(self)

    def get_initial_pop_size(self):
        return self._initial_pop_size

    def get_patients(self):
        return self._patients


class CohortOutputs:
    def __init__(self, simulated_cohort):
        """ extracts outputs from a simulated cohort
        :param simulated_cohort: a cohort after being simulated
        """

        self._infectionTimes = []        # patients' infection times
        self._count_infections = []
        self._utilities = []
        self._costs = []
        self._count_treated = []

        #infection curve
        self._infectionCurve = \
            PathCls.SamplePathBatchUpdate('Population size over time', id, simulated_cohort.get_initial_pop_size())

        # find patients' infection times
        for patient in simulated_cohort.get_patients():

            # get the patient infection time
            infection_time = patient.get_infection_duration()
            if not (infection_time is None):
                self._infectionTimes.append(infection_time)           # store the infection time of this patient
                self._infectionCurve.record(infection_time, -1)       # update the infection curve

            count_infections = patient.get_number_of_infections()
            count_treated = patient.get_number_of_infections()

            self._count_infections.append(count_infections)
            self._infectionTimes.append(infection_time)
            self._costs.append(patient.get_total_discounted_cost())
            self._utilities.append(patient.get_total_discounted_utility())
            self._count_treated.append(count_treated)

        # summary statistics
        self._sumStat_infectionTime = StatCls.SummaryStat('Patient infection time', self._infectionTimes)
        self._sumStat_number_infections = StatCls.SummaryStat('Time until infection', self._count_infections)
        self._sumStat_cost = StatCls.SummaryStat('Patient discounted cost', self._costs)
        self._sumStat_utility = StatCls.SummaryStat('Patient discounted utility', self._utilities)
        self._sumStat_treated= StatCls.SummaryStat('Number of Infections Treated', self._count_treated)

    def get_if_developed_infection(self):
        return self._count_infections

    def get_if_treated(self):
        return self._count_treated

    def get_infection_durations(self):
        return self._infectionTimes

    def get_sumStat_infection_times(self):
        return self._sumStat_infectionTime

    def get_infection_curve(self):
        return self._infectionCurve

    def get_sumStat_count_infections(self):
        return self._sumStat_number_infections

    def get_sumStat_count_treated(self):
        return self._sumStat_treated

    def get_costs(self):
        return self._costs

    def get_utilities(self):
        return self._utilities

    def get_sumStat_discounted_utility(self):
        return self._sumStat_utility

    def get_sumStat_discounted_cost(self):
        return self._sumStat_cost
