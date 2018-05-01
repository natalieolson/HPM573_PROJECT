import Hookworm_ParameterClasses as P
import Hookworm_MarkovModel as MarkovCls
import Hookworm_SupportMarkov as SupportMarkov

# ANNUAL TREATMENT
# create a cohort
cohort_annual = MarkovCls.Cohort(id=0, therapy=P.Therapies.ANNUAL)
simOutputs_annual = cohort_annual.simulate()

# SEMIANNUAL
# create a cohort
cohort_semi = MarkovCls.Cohort(id=1, therapy=P.Therapies.SEMI)
simOutputs_semi = cohort_semi.simulate()

# draw survival curves and histograms
SupportMarkov.draw_infection_curves_and_histograms(simOutputs_annual, simOutputs_semi)

# print the estimates
SupportMarkov.print_outcomes(simOutputs_annual, "Annual MDA")
SupportMarkov.print_outcomes(simOutputs_semi, "Semi-Annual MDA")

# print comparative outcomes
SupportMarkov.print_comparative_outcomes(simOutputs_annual, simOutputs_semi)

# report the CEA results
SupportMarkov.report_CEA_CBA(simOutputs_annual, simOutputs_semi)

