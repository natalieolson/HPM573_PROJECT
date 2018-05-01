import Hookworm_ParameterClasses as P
import Hookworm_MarkovModel as MarkovCls
import scr.SamplePathClasses as PathCls
import scr.FigureSupport as Figs

# create and cohort
cohort = MarkovCls.Cohort(
    id=0,
    therapy=P.Therapies.ANNUAL)

simOutputs = cohort.simulate()

# graph infection curve
PathCls.graph_sample_path(
    sample_path=simOutputs.get_infection_curve(),
    title='infection curve',
    x_label='Simulation time step',
    y_label='Number of infected patients')

# graph histogram of infection durations
Figs.graph_histogram(
    data=simOutputs.get_infection_durations(),
    #data=[1, 2, 4, 4, 2],
    title='infection times of patients',
    x_label='infection time (years)',
    y_label='Counts',
    bin_width=2
)

# graph histogram of number of infections
Figs.graph_histogram(
    data=simOutputs.get_if_developed_infection(),
    #data=[1, 5, 2, 3, 5],
    title='Number of infections per Patient',
    x_label='infections',
    y_label='Counts',
    bin_width=1
)

Figs.graph_histogram(
    data=simOutputs.get_if_treated(),
    title='Number of Infections Treated',
    x_label = 'Treatments',
    y_label="counts",
    bin_width = 1
)
