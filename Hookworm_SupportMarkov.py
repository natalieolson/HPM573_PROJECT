import Hookworm_InputData as Settings
import scr.FormatFunctions as F
import scr.StatisticalClasses as Stat
import scr.EconEvalClasses as Econ
import scr.SamplePathClasses as PathCls
import scr.FigureSupport as Figs


def print_outcomes(simOutput, therapy_name):
    """ prints the outcomes of a simulated cohort
    :param simOutput: output of a simulated cohort
    :param therapy_name: the name of the selected therapy
    """
    # mean and confidence interval text of patient infection time
    infection_mean_CI_text = F.format_estimate_interval(
        estimate=simOutput.get_sumStat_infection_times().get_mean(),
        interval=simOutput.get_sumStat_infection_times().get_t_CI(alpha=Settings.ALPHA),
        deci=2)

    # mean and confidence interval text of time to infection
    infection_mean_CI_text = F.format_estimate_interval(
        estimate=simOutput.get_sumStat_count_infections().get_mean(),
        interval=simOutput.get_sumStat_count_infections().get_t_CI(alpha=Settings.ALPHA),
        deci=2)

    treatment_mean_CI_text = F.format_estimate_interval(
        estimate = simOutput.get_sumStat_count_treated().get_mean(),
        interval = simOutput.get_sumStat_count_treated().get_t_CI(alpha=Settings.ALPHA),
        deci=2)

    cost_mean_CI_text = F.format_estimate_interval(
        estimate=simOutput.get_sumStat_discounted_cost().get_mean(),
        interval=simOutput.get_sumStat_discounted_cost().get_t_CI(alpha=Settings.ALPHA),
        deci=2)

    utility_mean_CI_text = F.format_estimate_interval(
        estimate=simOutput.get_sumStat_discounted_utility().get_mean(),
        interval=simOutput.get_sumStat_discounted_utility().get_t_CI(alpha=Settings.ALPHA),
        deci=2)

    # print outcomes
    print(therapy_name)
    print("  Estimate of mean and {:.{prec}%} CI of infection time:".format(1 - Settings.ALPHA, prec=0),
          infection_mean_CI_text)
    print("  Estimate of mean and {:.{prec}%} CI of time to infection:".format(1 - Settings.ALPHA, prec=0),
          infection_mean_CI_text)
    print("  Estimate of discounted cost and {:.{prec}%} CI:".format(1 - Settings.ALPHA, prec=0),
          cost_mean_CI_text)
    print("  Estimate of discounted utility and {:.{prec}%} CI:".format(1 - Settings.ALPHA, prec=0),
          utility_mean_CI_text)
    print("")


def draw_infection_curves_and_histograms(simOutputs_ANNUAL, simOutputs_SEMI):
    """ draws the infection curves and the histograms of time until HIV deaths
    :param simOutputs_mono: output of a cohort simulated under mono therapy
    :param simOutputs_combo: output of a cohort simulated under combination therapy
    """

    # histograms of infection times
    set_of_infection_times = [
        simOutputs_ANNUAL.get_infection_durations(),
        simOutputs_SEMI.get_infection_durations()
    ]

    # graph histograms
    Figs.graph_histograms(
        data_sets=set_of_infection_times,
        title='Histogram of Infection Duration',
        x_label='Infection Duration (Year)',
        y_label='Counts',
        bin_width=1,
        legend=['Annual Treatment', 'Semi-Annual Treatment'],
        transparency=0.6
    )



    #get infection curves of both treatments
    infection_curves = [
       simOutputs_ANNUAL.get_infection_curve(), simOutputs_SEMI.get_infection_curve()]

    # graph infection curve
    PathCls.graph_sample_paths(
        sample_paths=infection_curves,
        title='infection curve',
        x_label='Simulation time step (year)',
        y_label='Number of infected patients',
        legends=['ANNUAL', '6-month'])


def print_comparative_outcomes(simOutputs_ANNUAL, simOutputs_SEMI):
    """ prints average increase in infection time, discounted cost, and discounted utility
    under combination therapy compared to mono therapy
    :param simOutputs_mono: output of a cohort simulated under mono therapy
    :param simOutputs_combo: output of a cohort simulated under combination therapy
    """


    decrease_infection_time = Stat.DifferenceStatIndp(name="decrease in infection time",
                                                     x=simOutputs_SEMI.get_infection_durations(),
                                                     y_ref=simOutputs_ANNUAL.get_infection_durations())

    estimate_CI = F.format_estimate_interval(estimate=decrease_infection_time.get_mean(),
                                             interval=decrease_infection_time.get_t_CI(alpha=Settings.ALPHA),
                                             deci=2)

    estimate_CI =F.format_estimate_interval(
        estimate=decrease_infection_time.get_mean(),
        interval=decrease_infection_time.get_t_CI(alpha=Settings.ALPHA),
        deci=2)
    print("Average decrease in infection duration "
          "and {:.{prec}%} CI:".format(1 - Settings.ALPHA, prec=0),
          estimate_CI)

    # increase in discounted total cost under combination therapy with respect to mono therapy
    increase_discounted_cost = Stat.DifferenceStatIndp(
        name='Increase in discounted cost',
        x=simOutputs_SEMI.get_costs(),
        y_ref=simOutputs_ANNUAL.get_costs())

    # estimate and CI
    estimate_CI = F.format_estimate_interval(
        estimate=increase_discounted_cost.get_mean(),
        interval=increase_discounted_cost.get_t_CI(alpha=Settings.ALPHA),
        deci=0,
        form=F.FormatNumber.CURRENCY)
    print("Average increase in discounted cost "
          "and {:.{prec}%} CI:".format(1 - Settings.ALPHA, prec=0),
          estimate_CI)

    # increase in discounted total utility under combination therapy with respect to mono therapy
    increase_discounted_utility = Stat.DifferenceStatIndp(
        name='Increase in discounted cost',
        x=simOutputs_SEMI.get_utilities(),
        y_ref=simOutputs_ANNUAL.get_utilities())

    # estimate and CI
    estimate_CI = F.format_estimate_interval(
        estimate=increase_discounted_utility.get_mean(),
        interval=increase_discounted_utility.get_t_CI(alpha=Settings.ALPHA),
        deci=2)
    print("Average increase in discounted utility "
          "and {:.{prec}%} CI:".format(1 - Settings.ALPHA, prec=0),
          estimate_CI)


def report_CEA_CBA(simOutputs_ANNUAL, simOutputs_SEMI):
    annual_MDA=Econ.Strategy(
        name="Annual MDA",
        cost_obs=simOutputs_ANNUAL.get_costs(),
        effect_obs=simOutputs_ANNUAL.get_utilities())

    semiannual_MDA=Econ.Strategy(
        name="Semi-Annual MDA",
        cost_obs=simOutputs_SEMI.get_costs(),
        effect_obs=simOutputs_SEMI.get_utilities())

    listofStrategies = [annual_MDA, semiannual_MDA]

    CEA = Econ. CEA(listofStrategies, if_paired=False)

    CEA.show_CE_plane(
        title='Cost-Effectiveness Analysis',
        x_label='Additional discounted utility',
        y_label='Additional discounted cost',
        show_names=True,
        show_clouds=True,
        show_legend=True,
        figure_size=6,
        transparency=0.3
    )
    # report the CE table
    CEA.build_CE_table(
        interval=Econ.Interval.CONFIDENCE,
        alpha=Settings.ALPHA,
        cost_digits=2,
        effect_digits=2,
        icer_digits=2,
    )

    CBA = Econ.CBA(listofStrategies, if_paired=False)

    NBA = Econ.CBA(
        strategies = [annual_MDA, semiannual_MDA],
        if_paired= False
    )

    NBA.graph_deltaNMB_lines(
        min_wtp=0,
        max_wtp=900,
        x_label="Willingness-to-pay for one unit reduction in DALY ($)",
        y_label="Incremental Net Monetary Benefit ($)",
        interval=Econ.Interval.CONFIDENCE,
        transparency=0.4,
        show_legend=True,
        figure_size=6,
        title='Cost Benefit Analysis')
