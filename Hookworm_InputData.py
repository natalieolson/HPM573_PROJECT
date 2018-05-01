import numpy as np

# simulation settings
POP_SIZE = 500   # cohort population size
SIM_LENGTH = 10 # length of simulation (years)
ALPHA = 0.05        # significance level for calculating confidence intervals
DISCOUNT = 0.03     # daily discount rate

DELTA_T = 1/12

PSA_ON = False

#Transmission Rate (Rate of infection S-->I)
a = -np.log(1-(6.5/100))

#Natural Recovery Rate
b = 1/2

#Probability of Treatment(coverage)
c = 0.562 * 1

#Probability of Effective Treatment

d = 0.784 * 52

e = 52 * (1-0.784)

# transition matrix
TRANS_MATRIX = [
    [None,           a,        0],   #Susceptible
    [b,           None,           c],   #Infected
    [d,           e,        None] #Treatment
 ]

TRANS_MATRIX_SEMI = [
    [None,        a,          0],   #Susceptible
    [b,        None,           2*c],   #Infected
    [d,        e,          None] #Treatment
 ]


# annual cost of each health state
ANNUAL_STATE_COST = [
    0,   # Susceptible
    0,   # Infected
    .50,   #Treatment
    ]

# annual health utility of each health state
ANNUAL_STATE_UTILITY = [
    1.00,   # Susceptible
    0.902,   # Infected
    1.00,    #Treatment
    ]

# annual drug costs
MDA_COST = 0.5 #per person

RR_TREAT = 0.562
RR_RECOVERY = 0.784
