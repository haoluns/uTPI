# uTPI (Utility-Based Toxicity Probability Interval Design)
R codes to implement the utility-based toxicity probability interval design for dose finding in phase I/II trials.
# Description
We propose a utility-based toxicity probability interval (uTPI) design for finding the optimal biological dose (OBD) in phase I/II trials. The proposed design accounts for both toxicity and efficacy outcomes, and uses quasi-binomial likelihood to simplify the modeling of utility. The dose-assigment decisions of the uTPI design are adaptively made by maximizing the toxicity--efficacy tradeoffs.
The uTPI design is model-assisted in nature, which simply models the utility outcomes observed at the current dose level based on a quasi binomial likelihood.  Toxicity probability intervals are used to screen out overly toxic dose levels, and  then the dose escalation/de-escalation decisions are  made adaptively by comparing the posterior utility distributions of the adjacent levels of the current dose. The uTPI design is  flexible in accommodating various utility functions but it only needs minimum design parameters. A prominent feature of the uTPI design is that its dose-assignment decision table can be pre-calculated before the start of trial, which greatly simplifies  practical implementation of the design. 
#Functions
The repository includes two functions:
* uTPI-sim-fixed.R: The R code that includes the function ```get.oc``` to conduct simulations and obtain the operating characteristics of the uTPI design under certain scenarios.
```rscript
get.oc(targetE,targetT, pE.true,pT.true,u1,u2,ncohort, cohortsize,cutoff.eliT, cutoff.eliE, nstar,ntrial)
```
* uTPI-decisiontable.R: The R code that includes the function ```get.decision.table``` to generate the decision table for selecting the next dose level for the new patents by the uTPI design.
```rscipt
get.decision.table(u1,u2,targetE,targetT,cutoff.eliE,cutoff.eliT,nstar)
```


#Inputs
* ```targetE```: The clinically uninteresting efficacy probability, e.g., ```targetE<-0.25```.
* ```targetT```: The target toxicity probability, e.g., ```targetT<-0.3```.
* ```pE.true```: A vector of length *J* that stores the efficacy probability for each dose, where *J* is the total number of dose levels.
* ```pT.true```: A vector of length *J* that stores the toxicity probability for each dose, where *J* is the total number of dose levels.
* ```u1```: A value between 0 and 100 that represent the utility of the bivariate outcome of both toxicity and efficacy, e.g., ```u1<-70```.
* ```u2```: A value between 0 and 100 that represent the utility of the bivariate outcome of both non-toxicity and non-efficacy, e.g., ```u2<-30```.
* ```ncohort```: The number of cohorts in the trial, e.g., ```ncohort<-12```.
* ```cohortsize```: The size of each cohort,  e.g., ```cohortsize<-3```.
* ```cutoff.eliT```: The probability cutoff cutoff for eliminating over-toxic doses, e.g., ```cutoff.eliT<-0.95```.
* ```cutoff.eliE```: The probability cutoff cutoff for eliminating futile doses, e.g., ```cutoff.eliE<-0.9```.
* ```nstar```: The cutoff sample size for aggresive dose exploration, as a default, ```nstar<-9```.
* ```ntrial```: Number of replicated trials in the simulation, e.g., ```ntrial<-10000```.


#Example
We apply the uTPI design to the B-cell epitope vaccine immunotherapy trial.
*  Given the target toxicity rate of 0.3 and efficacy rate of 0.25, and the utility strucuture is 0.7 and 0.3 for (toxicity,efficacy) and (non-toxicity,non-efficacy) outcomes, respectively. The probability cutoff cutoff for eliminating over-toxic doses and futile doses are 0.95 and 0.9, respectively. The decision table of the uTPI design based on the default setting can be pre-calculated.

```rscript
u1 = 70
u2 = 30
targetE = 0.25
targetT =0.30
cutoff.eliE = 0.90
cutoff.eliT = 0.95
get.decision.table(u1,u2,targetE,targetT,cutoff.eliE,cutoff.eliT)
```
The output is given by 
```rscript
     Num. Patients Num. Tox Num. Eff Tox Interval Desirability Score
 [1,] 0             0        0        0            40                
 [2,] 3             0        0        1            12                
 [3,] 3             0        1        1            36                
 [4,] 3             0        2        1            56                
 [5,] 3             0        3        1            81                
 [6,] 3             1        0        4            12                
 [7,] 3             1        1        4            36                
 [8,] 3             1        2        4            56                
 [9,] 3             1        3        4            81                
[10,] 3             2        0        7            12     
```
* To obtain the operating characteristics under certain prespecified scenario, we may utilize the ```get.noc``` function to generate the averaged operating metrics, as given by
```rscript 
cohortsize = 3
ncohort = 12
targetE=0.25
targetT=0.3
u1 = 70
u2 = 30
pT.true<-c(0.20,0.40,0.45,0.50,0.55)
pE.true<-c(0.40,0.50,0.60,0.70,0.80)

oc=get.oc(targetE,targetT, pE.true,pT.true,u1,u2, ncohort, cohortsize, cutoff.eliT=0.95, cutoff.eliE=0.90, nstar=9,ntrial=10)
print(oc)

-----------------------output------------------------
The selection percentages of all the doses are 
50 30  0  0  0
The average number of patients allocated at all the doses are 14.1
The average number of patients allocated at all the doses are 12.9
The average number of patients allocated at all the doses are 2.1
The average number of patients allocated at all the doses are 0.9
The average number of patients allocated at all the doses are 0
14.1 12.9  2.1  0.9  0.0
The average number of efficacy outcomes at all the doses are 
3.2 4.3 1.1 0.3 0.0
The average number of toxicity outcomes at all the doses are 
6.5 6.8 1.1 0.4 0.0
The percentage of early stopped trials is 20
The average total number of of toxicity outcomes is 8.9
The average total number of of toxicity outcomes is 14.8
The selection percentage of the best dose is 50
The selection percentage of the optimal dose is 50
The number of overdoses is 15.9
The percentage of poor allocation is 40

```

# Authors and Reference
* Haolun Shi, Jiguo Cao, Ying Yuan, and Ruitao Lin
* Shi, H., Cao, J., Yuan, Y., Lin, R. (2019) uTPI: A Utility-Based Toxicity Probability Interval Design for Dose Finding in Phase I/II Trials.
