#ifndef RELATIVE_FUNCTIONS_H_INCLUDED
#define RELATIVE_FUNCTIONS_H_INCLUDED 1

#include <stdio.h>
#include <float.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_fit.h>
#include <err.h>
#include "dependencies/pcg_basic.h"
#include "sharedfunc_flag.h"
#include "main.h"
#include <tskit.h>
#include <tskit/tables.h>
#include <kastore.h>
#include <tskit/core.h>
#include <tskit/trees.h>

double RunSimulationRel(int tskitstatus, bool isabsolute, bool ismodular, int elementsperlb, char * Nxtimestepsname, char * popsizename, char * delmutratename, char * chromsizename, char * chromnumname, char * mubname, char * Sbname, int typeofrun, int Nxtimesteps, int popsize, int chromosomesize, int numberofchromosomes, double deleteriousmutationrate, double beneficialmutationrate, double Sb, int beneficialdistribution, double Sd, int deleteriousdistribution, gsl_rng * randomnumbergeneratorforgamma, FILE *miscfilepointer, FILE *veryverbosefilepointer, int rawdatafilesize, long double *haploidProductOfVmaxs, double michaelisConstant, int recessivityRunFlag, double initializationValRel, int addNonNeutral);

void PerformOneTimeStepRel(int tskitstatus, bool isabsolute, int isburninphaseover, bool ismodular, int elementsperlb, tsk_table_collection_t *treesequencetablecollection, tsk_id_t * wholepopulationnodesarray, tsk_id_t * wholepopulationsitesarray, int popsize, int totaltimesteps, double currenttimestep, long double *wholepopulationwistree, long double *wholepopulationwisarray, double *wholepopulationgenomes, long double * psumofwis, int chromosomesize, int numberofchromosomes, int totalindividualgenomelength, double deleteriousmutationrate, double beneficialmutationrate, double Sb, int beneficialdistribution, double Sd, int deleteriousdistribution, double *parent1gamete, double *parent2gamete, gsl_rng * randomnumbergeneratorforgamma, FILE *miscfilepointer, long double haploidProductOfVmaxs, double michaelisConstant, int recessivityRunFlag, double initializationValRel, int addNonNeutral);

void InitializePopulationRel(int tskitstatus, tsk_table_collection_t * treesequencetablecollection, tsk_id_t * wholepopulationnodesarray, tsk_id_t * wholepopulationsitesarray, long double *wholepopulationwistree, long double *wholepopulationwisarray, int popsize, double *wholepopulationgenomes, int totalpopulationgenomelength, int totaltimesteps, long double * psumofwis, double initializationValRel, long double haploidProductOfVmaxs, double michaelisConstant, int chromosomeSize, int recessivityRunFlag, int numberOfChromosomes);

int ChooseVictim(int populationsize);
int ChooseParentWithTree(long double *wholepopulationwistree, int popsize, long double sumofwis, FILE *miscfilepointer);

double CalculateWi(double *parent1gamete, double *parent2gamete, int totalindividualgenomelength, long double haploidProductOfVmaxs, double michaelisConstant, int recessivityRunFlag, double initializationValRel);
double FindFitnessContribution(double linkageBlockActivity, long double haploidProductOfVmaxs, double michaelisConstant);
double CalculateWiWithRecessivity(double *parent1gamete, double *parent2gamete, int totalindividualgenomelength, long double haploidProductOfVmaxs, double michaelisConstant, double initializationValRel);
double CalculateWiNoRecessivity(double *parent1gamete, double *parent2gamete, int totalindividualgenomelength);

void captureGeneActivitiesAndFitnessContributions(double *wholepopulationgenomes, char *captureFilename, int totalindividualgenomelength, int totalpopulationgenomelength, int popsize, long double haploidProductOfVmaxs, double michaelisConstant);
void captureSingleGeneActivityAndFitnessContribution(double *wholepopulationgenomes, char *captureFilename, int totalindividualgenomelength, int totalpopulationgenomelength, int popsize, long double haploidProductOfVmaxs, double michaelisConstant);
void renormalizeFitness(double *wholepopulationgenomes, long double *wholepopulationwisarray, long double *sortedwisarray, long double *psumofwis, long double *wholepopulationwistree, int totalpopulationgenomelength, long double *haploidProductOfVmaxs, double michaelisConstant, int popsize, FILE *miscfilepointer, int totalindividualgenomelength);

#endif // RELATIVE_FUNCTIONS_H_INCLUDED

