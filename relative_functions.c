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
#include "relative_functions.h"
#include "sharedfunc_flag.h"
#include "main.h"
#include <tskit.h>
#include <tskit/tables.h>
#include <kastore.h>
#include <tskit/core.h>
#include <tskit/trees.h>

double RunSimulationRel(int tskitstatus, bool isabsolute, bool ismodular, int elementsperlb, char * Nxtimestepsname, char * popsizename, char * delmutratename, char * chromsizename, char * chromnumname, char * mubname, char * Sbname, int typeofrun, int Nxtimesteps, int popsize, int chromosomesize, int numberofchromosomes, double deleteriousmutationrate, double beneficialmutationrate, double Sb, int beneficialdistribution, double Sd, int deleteriousdistribution, gsl_rng * randomnumbergeneratorforgamma, FILE *miscfilepointer, FILE *veryverbosefilepointer, int rawdatafilesize, long double *maxRateOfReaction, double michaelisConstant, int recessivityRunFlag, double initializationValRel, int addNonNeutral)
{
    // Additional bottleneck parameters: , int shouldBottleneck, int bottleneckPopsize, int bottleneckLength, int postBottleneckGenerations

    if(isabsolute){
        fprintf(miscfilepointer, "\n Trying to use RunSimulationRel within an absolute fitness program. \n");
        exit(0);
    }
    
    FILE *rawdatafilepointer;
    FILE *summarydatafilepointer;
    //tree sequence data files
    FILE *nodefilepointer;
    FILE *edgefilepointer;
    FILE *sitefilepointer;
    FILE *mutationfilepointer;
    
    int i, j, k;
    
    char * rawdatafilename = (char *) malloc(200);
    strcpy(rawdatafilename, "rawdatafor"); //starting the string that will be the name of the data file.

    strcat(rawdatafilename, "Nxtimesteps"); //for adding values of generations to the data name.
    strcat(rawdatafilename, Nxtimestepsname);

    strcat(rawdatafilename, "popsize"); //for adding values of starting population sizes to the data name.
    strcat(rawdatafilename, popsizename);

    strcat(rawdatafilename, "mutrate"); //for adding values of mutation rate to the data name (remember that mutation rate is currently the per-locus rate, not per-genome).
    strcat(rawdatafilename, delmutratename);

    strcat(rawdatafilename, "chromsize"); //for adding values of chromosome size to the data name.
    strcat(rawdatafilename, chromsizename);

    strcat(rawdatafilename, "chromnum"); //for adding values of the number of chromosomes to the data name.
    strcat(rawdatafilename, chromnumname);
    
    strcat(rawdatafilename, "benmutrate"); //for adding values of the beneficial mutation rate to the data name.
    strcat(rawdatafilename, mubname);
    
    strcat(rawdatafilename, "Sb"); //for adding values of the beneficial mutation effect size to the data name.
    strcat(rawdatafilename, Sbname);
    
    strcat(rawdatafilename, ".txt");

    rawdatafilepointer = fopen(rawdatafilename, "w"); //opens the file to which to print data to be plotted.
    fprintf(rawdatafilepointer, "Nxtimesteps,Sum.of.wis,Variance.in.fitness,Two.times.standard.deviation.fitness,Mean.fitness\n");
    
    char * summarydatafilename = (char *) malloc(100);
    strcpy(summarydatafilename, "summarydatafor");
    strcat(summarydatafilename, "Sb");
    strcat(summarydatafilename, Sbname);
    strcat(summarydatafilename, "mub");
    strcat(summarydatafilename, mubname);
    strcat(summarydatafilename, ".txt");
    summarydatafilepointer = fopen(summarydatafilename, "w"); //opens the file to which to print summary data.
    
    nodefilepointer = fopen("nodetable.txt", "w");
    edgefilepointer = fopen("edgetable.txt", "w");
    sitefilepointer = fopen("sitetable.txt", "w");
    mutationfilepointer = fopen("mutationtable.txt", "w");
    
    int totaltimesteps = Nxtimesteps * popsize;
    double currenttimestep = 0.0;
    double *wholepopulationgenomes;
    int totalpopulationgenomelength;
    int totalindividualgenomelength;
    totalpopulationgenomelength = popsize * numberofchromosomes * 2 * chromosomesize;
    totalindividualgenomelength = numberofchromosomes * 2 * chromosomesize;
    wholepopulationgenomes = malloc(sizeof(double) * totalpopulationgenomelength);
    long double sumofwis;
    long double *psumofwis = &sumofwis;
    long double *wholepopulationwistree;
    wholepopulationwistree = malloc(sizeof(long double) * popsize);
    
    //Following lines are from the tskit library.
    //Initializes the tables that make up the tree sequence recording.
    tsk_table_collection_t treesequencetablecollection;
    tsk_table_collection_t * tablepointer = &treesequencetablecollection;
    int returnvaluefortskfunctions = tsk_table_collection_init(&treesequencetablecollection, 0);
    check_tsk_error(returnvaluefortskfunctions);
    
    tsk_id_t *wholepopulationnodesarray;
    wholepopulationnodesarray = malloc(sizeof(tsk_id_t) * 2 * popsize);
    //The extant nodes need to have explicit identification in order to add edges between parents and children nodes.
    //Each node is only a single set of chromosomes, so the 2 here assumes diploidy.
    
    tsk_id_t wholepopulationsitesarray[totalindividualgenomelength / 2];
    //The number of sites is the number of linkage blocks in a single set of chromosomes (haploid), and won't change over the course of the simulation.
       
    long double *wholepopulationwisarray;
    wholepopulationwisarray = malloc(sizeof(long double) * popsize);
    //The Fenwick tree does not store each individual's wi, but rather a collection of partial sums.
    //For debugging purposes and data that requires summarizing wis, storing the wis in an array is necessary.

    long double *sortedwisarray;
    sortedwisarray = malloc(sizeof(long double) * popsize);
    //In order to visualize the distribution of fitness in the population,
    //the individual fitnesses need to be sorted, which requires a separate array.

    if (VERYVERBOSE == 1) {
        fprintf(veryverbosefilepointer, "Entered simulation run.\n");
        fflush(veryverbosefilepointer);
    }
    
    InitializePopulationRel(tskitstatus, &treesequencetablecollection, wholepopulationnodesarray, wholepopulationsitesarray, wholepopulationwistree, wholepopulationwisarray, popsize, wholepopulationgenomes, totalpopulationgenomelength, totaltimesteps, psumofwis, initializationValRel, *maxRateOfReaction, michaelisConstant, chromosomesize, recessivityRunFlag, numberofchromosomes);
    
    /*Sets the initial population to have zeroes in all their linkage blocks,
    death rates equal to the baseline wi, and an identifier number.
    It also sums the wis and returns the sum. Note that, if a run with recessivity
    is being performed, then gene activities are initialized to the log of the 
    a0 value given at the command line*/
    
    if (VERYVERBOSE == 1) {
        fprintf(veryverbosefilepointer, "Population initialized.\n");
        fflush(veryverbosefilepointer);
    }
    
    double *logaveragefitnesseachNtimesteps;
    logaveragefitnesseachNtimesteps = malloc(sizeof(double) * Nxtimesteps);
    //In order to calculate the slope of degradation of fitness,
    //I need to store the average fitness each generation.
    
    long double currentfittestindividualswi;
    double parent1gamete[numberofchromosomes*chromosomesize], parent2gamete[numberofchromosomes*chromosomesize];
    
    //Following array is to store the variance in log(fitness) for twenty generations at a time for estimating the end of the burn-in phase.
    //The choice of twenty is completely arbitrary and could eventually be an input if it seems worth it.
    size_t step = 1;
    double *last200Ntimestepsvariance;
    double *literallyjustlast200Ntimesteps;
    literallyjustlast200Ntimesteps = malloc(sizeof(double) * 200);
    last200Ntimestepsvariance = malloc(sizeof(double) * 200);
    for (k = 0; k < 200; k++) {
        literallyjustlast200Ntimesteps[k] = 0.0;
        last200Ntimestepsvariance[k] = 0.0;
    }
    double slopeofvariance;
    int isburninphaseover = 0;
    int didpopulationcrash = 0;
    int endofburninphase;
    int endofdelay = Nxtimesteps-1;
    int endofsimulation = Nxtimesteps-1;
    int Nxtimestepsafterburnin = 0;
    double arbitrarynumber;
    arbitrarynumber = (-1 * 0.007 / popsize); //using a number somewhere close to the mean of the DFE for deleterious mutations.
    double slopeoflogfitness;    
    double variancesum;
    
    if (VERYVERBOSE == 1) {
        fprintf(veryverbosefilepointer, "Variables initialized, preparing to begin simulation.\n");
        fflush(veryverbosefilepointer);
    }

    // Determines the generations in which we should capture information about the genome, for the 
    // purposes of ensuring that the fitness functionality is working properly
    int numberOfCheckpoints = 1000;
    int generationModForCheckpoints = Nxtimesteps / numberOfCheckpoints; 

    // Declares the char arrays that will contain the current filenames to write to, both for capturing the 
    // average activity for each gene across all individuals, and for capturing the gene activity and fitness 
    // contribution for a single gene (the first one) for all individuals in the population 
    char checkpointFilename[] = "checkpointFileActivitiesContributions_gen";
    char singleCheckpointFilename[] = "checkpointFileSingleIndividualActivitiesContributions_gen";
    
    //BEGIN THE SIMULATION FOR LOOP 
    for (i = 0; i < Nxtimesteps; i++) {

        // Checks to see wheter this is a generation in which we should re-normalize fitness
        if (i % RENORM_FREQUENCY == 0 && i != 0) {
            // If so, we call our renormalization method
            renormalizeFitness(wholepopulationgenomes, wholepopulationwisarray, sortedwisarray, psumofwis, wholepopulationwistree, totalpopulationgenomelength, maxRateOfReaction, michaelisConstant, popsize, miscfilepointer, totalindividualgenomelength);
        }

        // Checks to see whether this is a generation in which we should capture information on the genome,
        // capturing this information if desired
        if (i % generationModForCheckpoints == 0 && i != 0) {
            // Constructs the filenames to actually write to
            char currentCheckpointFilename[90];
            char currentSingleCheckpointFilename[90];
            sprintf(currentCheckpointFilename, "%s%d.txt", checkpointFilename, i + 1);
            sprintf(currentSingleCheckpointFilename, "%s%d.txt", singleCheckpointFilename, i + 1);

            // Writes to this filename via our helper method
            captureGeneActivitiesAndFitnessContributions(wholepopulationgenomes, currentCheckpointFilename, totalindividualgenomelength, totalpopulationgenomelength, popsize, *maxRateOfReaction, michaelisConstant);
            captureSingleGeneActivityAndFitnessContribution(wholepopulationgenomes, currentSingleCheckpointFilename, totalindividualgenomelength, totalpopulationgenomelength, popsize, *maxRateOfReaction, michaelisConstant);
        }
        
        //Following code performs N rounds of paired births and deaths.
        for (j = 0; j < popsize; j++) {
            if (i > 1000) {
                fprintf(miscfilepointer, "Current sum of wis: %Lf\n", *psumofwis);
                fprintf(miscfilepointer, "Current Vmax: %Lf\n", *maxRateOfReaction);
            }
            currenttimestep += 1.0;            
            PerformOneTimeStepRel(tskitstatus, isabsolute, isburninphaseover, ismodular, elementsperlb, &treesequencetablecollection, wholepopulationnodesarray, wholepopulationsitesarray, popsize, totaltimesteps, currenttimestep,wholepopulationwistree, wholepopulationwisarray, wholepopulationgenomes, psumofwis, chromosomesize, numberofchromosomes, totalindividualgenomelength, deleteriousmutationrate, beneficialmutationrate, Sb, beneficialdistribution, Sd, deleteriousdistribution, parent1gamete, parent2gamete, randomnumbergeneratorforgamma, miscfilepointer, *maxRateOfReaction, michaelisConstant, recessivityRunFlag, initializationValRel, addNonNeutral);  
        }
        if (i > 1000) {
            fflush(miscfilepointer);
        }


        //Following code calculates the variance in log(fitness) of the population after this generation of births and deaths.
        //May use an imprecise algorithm -- check before using as data.
        variancesum = CalculateVarianceInLogFitness(popsize, wholepopulationwisarray, *psumofwis);
        
        
        //This is the main data output, currently the summed fitness and variance in log(fitness) in the population.
        // Note that data is only written to the raw data file when i is a multiple of the WRITEFREQUENCY
        // global variable in order to limit runtimes. The final generation is also always written to ensure that
        // the final state is consistently known.

        // In order to check to see whether the population is actually undergoing deleterious dynamics
        // when we expect it to do so, an additional piece of information is now getting written to 
        // the output file, which is the standard deviation in log fitness, determined as the square
        // root of the variance in log fitness. The mean fitness is also calculated as the sum of 
        // wi's divided by the number of individuals. 
        if (i % WRITEFREQUENCY == 0 || i == Nxtimesteps - 1) {
            // double stdeviation = pow(variancesum, 0.5);
            // double meanFitness = *psumofwis / popsize;

            fprintf(rawdatafilepointer, "%d,%Lf\n", i+1, *psumofwis);
            fflush(rawdatafilepointer);
        }

        //fprintf(rawdatafilepointer, "%d \n", i+1);
        //fflush(rawdatafilepointer);

        //Following code periodically simplifies the tree sequence tables, which requires the tables to be sorted each time.
        //The tskit API calls out sorting each time as inefficient, but they haven't yet uploaded an example of how to do it differently. The API currently suggests that they would use the C++ API, which I don't want to deal with.
        //Simplify interval is currently set to 5. Should be either an input parameter or global variable at some point, probably.
        if (tskitstatus > 0){
            if (i % 10 == 0) {
                returnvaluefortskfunctions = tsk_table_collection_sort(&treesequencetablecollection, NULL, 0);
                check_tsk_error(returnvaluefortskfunctions);
            
                returnvaluefortskfunctions = tsk_table_collection_simplify(&treesequencetablecollection, wholepopulationnodesarray, (2*popsize), 0, NULL);
                check_tsk_error(returnvaluefortskfunctions);
            
                for (k = 0; k < (2*popsize); k++) {
                    wholepopulationnodesarray[k] = k;
                }   
            }
        }

        double c0, cov00, cov01, cov11, sumsq;

        if (isburninphaseover == 0) {
            UpdateLast200NTimeSteps(last200Ntimestepsvariance, variancesum);
            UpdateLast200NTimeSteps(literallyjustlast200Ntimesteps, i+1);
            if (i > 199) {           //to avoid calling the end of the burn-in phase at generation one
                                    //because of setting pre-simulation generations to zeroes
                                    //I just won't start looking for the end of the burn-in phase until 200 generations
                                    //This would be a mild problem if a simulation should end in 200 generations, but that shouldn't ever happen with the DFE I'm using.
                slopeofvariance = 0.0;
                gsl_fit_linear(literallyjustlast200Ntimesteps, step, last200Ntimestepsvariance, step, 200, &c0, &slopeofvariance, &cov00, &cov01, &cov11, &sumsq);
                if (slopeofvariance < arbitrarynumber) {
                    endofburninphase = i;
                    endofdelay = endofburninphase + 500;
                    isburninphaseover = 1;
                    fprintf(miscfilepointer, "Burn-in phase called as ending in generation %d\n", i+1);
                    fprintf(summarydatafilepointer, "Burn-in phase called as ending in generation %d\n", i+1);
                    if (VERBOSE == 1) {
                        fflush(miscfilepointer);
                        fflush(summarydatafilepointer);
                    }
                    
                }
            }
        }        
        
        //This is to produce a histogram of the wis of the entire population from a single generation.
        //It's terrible and completely non-modular, but I just can't bring myself to add in two more user-input arguments.
        if (typeofrun == 1) {
            if (i == 1999) {
                
                if (INDIVIDUALWIDATA == 1) {
                    if (VERYVERBOSE == 1) {
                        fprintf(veryverbosefilepointer, "Just before individual wi data lines.\n");
                        fflush(veryverbosefilepointer);
                    }                    
                    fprintf(summarydatafilepointer, "Individual, Wi\n");
                    for (k = 0; k < popsize; k++) {
                        fprintf(summarydatafilepointer, "%d,%Lf\n", k+1, wholepopulationwisarray[k]);
                    }
                }
            }
        }
        
        //If the burn-in phase has been called, wait 500 generations to start recording fitnesses.
        //This is to be sure that even when the beneficial rates/sizes are large, the only generations recorded will be from the uniformly sloping part of the simulation.
        //The average fitness from any generation after this delay period is recorded in the array of average fitnesses.
        if (i > endofdelay) {
            logaveragefitnesseachNtimesteps[Nxtimestepsafterburnin] = log((double) *psumofwis / (double) popsize);
            if (VERYVERBOSE == 1) {
                fprintf(veryverbosefilepointer, "log average fitness in generation %d, %d generations after burn-in, is: %f\n", i, Nxtimestepsafterburnin, logaveragefitnesseachNtimesteps[Nxtimestepsafterburnin]);
                fflush(veryverbosefilepointer);
            }
            Nxtimestepsafterburnin += 1;
        }
        
        //These lines ensure that the magnitude of fitness hasn't declined by too much.
        //At extremely small fitness values, floating-point math becomes imprecise.
        //These lines end the simulation if fitness declines below 10^-10, which should represent a completely degraded population.
        currentfittestindividualswi = FindFittestWi(wholepopulationwisarray, popsize);
        if (currentfittestindividualswi < pow(10.0, -10.0)) {
            fprintf(miscfilepointer, "\nFitness declined to less than 10^-10 during generation %d.", i+1);
            fprintf(summarydatafilepointer, "Fitness declined to catastrophic levels in generation %d.\n", i+1);
            endofsimulation = i;
            i = Nxtimesteps;
            didpopulationcrash = 1;
        }
    }
    
    //END OF SIMULATION FOR LOOP
    //I don't guarantee that the simplify interval matches up with the number of Nxtimesteps
    //Tree sequence recording requires that tables are sorted on the back end,
    //so I sort once again here at the end to ensure that all tables are sorted before they're read to file.
    //This might be inefficient, I'm not sure.
    if(tskitstatus > 0){
        returnvaluefortskfunctions = tsk_table_collection_sort(&treesequencetablecollection, NULL, 0);
        check_tsk_error(returnvaluefortskfunctions);

        returnvaluefortskfunctions = tsk_table_collection_simplify(&treesequencetablecollection, wholepopulationnodesarray, (2*popsize), 0, NULL);
        check_tsk_error(returnvaluefortskfunctions);
    
        for (k = 0; k < (2*popsize); k++) {
            wholepopulationnodesarray[k] = k;
        }
    
    //Printing out the node table in a way readable by python on the back end.
        fprintf(nodefilepointer, "is_sample time\n");
        for (k = 0; k < tablepointer->nodes.num_rows; k++) {
            if (k < (2*popsize)) {
                fprintf(nodefilepointer, "1 %f\n", tablepointer->nodes.time[k]);
            } else {
                fprintf(nodefilepointer, "0 %f\n", tablepointer->nodes.time[k]);
            }
        }
    
    //Printing out the edge table in a way readable by python on the back end.
        fprintf(edgefilepointer, "left right parent child\n");
        for (k = 0; k < tablepointer->edges.num_rows; k++) {
            fprintf(edgefilepointer, "%f %f %d %d\n", tablepointer->edges.left[k], tablepointer->edges.right[k], tablepointer->edges.parent[k], tablepointer->edges.child[k]);
        }
    
    //Printing out the site table in a way readable by python.
        fprintf(sitefilepointer, "position ancestral_state\n");
        for (k = 0; k < tablepointer->sites.num_rows; k++) {
            fprintf(sitefilepointer, "%f 0.0\n", tablepointer->sites.position[k]);
        }
    
    //Printing out the mutation table in a way readable by python.
        fprintf(mutationfilepointer, "site node derived_state\n");
        for (k = 0; k < tablepointer->mutations.num_rows; k++) {
            fprintf(mutationfilepointer, "%d %d %.12s\n", tablepointer->mutations.site[k], tablepointer->mutations.node[k], (tablepointer->mutations.derived_state + k*12));
        }
    }

    if (VERYVERBOSE == 1) {
        fprintf(veryverbosefilepointer, "Finished simulation with mean sb %.6f. Final population summed fitness was: %Lf\n", Sb, *psumofwis);
        fflush(veryverbosefilepointer);
    }

    if (didpopulationcrash == 0) {
        endofsimulation = i;
    }

    if (isburninphaseover == 1) {
        if (VERYVERBOSE == 1) {
            fprintf(veryverbosefilepointer, "Calculating slope of log fitness with the following parameters: endofsimulation = %d, endofdelay = %d, generationsafterburnin = %d\nLog fitness each generation: ", endofsimulation, endofdelay, Nxtimestepsafterburnin);
            for (j = 0; j < (endofsimulation - endofdelay); j++) {
                fprintf(veryverbosefilepointer, "%f ", logaveragefitnesseachNtimesteps[j]);
            }
            fprintf(veryverbosefilepointer, "\n");
            fflush(veryverbosefilepointer);
        }
        slopeoflogfitness = CalculateSlopeOfLogFitness(endofsimulation, endofdelay, logaveragefitnesseachNtimesteps);
        fprintf(summarydatafilepointer, "Slope of log(fitness) after the burn-in phase: %f\n", slopeoflogfitness);
        
        if (VERBOSE == 1) {
            fflush(summarydatafilepointer);
            fflush(rawdatafilepointer);
        }
        
        fclose(rawdatafilepointer); 
        fclose(summarydatafilepointer);
        fclose(nodefilepointer);
        fclose(edgefilepointer);
        fclose(sitefilepointer);
        fclose(mutationfilepointer);

        free(rawdatafilename);
        free(summarydatafilename);
        free(logaveragefitnesseachNtimesteps);
        free(literallyjustlast200Ntimesteps);
        free(last200Ntimestepsvariance);
        free(wholepopulationgenomes);
        free(wholepopulationwistree);
        free(wholepopulationwisarray);
        free(wholepopulationnodesarray);
        free(sortedwisarray);
        
        tsk_table_collection_free(&treesequencetablecollection);
        
        return slopeoflogfitness;
    }

    if (isburninphaseover == 0) {
        fprintf(summarydatafilepointer, "End of burn-in phase not reached.");
        
        fclose(rawdatafilepointer); 
        fclose(summarydatafilepointer);
        fclose(nodefilepointer);
        fclose(edgefilepointer);
        fclose(sitefilepointer);
        fclose(mutationfilepointer);

        free(rawdatafilename);
        free(summarydatafilename);
        free(logaveragefitnesseachNtimesteps);
        free(literallyjustlast200Ntimesteps);
        free(last200Ntimestepsvariance);
        free(wholepopulationgenomes);
        free(wholepopulationwistree);
        free(wholepopulationwisarray);
        free(wholepopulationnodesarray);
        free(sortedwisarray);
        
        tsk_table_collection_free(&treesequencetablecollection);
        
        return -1.0;
    }
}


void PerformOneTimeStepRel(int tskitstatus, bool isabsolute, int isburninphaseover, bool ismodular, int elementsperlb, tsk_table_collection_t *treesequencetablecollection, tsk_id_t * wholepopulationnodesarray, tsk_id_t * wholepopulationsitesarray, int popsize, int totaltimesteps, double currenttimestep, long double *wholepopulationwistree, long double *wholepopulationwisarray, double *wholepopulationgenomes, long double * psumofwis, int chromosomesize, int numberofchromosomes, int totalindividualgenomelength, double deleteriousmutationrate, double beneficialmutationrate, double Sb, int beneficialdistribution, double Sd, int deleteriousdistribution, double *parent1gamete, double *parent2gamete, gsl_rng * randomnumbergeneratorforgamma, FILE *miscfilepointer, long double maxRateOfReaction, double michaelisConstant, int recessivityRunFlag, double initializationValRel, int addNonNeutral)
{
    int currentparent1, currentparent2, currentvictim;

    currentvictim = ChooseVictim(popsize);
    currentparent1 = ChooseParentWithTree(wholepopulationwistree, popsize, *psumofwis, miscfilepointer);
    currentparent2 = ChooseParentWithTree(wholepopulationwistree, popsize, *psumofwis, miscfilepointer);
    while (currentparent1 == currentparent2) { //probably not ideal, since it'll never break with population sizes of zero or one.
        currentparent2 = ChooseParentWithTree(wholepopulationwistree, popsize, *psumofwis, miscfilepointer);
    }
    
    tsk_id_t childnode1, childnode2;
   
    RecombineChromosomesIntoGamete(isabsolute, tskitstatus, ismodular, elementsperlb, isburninphaseover, treesequencetablecollection, wholepopulationnodesarray, &childnode1, totaltimesteps, currenttimestep, currentparent1, chromosomesize, numberofchromosomes, parent1gamete, wholepopulationgenomes, totalindividualgenomelength);
    ProduceMutatedGamete(tskitstatus, isburninphaseover, treesequencetablecollection, wholepopulationnodesarray, wholepopulationsitesarray, &childnode1, totaltimesteps, currenttimestep, currentparent1, isabsolute, totalindividualgenomelength, deleteriousmutationrate, beneficialmutationrate, Sb, beneficialdistribution, Sd, deleteriousdistribution, parent1gamete, randomnumbergeneratorforgamma, miscfilepointer, addNonNeutral);
        
    RecombineChromosomesIntoGamete(isabsolute, tskitstatus, ismodular, elementsperlb, isburninphaseover, treesequencetablecollection, wholepopulationnodesarray, &childnode2, totaltimesteps, currenttimestep, currentparent2, chromosomesize, numberofchromosomes, parent2gamete, wholepopulationgenomes, totalindividualgenomelength);
    ProduceMutatedGamete(tskitstatus, isburninphaseover, treesequencetablecollection, wholepopulationnodesarray, wholepopulationsitesarray, &childnode2, totaltimesteps, currenttimestep, currentparent2, isabsolute, totalindividualgenomelength, deleteriousmutationrate, beneficialmutationrate, Sb, beneficialdistribution, Sd, deleteriousdistribution, parent2gamete, randomnumbergeneratorforgamma, miscfilepointer, addNonNeutral);
               
    //next variables are only used for absolute simulations, however since InitializePopulation is a shared function I must initialize them here, even if I don't use them again
    long double *pInverseSumOfWis;
    long double *pInverseSumOfWissquared;
    long double *psumofload;
    long double *psumofloadsquared;
    bool *wholepopulationisfree;
    int *wholepopulationindex;
    long double *wholepopulationdeathratesarray;
    int *pPopSize;
	double b_0, r, s;
	int i_init;
    
    PerformDeath(isabsolute, tskitstatus, isburninphaseover, popsize, pPopSize, currentvictim, deleteriousdistribution, wholepopulationwistree, wholepopulationwisarray, wholepopulationdeathratesarray, wholepopulationindex, wholepopulationisfree, psumofwis, pInverseSumOfWis, pInverseSumOfWissquared, b_0, r, i_init, s, psumofload, psumofloadsquared, wholepopulationnodesarray, miscfilepointer);
    
    PerformBirth(tskitstatus, isburninphaseover, ismodular, elementsperlb, treesequencetablecollection, wholepopulationnodesarray, childnode1, childnode2, isabsolute, parent1gamete, parent2gamete, popsize, pPopSize, currentvictim, wholepopulationgenomes, totalindividualgenomelength, deleteriousdistribution, wholepopulationwistree, wholepopulationwisarray, wholepopulationdeathratesarray,wholepopulationindex, wholepopulationisfree, psumofwis, pInverseSumOfWis, pInverseSumOfWissquared, b_0, r, i_init, s, psumofload, psumofloadsquared, miscfilepointer, maxRateOfReaction, michaelisConstant, recessivityRunFlag, initializationValRel);
    
}

void InitializePopulationRel(int tskitstatus, tsk_table_collection_t * treesequencetablecollection, tsk_id_t * wholepopulationnodesarray, tsk_id_t * wholepopulationsitesarray, long double *wholepopulationwistree, long double *wholepopulationwisarray, int popsize, double *wholepopulationgenomes, int totalpopulationgenomelength, int totaltimesteps, long double * psumofwis, double initializationValRel, long double maxRateOfReaction, double michaelisConstant, int chromosomeSize, int recessivityRunFlag, int numberOfChromosomes) 
{
    int i, j;
    
    double haploidgenomelength = (double) ((totalpopulationgenomelength / popsize) / 2);
    
    // To limit computational issues for small a0's, we multiply by 1/W(no mutations) each time
    // fitness is calculated to effectively "normalize" values while maintaining relative relationships
    // In this case, to avoid unnecessary computation, all wis are simply initialized to 1.0
    for (i = 0; i < popsize; i++){
            // Note that, due to normalization, whether or not we have recessivity turned on, all
            // individuals are initialized as mutationless with a wi of 1.0
            wholepopulationwistree[i] = 1.0;
            wholepopulationwisarray[i] = 1.0;
    }
    //this for loop taken from the Fen_init function in sample implementation from 'Fenwick tree' Wikipedia page.
    for (i = 0; i < popsize; i++) {
        j = i + LSB(i+1);
        if (j < popsize) {
            wholepopulationwistree[j] += wholepopulationwistree[i];
        }
    }
    
    // Again, due to normalization, whether or not recessivity is turned on, all individuals
    // are initialized with a wi of 1.0
    *psumofwis = (long double)popsize;
    
    for (i = 0; i < totalpopulationgenomelength; i++){
        if (recessivityRunFlag == 1) {
            // Case for recessivity- initialize with a particular value of a0
            wholepopulationgenomes[i] = log(initializationValRel);
        } else {
            // Case for no recessivity- gene activity initialized to 0 (log(1))
            wholepopulationgenomes[i] = 0.0;
        }
    }
    //The following lines initialize the node table for tree sequence recording.
    //Note that nodes here are single sets of chromosomes, so the 2x popsize here assumes diploidy.
    if (tskitstatus > 0){
	treesequencetablecollection->sequence_length = haploidgenomelength;

        for (i = 0; i < (2 * popsize); i++) {
            wholepopulationnodesarray[i] = tsk_node_table_add_row(&treesequencetablecollection->nodes, 0, totaltimesteps, TSK_NULL, TSK_NULL, NULL, 0);
            check_tsk_error(wholepopulationnodesarray[i]);
        }
    
    // Nothing changes in the nodes, but the 0 represents the initialization- will need to change this!
    //The following lines add a site to the tree sequence recording site table corresponding to each linkage block, with ancestral state of 0.
        for (i = 0; i < haploidgenomelength; i++) {
            wholepopulationsitesarray[i] = tsk_site_table_add_row(&treesequencetablecollection->sites, i, "0.00000000", 10, NULL, 0);
            check_tsk_error(wholepopulationsitesarray[i]);
        }
    }
}

//In this model, individuals die at random. There's no selection happening here.
int ChooseVictim(int populationsize)
{
    int randomindividual = pcg32_boundedrand(populationsize);
    return randomindividual;
}

//The tree in the name is the Fenwick tree, which stores the fitnesses of individuals in the population.
//This function is where selection occurs -- individuals with higher-than-average fitness will be chosen more often as parents.
int ChooseParentWithTree(long double *wholepopulationwistree, int popsize, long double sumofwis, FILE *miscfilepointer)
{
    long double randomnumberofbirth;
    int newparent = 0;
    
    randomnumberofbirth = (ldexp(pcg32_random(), -32)) * sumofwis;
    //Above line generates a random integer between 0 and 2^32, then multiplies by 2^-32
    //to generate a float between 0 and 1 and then multiplies by the sum of wis
    //to get a number between 0 and the sum of wis.
    
    int leftbound, rightbound;
    leftbound = 0;
    rightbound = popsize;
    if (leftbound > rightbound) {
        fprintf(miscfilepointer, "\nError: population size is %d.", popsize);
        return -1;
    }
    //Above lines initialize the variables necessary for the SearchTree function and check for an extinct population.
    
    newparent = (SearchTree(leftbound, rightbound, randomnumberofbirth, wholepopulationwistree));
    return newparent;
}

// Method for calculating Wi without recessivity
double CalculateWiNoRecessivity(double *parent1gamete, double *parent2gamete, int totalindividualgenomelength)
{
    double newwi = 0.0;
    long double currentlinkageblockssum = 0.0;
    int i;

    for (i = 0; i < (totalindividualgenomelength/2); i++) {
        currentlinkageblockssum += parent1gamete[i];
        currentlinkageblockssum += parent2gamete[i];
    }
    newwi = exp(currentlinkageblockssum);
    return newwi;
}

// Method for calculating Wi with recessivity
// Note that this now involves multiplying each wi by 1/W(no mutations) in order to normalize our
// wi values and limit the computational errors resulting from these values being too small in magnitude
// NOTE: IT IS ASSUMED THAT THE VALUE PASSED IN FOR VMAX IS MEANT TO BE OUTSIDE OF THE MULTIPLICATION
double CalculateWiWithRecessivity(double *parent1gamete, double *parent2gamete, int totalindividualgenomelength, long double maxRateOfReaction, double michaelisConstant, double initializationValRel) {
    // Initializes the new fitness contribution value
    double newWiFitnessContribution = 0.0;
    long double currentFitnessContributionLogSum = 0.0;
    double FindFitnessContribution(double linkageBlockActivity, long double maxRateOfReaction, double michaelisConstant);

    double nonExp = (2.0 * initializationValRel * maxRateOfReaction) / (michaelisConstant + (2.0 * initializationValRel));

    // Loops through the given linkage block array
    int loopIndex;
    for (loopIndex = 0; loopIndex < (totalindividualgenomelength/2); loopIndex++) {
        // Gets the current log values for linkage block activity
        double logLinkageBlockActivityOne = parent1gamete[loopIndex];
        double logLinkageBlockActivityTwo = parent2gamete[loopIndex];

        // Exponentiates these values to get a
        double linkageBlockActivityOne = exp(logLinkageBlockActivityOne);
        double linkageBlockActivityTwo = exp(logLinkageBlockActivityTwo);

        // Calculates the fitness contribution for this pair of linkage blocks
        // Note that Vmax is currently being pulled out of this calculation, and instead is 
        // multiplied onto the final result for fitness
        double bothLinkageBlockFitnessContribution = FindFitnessContribution(linkageBlockActivityOne + linkageBlockActivityTwo, maxRateOfReaction, michaelisConstant);

        // Adds the log of this value to the current fitness contribution log sum
        currentFitnessContributionLogSum += log(bothLinkageBlockFitnessContribution);
    }

    // Exponentiates the sum of all of the log fitness contributions, then multiplies
    // by Vmax, to get the new Wi
    newWiFitnessContribution = exp(currentFitnessContributionLogSum);
    newWiFitnessContribution = newWiFitnessContribution * maxRateOfReaction;

    // Multiplies this by 1/W(no mutations) to normalize wi
    // Note that the exponent represents chromosomesize * numberofchromosomes
    newWiFitnessContribution = pow(nonExp, -(totalindividualgenomelength / 2)) * newWiFitnessContribution;

    // Returns this value
    return(newWiFitnessContribution);

}

// Helper method for calculating the expected fitness contribution for a pair of linkage blocks
double FindFitnessContribution(double bothLinkageBlockActivity, long double maxRateOfReaction, double michaelisConstant) {
    // Returns linkage block activity times Vmax, divided by (Km + linkage block activity)
    // Currently, Vmax is not being used in this calculation- instead, it is pulled out of the multiplication
    // all together to simplify the fitness renormalization process
    return((bothLinkageBlockActivity)/(michaelisConstant + bothLinkageBlockActivity));
}

// Wrapper method for calculating Wi. Based on the value of the recessivity flag (with 0 indicating
// a run not including recessivity and 1 indicating a run including recessivity), calculates
// and returns the new value of Wi, either with or without recessivity
double CalculateWi(double *parent1gamete, double *parent2gamete, int totalindividualgenomelength, long double maxRateOfReaction, double michaelisConstant, int recessivityRunFlag, double initializationValRel) {
    double CalculateWiWithRecessivity(double *parent1gamete, double *parent2gamete, int totalindividualgenomelength, long double maxRateOfReaction, double michaelisConstant, double initializationValRel);
    double CalculateWiNoRecessivity(double *parent1gamete, double *parent2gamete, int totalindividualgenomelength);

    if (recessivityRunFlag == 0) {
        return(CalculateWiNoRecessivity(parent1gamete, parent2gamete, totalindividualgenomelength));
    } else {
        return(CalculateWiWithRecessivity(parent1gamete, parent2gamete, totalindividualgenomelength, maxRateOfReaction, michaelisConstant, initializationValRel));
    }
}

// Helper method which can be used to re-normalize fitness. Takes in pointers to everything that needs to be 
// updated with the re-normalized fitnesses, and uses the current mean gene activity to calculate the factor
// used for the re-normalization.
// NOTE: I can't find any uses for sortedwisarray, so I'm not currently updating this- may cause issues, but I can't find 
// any points in the code where it's actually being used. 
void renormalizeFitness(double *wholepopulationgenomes, long double *wholepopulationwisarray, long double *sortedwisarray, long double *psumofwis, long double *wholepopulationwistree, int totalpopulationgenomelength, long double *maxRateOfReaction, double michaelisConstant, int popsize, FILE *miscfilepointer, int totalindividualgenomelength) {
    // First, obtains the current mean fitness
    long double mean_fitness = *psumofwis / (long double) popsize;

    // Divides all of the values in the array of wi's by this value, reconstructing the sum of wi's
    // and creating an empty Fenwick tree as we do so
    *psumofwis = 0;

    // long double actual_range1 = 0;
    // long double actual_range2 = 0;
    for (int i = 0; i < popsize; i++) {
        wholepopulationwisarray[i] = wholepopulationwisarray[i] / mean_fitness;
        // wholepopulationwistree[i] = 0;
        *psumofwis += wholepopulationwisarray[i];

        /*
        if (i >= 676 && i < 823) {
            actual_range1 += wholepopulationwisarray[i];
        }

        if (i >=2 && i < 587) {
            actual_range2 += wholepopulationwisarray[i];
        }
        */
    }

    /*
    // Initializes a new Fenwick tree with the updated values
    for (int i = 0; i < popsize; i++) {
        Fen_set(wholepopulationwistree, popsize, wholepopulationwisarray[i], i);
    }
    

    // CURRENTLY, FOR DEBUGGING PURPOSES, SUMS FENWICK TO SEE IF SUM IS EQUAL TO SUM OF
    // FITNESSES
    long double tree_sum = Fen_sum(wholepopulationwistree, popsize);
    fprintf(miscfilepointer, "Current fitness sum: %Lf; Current Fenwick sum: %Lf\n", *psumofwis, tree_sum);

    // ALSO CHECKS TO MAKE SURE THAT A FEW RANGES MATCH
    long double tree_range1 = Fen_range(wholepopulationwistree, 676, 823);
    fprintf(miscfilepointer, "Current range1 (676 to 823) sum: %Lf; Current Fenwick range1 sum: %Lf\n", actual_range1, tree_range1);

    long double tree_range2 = Fen_range(wholepopulationwistree, 2, 587);
    fprintf(miscfilepointer, "Current range2 (2 to 587) sum: %Lf; Current Fenwick range2 sum: %Lf\n", actual_range2, tree_range2);
    */

    // Updates Vmax, simply dividing it by the calculated mean fitness, since it is assumed that 
    // Vmax is outside of the multiplication loop for calculating fitness
    *maxRateOfReaction = *maxRateOfReaction / mean_fitness;
}


// Helper method which can be used to capture the average activity and the average fitness contribution
// for all genes (linkage blocks) across all individuals. The method takes in a string representing the
// name of the file to write to, and writes the specified information to a file with this name. Note that
// this method is specifically used to check that the fitness functionalities of the program are working
// as expected, and to provide an alternate method for debugging.
void captureGeneActivitiesAndFitnessContributions(double *wholepopulationgenomes, char *captureFilename, int totalindividualgenomelength, int totalpopulationgenomelength, int popsize, long double maxRateOfReaction, double michaelisConstant) {
    // Opens the given file in write mode, giving an error if unable to do so
    FILE *fp = fopen(captureFilename, "w");
    if (fp == NULL) {
        perror(captureFilename);
        exit(1);
    }

    // Adds the column labels to the csv file that will be output
    fprintf(fp, "Gene,Average.activity,Average.fitness.contribution\n");

    // Finds the haploid genome length
    int haploidGenomeLength = totalindividualgenomelength / 2;

    // Initializes an array which will store the sum of the activities associated with each gene, which 
    // will be used for finding the average activity associated with each gene. This array will be filled
    // with 0's initially for simplicity
    double *geneActivitySums = malloc(sizeof(double) * haploidGenomeLength);
    for (int i = 0; i < haploidGenomeLength; i++) {
        geneActivitySums[i] = 0;
    }

    // Iterates through the full wholepopulationgenomes array
    int i = 0;
    int nextHaploidIndex, curGeneIndex;
    double geneCopyActivitySum;
    while (i < totalpopulationgenomelength) {
        // Gets the index of the start of the second set of chromosomes for the current individual
        // within the array
        nextHaploidIndex = i + haploidGenomeLength;

        // Iterates through the current individual
        curGeneIndex = 0;
        while (curGeneIndex < haploidGenomeLength) {
            // Exponentiates the gene copy activities, then sums them
            geneCopyActivitySum = exp(wholepopulationgenomes[i]) + exp(wholepopulationgenomes[nextHaploidIndex]);

            // Adds this sum of gene copy activities to the appropriate spot in our array of gene activity
            // sums
            geneActivitySums[curGeneIndex] += geneCopyActivitySum;

            // Increments the current gene index
            curGeneIndex++;

            // Increments the indices for the wholepopulationgenomes array
            i++, nextHaploidIndex++;
        }

        // Once all of the genes for the current individual have been considered, we then move onto the
        // next individual, setting i to be nextHaploidIndex, since this current represents the index for
        // the start of the next individual
        i = nextHaploidIndex;
    }

    // Now that we have our array of gene sums, we divide each entry in the array to get the average activity
    // for that gene across all individuals, and plug this value into our Michaelis-Menten function in order to
    // obtain the average fitness contribution for that gene across all individuals
    double averageGeneActivity, averageFitnessContribution;
    for (i = 0; i < haploidGenomeLength; i++) {
        // Calculates these values
        averageGeneActivity = geneActivitySums[i] / popsize;
        averageFitnessContribution = FindFitnessContribution(averageGeneActivity, maxRateOfReaction, michaelisConstant);

        // Adds them to our file
        fprintf(fp, "%d,%.18f,%.18f\n", i + 1, averageGeneActivity, averageFitnessContribution);
    }

    // Frees the space associated with our gene activity sums array
    free(geneActivitySums);
    
    // Closes the file, giving an error if we'e unable to do so
    int closeResult = fclose(fp);
    if (closeResult == EOF) {
        fprintf(stderr, "Error: unable to close capture file\n");
        exit(1);
    }
}


// Helper method for capturing the activity and fitness contribution of a particular gene across
// all individuals, which can help us to get an idea of how different versions of this gene are 
// "competing" with one another. Note that this method automatically captures this information for
// the first gene in each individual's genome.
void captureSingleGeneActivityAndFitnessContribution(double *wholepopulationgenomes, char *captureFilename, int totalindividualgenomelength, int totalpopulationgenomelength, int popsize, long double maxRateOfReaction, double michaelisConstant) {
    // First, we attempt to open the given file in write mode, giving an error and exiting if unable to 
    // do so
    FILE *fp = fopen(captureFilename, "w");
    if (fp == NULL) {
        perror(captureFilename);
        exit(1);
    }

    // Finds the haploid genome length
    int haploidGenomeLength = totalindividualgenomelength / 2;

    // Adds appropriate headers to the output file
    fprintf(fp, "Individual,GeneOneActivity,GeneOneFitnessContribution\n");

    // Iterates through the full wholepopulationgenomes array, capturing information about the gene of interest
    int i = 0;
    int curIndividualIndex = 0;
    double curActivity, curFitnessContribution;
    while (i < totalpopulationgenomelength) {
        // Gets the activity for the current gene, writing this value to the output
        double curActivity = exp(wholepopulationgenomes[i]) + exp(wholepopulationgenomes[i + haploidGenomeLength]);
        fprintf(fp, "%d,%.18f,", curIndividualIndex + 1, curActivity);

        // Finds the fitness contribution associated with this gene, writing this to the output as well
        curFitnessContribution = FindFitnessContribution(curActivity, maxRateOfReaction, michaelisConstant);
        fprintf(fp, "%.18f\n", curFitnessContribution);

        // Updates our indices
        i += totalindividualgenomelength;
        curIndividualIndex++;
    }

    // Finally, the file we've been writing to is closed
    int closeResult = fclose(fp);
    if (closeResult == EOF) {
        fprintf(stderr, "Error: file could not be closed\n");
        exit(1);
    }
}
