#ifndef GLOBAL_VARS_H_INCLUDED
#define GLOBAL_VARS_H_INCLUDED 1

#define VERBOSE 0
#define VERYVERBOSE 0
#define MISCELLANEOUS 1
#define INDIVIDUALWIDATA 1
#define OVERLAPPINGGENERATIONS 1 //probably should be an input argument at some point.
#define LSB(i) ((i) & -(i)) //isolates least significant single bit for fenwick tree
#define PI 3.141592654

// Defines the frequency with which data is written to the raw data output file, in terms
// of generations (i.e. 5 corresponds to recording data every 5 generations)
#define WRITEFREQUENCY 5

// Defines the frequency with which the population should be re-normalized. May need to 
// be adjusted based on particular parameter combination to prevent things from crashing
#define RENORM_FREQUENCY 1000

#define check_tsk_error(val)                                                            \
    if (val < 0) {                                                                      \
        errx(EXIT_FAILURE, "line %d: %s", __LINE__, tsk_strerror(val));                 \
    }//error checking for tree sequence recording

#endif // GLOBAL_VARS_H_INCLUDED


