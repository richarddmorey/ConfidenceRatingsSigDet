# ConfidenceRatingsSigDet
C code for analyzing confidence ratings data using signal detection


This repository contains the code needed to analyse confidence-ratings recognition memory data ([Morey, Pratte, & Rouder, 2008](http://www.sciencedirect.com/science/article/pii/S0022249608000242)). 

The files included are:

File name| Description
---------|---------
sd2.c    |      C source code for analysis
sd2.exe  |     Windows executable for analysis
cygwin1.dll |   Libraries required for the executable to run in Windows
an2.dat     |  Example data (simulated data)
an2.chn     |  Example chains; the result of the analysis
analysis.R  |  R code to analyze the chains that the executable outputs
chains.pdf  |  PDF file of the chains from the example data, from the R analysis
items.pdf   |  PDF file of item effects from the example data, from the R analysis
----------------


### Data input format
The data files to be analyzed must be in a specific format to work. See the an2.dat file for an example.

Each ROW of the input file is one participant.
* ODD columns are the response (1...K) that the participant made to a particular item. 
* EVEN columns indicate whether the item in the column immediately to the left was studied or unstudied (1 or 0, respectively) 

Each column represents the same item for each participant. For instance, if the column 1 for participant 1 is their response for "COMPUTER," then column 1 for participant 2 (row 2) contains participant 2's response to "COMPUTER".

**Missing values: if there is no data for a participant-by-item combination, make the response for that combination -1 in the data file. Studied or unstudied can be anything.**

### Analysis output format

In the output chains, item 1 correponds to column 1, item 2 to column 2, etc. The same holds for rows/participants.

### Running an analysis in Windows
`sd.exe` runs the analysis. The executable must be used with arguments, like this:

    sd2.exe 10000 4 10000 2 1 2 1 .2 an2

Arguments:

1. MCMC iterations
2. Number of response categories (>2)
3. Prior variance on `mu_s` and `mu_n` population means
4. `a` for the IG prior on process variances `sigma^2_n` and `sigma^2_s`
5. `b` for the IG prior on process variances `sigma^2_n` and `sigma^2_s`
6. `e` for the IG prior on random effect variances
7. `f` for the IG prior on random effect variances
8. Metropolis step standard deviation
9. Name of analysis. 
10. [optional] Data file to analyze. If this is not specified, argument 9 determines the name (e.g., `[NAME].dat`)


Any questions can be directed to Richard Morey at richarddmorey@gmail.com.
