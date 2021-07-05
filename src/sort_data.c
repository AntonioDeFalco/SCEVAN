#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <R.h>

/* Maximum input line for the data file (1024*100 allows to analyze about 2000
 samples) You can increase this constant for more large datasets */
#define MAX_INPUT_LINE 102400

// Number of 
#define ROUND_CNST 1000

/*
 * Input parameters 
 */
char *dataset;
char *output_file_name;

/*
 * Probes have to be sorted 
 */
struct probe **probes;
char *header;

/*
 * Dataset properties 
 */
int num_cols = 0;
int num_rows = 0;

/*
 * Support Types 
 */
typedef float *data_ptr;
typedef struct probe {
    char *snp;
    int chromosome;
    int position;
    struct probe *next;
    data_ptr data;
} probe;


// Function to sort data form file
void sort_dataset();

// Function to write data in output file
void write_data();

// Function to compare probes
int cmp_probes(const void *p1, const void *p2)
{
    struct probe *const *a = p1;
    struct probe *const *b = p2;

    if ((*a)->chromosome < (*b)->chromosome) {
        return -1;
    } else {
        if ((*a)->chromosome > (*b)->chromosome) {
            return 1;
        } else {
            if ((*a)->position < (*b)->position) {
                return -1;
            } else {
                if ((*a)->position > (*b)->position) {
                    return 1;
                } else {
                    return 0;
                }
            }
        }
    }
}

/*--------End Funtion Headers -------------------*/

int sort_data(char **data, char **output)
{

    Rprintf("\n------ START DATA SORTING  ------\n");
    /*
     * Takes dataset and output files 
     */
    dataset = *data;
    output_file_name = *output;

    /*
     * Load Data 
     */
    Rprintf("\nSTEP 1: Sorting Data\n");
    sort_dataset();
    Rprintf("\t- Number of Observations = %d\n\t- Number of Probes = %d\n",
            num_cols - 3, num_rows);
    Rprintf("\tData Correctly Sorted\n");

    /*
     * Write Data 
     */
    Rprintf("\nSTEP 2: Writing Results in '%s' Output File\n",
            output_file_name);
    write_data();
    Rprintf("\tFile '%s' Correctly Written\n", output_file_name);
    Rprintf("\n------ DATA SORTING SUCCESFULLY COMPLETED------\n");
    return 0;
}

void sort_dataset()
{

    char line[MAX_INPUT_LINE];
    char *elem;
    char *sep = "\t";
    probe *temp_probe;
    int nchar;
    num_cols = 0;
    num_rows = 0;
    FILE *file = fopen(dataset, "r");
    /*
       Read header of Data File 
     */
    if (fgets(line, sizeof (line), file) != 0) {
        header = malloc(sizeof (char) * (strlen(line) + 1));
        strcpy(header, line);

        elem = strtok(line, sep);

        while (elem) {
            num_cols++;
            elem = strtok('\0', sep);
        }
    }
    /*
     * Read Data 
     */
    int index;
    while (fgets(line, sizeof (line), file) != 0) {

        probes =
            (struct probe **) realloc(probes,
                                      (num_rows +
                                       1) * sizeof (struct probe *));

        temp_probe = (struct probe *) malloc(sizeof (struct probe));
        temp_probe->data = (data_ptr) malloc(sizeof (float) * (num_cols - 3));
        elem = strtok(line, sep);
        index = 0;
        while (elem) {
            if (index == 0) {
                nchar = strlen(elem);
                temp_probe->snp = malloc(sizeof (char) * (nchar + 1));
                strcpy(temp_probe->snp, elem);
            }

            if (index == 1) {
                if (strcmp(elem, "X") == 0) {
                    temp_probe->chromosome = 23;
                } else {
                    if (strcmp(elem, "Y") == 0) {
                        temp_probe->chromosome = 24;
                    } else {
                        if (strcmp(elem, "MT") == 0) {
                            temp_probe->chromosome = 25;
                        } else {
                            temp_probe->chromosome = atoi(elem);
                        }
                    }
                }
            }

            if (index == 2) {
                temp_probe->position = atoi(elem);
            }

            if (index > 2) {
                temp_probe->data[index - 3] = atof(elem);
            }
            index++;
            elem = strtok('\0', sep);
        }
        probes[num_rows] = temp_probe;
        num_rows++;
    }
    fclose(file);

    qsort(probes, num_rows, sizeof (struct probe *), cmp_probes);
}

void write_data()
{
    FILE *file = fopen(output_file_name, "w");
    int i, j;
    fprintf(file, "%s", header);

    for (i = 0; i < num_rows; i++) {
        fprintf(file, "%s\t%d\t%d\t", probes[i]->snp, probes[i]->chromosome,
                probes[i]->position);
        for (j = 0; j < (num_cols - 3); j++) {
            if (j < (num_cols - 4)) {
                fprintf(file, "%f\t", probes[i]->data[j]);
            } else {
                fprintf(file, "%f\n", probes[i]->data[j]);
            }
        }

    }
    fclose(file);
}
