#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <R.h>

/* Maximum input line for the data file (1024*200 allows to analyze about 10000
 samples) You can increase this constant for more large datasets */
#define MAX_INPUT_LINE 512000


#define ROUND_CNST 1000

/* Input parameters */
char *dataset;
char *output_file_name;
float beta;
int min_region_probe_size;
int min_region_bp_size;
float loss_threshold;
float gain_threshold;
int baf[1];
float loh_threshold;
float loh_frequency;
int bs;
float pval_threshold;

/* Dataset properties */
int num_samples = 0;
int num_probes = 0;
int num_chromosomes;
int num_seg_regions = 0;

/* Support Types */
typedef float *data_ptr;
typedef struct probe {
    int chromosome;
    int position;
    struct probe *next;
    data_ptr data;
} probe;

/* In each seg_element the segmentation for a chromosome is saved */
typedef int *data_ptr_int;
typedef struct seg_element {
    int chromosome;
    int num_regions;
    struct seg_element *next;
    data_ptr_int start;
    data_ptr_int end;
    data_ptr_int size;
    data_ptr l2_mean;

} seg_element;

/* Dataset Matrix and Position Array */
float **lrr_matrix;
float **baf_matrix;
int *position_matrix;
int *chromosome_matrix;
/* Positions in which chromosomes start and end */
int *chr_brks_start, *chr_brks_end;

/* Arrays to Store Segmentation Properties */
int *seg_chromosomes;
int *seg_start;
int *seg_end;
int *seg_size;
float *seg_l2_mean;
float *seg_loss_pval;
int *seg_loss_perc;
float *seg_gain_pval;
int *seg_gain_perc;
float *seg_loh_pval;
int *seg_loh_perc;
int *num_loss_sample;
int *num_gain_sample;
int *num_loh_sample;

float *mean_loss;
float *mean_gain;
float *mean_loh;

/*----------Start Function Headers ----------------*/

int read_params(char *filename);

void load_data();
void call_VegaMC();
void compute_matrices();
void compute_pvalue();
float calc_mean(float vals[], int n_elem);
float calc_std(float vals[], int n_elem);
int generate_binomial(float theta);
void write_segementation();
/*--------End Funtion Headers -------------------*/

void vegaMC(float **data, int *markers_start, int *start, int *end, int *size,
            float *mean, int *n, float *be, float *std, int *n_reg,
            int *n_samples, float *weight, float *weight_sum);

int run_vegaMC(char **data, char **out, double *b, int *mrbs, double *losst,
               double *gaint, int *ba, double *loht, double *lohf, int *bsp,
               int *ns, int *np, int *nc)
{
    //Rprintf("\n------ START COMPUTATION ------\n");

    /* Read Input Parameters */
    //Rprintf("STEP 1: Loading Parameters\n");
    beta = (float) *b;
    min_region_bp_size = *mrbs;
    loss_threshold = (float) *losst;
    gain_threshold = (float) *gaint;
    baf[0] = *ba;
    loh_threshold = (float) *loht;
    loh_frequency = (float) *lohf;
    bs = *bsp;

    //Rprintf("\tParameters Correctly Loaded\n");

    /* Takes dataset and output files */
    dataset = *data;
    output_file_name = *out;

    /* Load Data */
    //Rprintf("\nSTEP 2: Loading Data\n");
    load_data();
    //Rprintf
    //    ("\t- Number of Samples = %d\n\t- Number of Probes = %d\n\t- Number of Chromosomes = %d\n",
    //     num_samples, num_probes, num_chromosomes);
    //Rprintf("\tData Correctly Loaded\n");
    beta = beta * (float) num_samples;
  
  //Rprintf("\nSTEP 3: Performing Joint Segmentation\n");
    call_VegaMC();
    //Rprintf
    //    ("\tJoint Segmentation Successfully Completed with %d Segmented Regions\n",
    //     num_seg_regions);

    /* Run VegaMC on Each Chromosome */
    //Rprintf("\nSTEP 4: Computing Aberration Matrices\n");
    compute_matrices();
    //Rprintf("\tAberration Matrices Successfully Computed\n", num_seg_regions);

    /* Assess the p-value for each region */
    //Rprintf("\nSTEP 5: Assessing Statistical Significativity\n");
    compute_pvalue();
    //Rprintf("\tStatistical Analysis Successfully Completed\n",
    //        num_seg_regions);

    //Rprintf("\nSTEP 6: Writing Results in '%s' Output File\n",
    //        output_file_name);
    write_segementation();
    //Rprintf("\tFile '%s' Correctly Written\n", output_file_name);

    *ns = num_samples;
    *np = num_probes;
    *nc = num_chromosomes;
    //Rprintf("\n------ COMPUTATION SUCCESFULLY COMPLETED------\n");
    return 0;
}


void write_segementation()
{
    FILE *file = fopen(output_file_name, "w");
    int i, len;
    char *line =
"Chr\tStart\tEnd\tSize\tMean\tL pv\tG pv\tLOH pv\t% L\t%G\t%LOH\t\t\t\t\n";
    len = strlen(line);
    fwrite(line, len, 1, file);
    int chr, bpst, bpe, bps;
    float l2mean, lossp, gainp, lohp, lossprc, gainprc, lohprc;
    int ps;
    float ml, mg, mloh;
    for (i = 0; i < num_seg_regions; i++) {
        chr = seg_chromosomes[i];
        bpst = position_matrix[seg_start[i]];
        bpe = position_matrix[seg_end[i]];
        bps = bpe - bpst + 1;
        l2mean = seg_l2_mean[i];
        lossp = seg_loss_pval[i];
        lossprc = (float) ((float) seg_loss_perc[i] / (float) num_samples);
	gainp = seg_gain_pval[i];
        gainprc = (float) ((float) seg_gain_perc[i] / (float) num_samples);
        lohp = seg_loh_pval[i];
        lohprc = (float) ((float) seg_loh_perc[i] / (float) num_samples);
	ps = seg_size[i];
	ml = mean_loss[i];
	mg = mean_gain[i];
	mloh = mean_loh[i];
	
        if (bps > min_region_bp_size) {
            if (chr < 23) {
                fprintf(file,
	  "%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%f\t%f\t%f\n",
                        chr, bpst, bpe, bps, l2mean, lossp, gainp, lohp,
                        (lossprc), (gainprc), (lohprc), ps, ml, mg, mloh);
            }
            if (chr == 23) {
                fprintf(file,
	    "X\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%f\t%f\t%f\n",
			bpst, bpe, bps, l2mean, lossp, gainp, lohp,
			(lossprc), (gainprc), (lohprc), ps, ml, mg, mloh);
            }
            if (chr == 24) {
                fprintf(file,
	    "Y\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%f\t%f\t%f\n",
			bpst, bpe, bps, l2mean, lossp, gainp, lohp,
			(lossprc), (gainprc), (lohprc), ps, ml, mg, mloh);
            }
            if (chr == 25) {
                fprintf(file,
		"MT\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%f\t%f\t%f\n",
			bpst, bpe, bps, l2mean, lossp, gainp, lohp,
			(lossprc), (gainprc), (lohprc), ps, ml, mg, mloh);
            }
        }
    }
    fclose(file);
}

void compute_pvalue()
{
    float *null_dist_loss = malloc(sizeof(float) * (num_samples + 1));
    float *null_dist_gain = malloc(sizeof(float) * (num_samples + 1));
    float *null_dist_loh = malloc(sizeof(float) * (num_samples + 1));
    int i, j, k, count_loss, count_gain, count_loh;
    float theta;

    for (j = 0; j <= num_samples; j++) {
        null_dist_loss[j] = 0;
        null_dist_gain[j] = 0;
        null_dist_loh[j] = 0;
    }

    for (k = 0; k < bs; k++) {
        count_loss = 0;
        count_gain = 0;
        count_loh = 0;
        for (j = 0; j < num_samples; j++) {
            theta = (float) num_loss_sample[j] / num_seg_regions;
            count_loss += generate_binomial(theta);
            theta = (float) num_gain_sample[j] / num_seg_regions;
            count_gain += generate_binomial(theta);
            if (baf[0] == 1) {
                theta = (float) num_loh_sample[j] / num_seg_regions;
                count_loh += generate_binomial(theta);
            }
        }
        null_dist_loss[count_loss] += 1;
        null_dist_gain[count_gain] += 1;
        if (baf[0] == 1) {
            null_dist_loh[count_loh] += 1;
        }
    }

    for (j = num_samples-1; j >= 0; j--) {
        null_dist_loss[j] += null_dist_loss[j + 1];
        null_dist_gain[j] += null_dist_gain[j + 1];
        null_dist_loh[j] += null_dist_loh[j + 1];
    }
    for (j = 0; j <= num_samples; j++) {
        null_dist_loss[j] = null_dist_loss[j] / bs;
        null_dist_gain[j] = null_dist_gain[j] / bs;
        null_dist_loh[j] = null_dist_loh[j] / bs;
    }
    for (i = 0; i < num_seg_regions; i++) {
        count_loss = 0;
        count_gain = 0;
        count_loh = 0;
        seg_loss_pval[i] = null_dist_loss[seg_loss_perc[i]];
        seg_gain_pval[i] = null_dist_gain[seg_gain_perc[i]];
        if (baf[0] == 1) {
            seg_loh_pval[i] = null_dist_loh[seg_loh_perc[i]];
        } else {
            seg_loh_pval[i] = 1;
        }
    }
}

void compute_matrices()
{
    /* Allocate Memory for the array storing the num of aberrations for each sample */
    num_loss_sample = malloc(sizeof(int) * num_samples);
    num_gain_sample = malloc(sizeof(int) * num_samples);
    num_loh_sample = malloc(sizeof(int) * num_samples);
    
    mean_loss = malloc(sizeof(float) * num_seg_regions);
    mean_gain = malloc(sizeof(float) * num_seg_regions);
    mean_loh = malloc(sizeof(float) * num_seg_regions);
    
    float m;
    int s, i, j, k;
    for (j = 0; j < num_samples; j++) {
        num_loss_sample[j] = 0;
        num_gain_sample[j] = 0;
        num_loh_sample[j] = 0;
	
    }
    for (s = 0; s < num_seg_regions; s++) {
	mean_loss[s]=0;
	mean_gain[s]=0;
	mean_loh[s]=0;
	
        int start = seg_start[s];
        int end = seg_end[s];
        int size = seg_size[s];
        int n_loss = 0;
        int n_gain = 0;
        int n_loh = 0;
        for (j = 0; j < num_samples; j++) {
            float tmp_lrr[size];
            int n_baf = 0;
            int tot_baf = 0;
            k = 0;
            for (i = start; i <= end; i++) {
                tmp_lrr[k] = lrr_matrix[i][j];
                k++;
                if (baf[0] == 1) {
                    if (baf_matrix[i][j] != 2) {
                        tot_baf++;
                        if ((baf_matrix[i][j] < (1 - loh_threshold))
                            || (baf_matrix[i][j] > loh_threshold)) {
                            n_baf++;
                        }
                    }
                }
            }
            m = calc_mean(tmp_lrr, size);
            if (m < loss_threshold) {
                n_loss++;
                num_loss_sample[j]++;
		mean_loss[s]+=m;
            }
            if (m > gain_threshold) {
                n_gain++;
                num_gain_sample[j]++;
		mean_gain[s]+=m;
            }
            if (baf[0] == 1) {
                float p_baf = (float) n_baf / size;
                if (!(m < loss_threshold || m > gain_threshold)
                    && p_baf > loh_frequency) {
		    mean_loh[s]+=m;
                    n_loh++;
                }
            }
        }
     
        seg_loss_perc[s] = n_loss;
        seg_gain_perc[s] = n_gain;
        seg_loh_perc[s] = n_loh;
	if(n_loss > 0){
	    mean_loss[s] = (float) mean_loss[s] / n_loss;
	}
	if(n_gain > 0){
	    mean_gain[s] = (float) mean_gain[s] / n_gain;
	}
	if(n_loh > 0){
	    mean_loh[s] = (float) mean_loh[s] / n_loh;
	}
    }

}

void call_VegaMC()
{
    int c;
    int ntot = 0;
    seg_element *first_seg, *prev_seg;
    for (c = 0; c < num_chromosomes; c++) {
      
        int start_chr = chr_brks_start[c];
        int end_chr = chr_brks_end[c];
        int np = (end_chr - start_chr) + 1;

        float *l2_mean = malloc(sizeof(float) * (np));
        int *start = malloc(sizeof(int) * (np));
        int *end = malloc(sizeof(int) * (np));
        int *start_pos = malloc(sizeof(int) * (np));
        int *size = malloc(sizeof(int) * (np));
        float *std = malloc(sizeof(float) * (num_samples));
        int n_reg = 0;
        float *weight = malloc(sizeof(float) * (num_samples));;
        float weight_sum = num_samples;

        int i, j, k;

        float **tmp_matrix = malloc(sizeof(float *) * num_samples);
        float *tmp_v = malloc(sizeof(float) * np);

        //Rprintf("\t- Analyzing Chromosome %d of %d (composed by %d markers)",
        //       c + 1, num_chromosomes, np);
        for (j = 0; j < num_samples; j++) {
            tmp_matrix[j] = (float *) malloc(sizeof(float) * np);

            for (i = start_chr, k = 0; i <= end_chr; i++, k++) {
                tmp_v[k] = lrr_matrix[i][j];
                tmp_matrix[j][k] = lrr_matrix[i][j];
                if (j == 0) {
                    start_pos[k] = i;
                }
            }

            std[j] = calc_std(tmp_v, np);
            weight[j] = 1.0;

        }
        vegaMC(tmp_matrix, start_pos, start, end, size, l2_mean, &np, &beta,
               std, &n_reg, &num_samples, weight, &weight_sum);
        //Rprintf("\n\t- %d Segmented Regions for Chromsome %d\n\n", n_reg,
        //        c + 1);
        ntot += n_reg;
        
        seg_element *tmp;
        tmp = (struct seg_element *) malloc(sizeof(struct seg_element));
        tmp->chromosome = c;
        tmp->num_regions = n_reg;
        tmp->start = (data_ptr_int) malloc(sizeof(int) * (n_reg));
        tmp->end = (data_ptr_int) malloc(sizeof(int) * (n_reg));
        tmp->size = (data_ptr_int) malloc(sizeof(int) * (n_reg));
        tmp->l2_mean = (data_ptr) malloc(sizeof(float) * (n_reg));
        for (i = 0; i < n_reg; i++) {
            tmp->start[i] = start[i];
            tmp->end[i] = end[i];
            tmp->size[i] = size[i];
            tmp->l2_mean[i] = l2_mean[i];
            if (c == 0) {
                first_seg = tmp;
                prev_seg = first_seg;
            } else {
                prev_seg->next = tmp;
                prev_seg = tmp;
            }
        }
    }
    num_seg_regions = ntot;
    /* Allocation for Segmentation Properties */
    seg_chromosomes = malloc(sizeof(int) * num_seg_regions);
    seg_start = malloc(sizeof(int) * num_seg_regions);
    seg_end = malloc(sizeof(int) * num_seg_regions);
    seg_size = malloc(sizeof(int) * num_seg_regions);
    seg_l2_mean = malloc(sizeof(float) * num_seg_regions);
    seg_loss_pval = malloc(sizeof(float) * num_seg_regions);
    seg_loss_perc = malloc(sizeof(int) * num_seg_regions);
    seg_gain_pval = malloc(sizeof(float) * num_seg_regions);
    seg_gain_perc = malloc(sizeof(int) * num_seg_regions);
    seg_loh_pval = malloc(sizeof(float) * num_seg_regions);
    seg_loh_perc = malloc(sizeof(int) * num_seg_regions);


    /* Use the linked list to store the segmentation for all chromosomes */
    seg_element *tmp;
    tmp = first_seg;
    int i;
    int j = 0;
    for (c = 0; c < num_chromosomes; c++) {
        for (i = 0; i < (tmp->num_regions); i++) {
            seg_start[j] = tmp->start[i];
            seg_end[j] = tmp->end[i];
            seg_size[j] = tmp->size[i];
            seg_l2_mean[j] = tmp->l2_mean[i];
            seg_chromosomes[j] = (tmp->chromosome) + 1;
            j++;
        }
        tmp = tmp->next;
    }

}

void load_data()
{

    int ns = 0, np = 0;
    char line[MAX_INPUT_LINE];
    char *sep = "\t";
    char *elem;
    probe *first_probe, *prev_probe, *temp_probe;
    float *data_tmp;
  
    FILE *file = fopen(dataset, "r");
    num_chromosomes = 0;
  
    if (fgets(line, sizeof(line), file) != 0) {
        elem = strtok(line, sep);
        while (elem) {
            ns++;
            elem = strtok('\0', sep);
        }
    }
    num_samples = ns - 3;
    
    data_tmp = malloc(sizeof(float) * (num_samples));
    int prev_chr = 0;
    while (fgets(line, sizeof(line), file) != 0) {
        int i = 0, position, chr;
        elem = strtok(line, sep);
        while (elem) {
            if (i == 1) {
                if (strcmp(elem, "X") == 0) {
                    chr = 23;
                } else {
                    if (strcmp(elem, "Y") == 0) {
                        chr = 24;
                    } else {
                        if (strcmp(elem, "MT") == 0) {
                            chr = 25;
                        } else {
                            chr = atoi(elem);
                        }
                    }
                }
                if (prev_chr != chr) {
                    num_chromosomes++;
                    prev_chr = chr;
                }
            }
            if (i == 2) {
                position = atoi(elem);
            }
            if (i > 2 && i < (num_samples + 3)) {
                data_tmp[i - 3] = atof(elem);
            }
            elem = strtok('\0', sep);
            i++;
        }
        
        int j;
        /* Create a new probe object */
        temp_probe = (struct probe *) malloc(sizeof(struct probe));
        temp_probe->data = (data_ptr) malloc(sizeof(float) * (num_samples));
        temp_probe->chromosome = chr;
        temp_probe->position = position;
        temp_probe->next = NULL;
        for (j = 0; j < num_samples; j++) {
            temp_probe->data[j] = data_tmp[j];
        }
        if (np == 0) {
            first_probe = temp_probe;
            prev_probe = first_probe;
        } else {
            prev_probe->next = temp_probe;
            prev_probe = temp_probe;
        }
        np++;
        
    }
    
    free(data_tmp);
    fclose(file);

    num_probes = np;
    /* Initialize the data structures */
    lrr_matrix = (float **) malloc(sizeof(float *) * num_probes);
    if (baf[0] == 1) {
        num_samples = num_samples / 2;
        baf_matrix = (float **) malloc(sizeof(float *) * num_probes);
    }
    position_matrix = malloc(sizeof(int) * num_probes);
    chromosome_matrix = malloc(sizeof(int) * num_probes);
    chr_brks_start = malloc(sizeof(int) * num_chromosomes);
    chr_brks_end = malloc(sizeof(int) * num_chromosomes);
    
    temp_probe = first_probe;
    prev_chr = temp_probe->chromosome;
    int i, j, k;
    int start = 0, chr_index = 0;
    if (baf[0] == 1) {
        /* Fill the data_matrix, the positions_matrix, the  chr_brks_start and chr_brks_end */
        for (i = 0; i < num_probes; i++) {
            lrr_matrix[i] = (float *) malloc(sizeof(float) * num_samples);
            baf_matrix[i] = (float *) malloc(sizeof(float) * num_samples);
            position_matrix[i] = temp_probe->position;
            chromosome_matrix[i] = temp_probe->chromosome;
            for (j = 0, k = 0; j < (num_samples * 2); j += 2, k++) {
                lrr_matrix[i][k] = temp_probe->data[j];
                baf_matrix[i][k] = temp_probe->data[j + 1];

            }

            if (prev_chr != temp_probe->chromosome) {
                chr_brks_start[chr_index] = start;
                chr_brks_end[chr_index] = i - 1;
                start = i;
                prev_chr = temp_probe->chromosome;
                chr_index++;
            }
            temp_probe = temp_probe->next;
        }
        chr_brks_start[num_chromosomes - 1] = start;
        chr_brks_end[num_chromosomes - 1] = num_probes - 1;
    } else {
        for (i = 0; i < num_probes; i++) {
            lrr_matrix[i] = (float *) malloc(sizeof(float) * num_samples);
            position_matrix[i] = temp_probe->position;
            chromosome_matrix[i] = temp_probe->chromosome;
            for (j = 0; j < num_samples; j++) {
                lrr_matrix[i][j] = temp_probe->data[j];
            }
            if (prev_chr != temp_probe->chromosome) {
                chr_brks_start[chr_index] = start;
                chr_brks_end[chr_index] = i - 1;
                start = i;
                prev_chr = temp_probe->chromosome;
                chr_index++;
            }
            temp_probe = temp_probe->next;
        }
        chr_brks_start[num_chromosomes - 1] = start;
        chr_brks_end[num_chromosomes - 1] = num_probes - 1;
    }
}

int read_params(char *filename)
{
    int index = 0;
    char line[MAX_INPUT_LINE];
    char *sep = "\t";
    char *elem, *brkt;
    FILE *file = fopen(filename, "r");

    while ((fgets(line, sizeof(line), file)) != 0) {
        index++;
        elem = strtok(line, sep);
        switch (index) {
        case 1:
            elem = strtok(NULL, sep);
            beta = atof(elem);
            break;
        case 2:
            elem = strtok(NULL, sep);
            min_region_probe_size = atoi(elem);
            break;
        case 3:
            elem = strtok(NULL, sep);
            min_region_bp_size = atoi(elem);
            break;
        case 4:
            elem = strtok(NULL, sep);
            loss_threshold = atof(elem);
            break;
        case 5:
            elem = strtok(NULL, sep);
            gain_threshold = atof(elem);
            break;
        case 6:
            elem = strtok(NULL, sep);
            baf[0] = atoi(elem);
            break;
        case 7:
            elem = strtok(NULL, sep);
            loh_threshold = atof(elem);
            break;
        case 8:
            elem = strtok(NULL, sep);
            loh_frequency = atof(elem);
            break;
        case 9:
            elem = strtok(NULL, sep);
            bs = atoi(elem);
            break;
        case 10:
            elem = strtok(NULL, sep);
            pval_threshold = atof(elem);
            break;
        }
    }
    fclose(file);
    return 0;
}

/*--- Start of Support mean and std Functions ---*/
float calc_mean(float vals[], int n_elem)
{
    int i;
    float sum = 0;
    for (i = 0; i < n_elem; i++) {
        sum += vals[i];
    }
    return sum / n_elem;
}

float calc_std(float vals[], int n_elem)
{
    int i;
    float sum = 0;
    float m = calc_mean(vals, n_elem);
    for (i = 0; i < n_elem; i++) {
        sum += pow(vals[i] - m, 2);
    }
    return sqrt(sum / (n_elem - 1));

}

int generate_binomial(float theta)
{
    float r = (float) rand() / RAND_MAX;
    int r2 = (int) (r * ROUND_CNST);
    r = (float) r2 / ROUND_CNST;
    if (r <= theta) {
        return 1;
    } else {
        return 0;
    }
}

/*--- End of Support mean and std Functions ---*/
