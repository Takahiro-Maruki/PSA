// Updated on 11/24/21

/* 

Program PSA.cpp to analyze population structure from high-throughput sequencing data of 
diploid individuals from multiple populations. The genotype frequencies necessary for the 
analysis are estimated beforehand by GFE in the F mode.  The population structure is analyzed 
only when a site is polymorphic in the total population.  Fixation indices are calculated only 
from populations with ML estimates and effective sample size equal to or greater than ten.  
Results at all sites are shown.  Statistical significance of the polymorphism in a deme is 
examined when deciding to add a new allele.  The input file contains the error-rate estimate.  
The mean within-population heterozygosity Hs and heterozygosity in the total population Ht are 
reported in the output.  Hs and Ht are also printed out at sites with more than two alleles.  
The total coverage is simply the sum of depths of coverage over the populations.  The minor-allele 
frequency estimate in the total population is reported in the output.  

*/

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <algorithm>
using namespace std;

int main(int argc, char *argv[])
{
	// Default values of the options
	const char* in_file_name = {"In_PSA.txt"};
	const char* out_file_name = {"Out_PSA.txt"};
	double min_Ni = 10.0;
	double cv = 5.991;
	int print_help = 0;
	
	int argz = 1; // argument counter

	// Read specified settings
	while( (argz<argc) && (argv[argz][0] == '-') ) {
		if (strcmp(argv[argz], "-h") == 0) {
			print_help = 1;
		} else if (strcmp(argv[argz], "-in") == 0) {
			in_file_name = argv[++argz];
		} else if (strcmp(argv[argz], "-out") == 0) {
			out_file_name = argv[++argz];
		} else if (strcmp(argv[argz], "-min_Ni") == 0) {
			sscanf(argv[++argz], "%lf", &min_Ni);
		} else if (strcmp(argv[argz], "-cv") == 0) {
			sscanf(argv[++argz], "%lf", &cv);
		} else {
			fprintf(stderr, "unknown option %s\n", argv[argz]);
			print_help = 1;
			break;
		}
		argz++;
	}
	if (print_help) { // print error/usage message ?
		fprintf(stderr, "USAGE: %s {<options>}\n", argv[0]);
		fprintf(stderr, "	options:\n");
		fprintf(stderr, "	-h: print the usage message\n");
		fprintf(stderr, "	-in <s>: specify the input file name\n");
		fprintf(stderr, "       -out <s>: specify the output file name\n");
		fprintf(stderr, "	-min_Ni <f>: specify the minimum effective number of sampled individuals required in a deme\n");
		fprintf(stderr, "       -cv <f>: specify the chi-square critical value for the polymorphism test\n");
		exit(1);
	}

	string line; // String buffer
	
	ifstream inputFile(in_file_name); // Try to open the input file
	if ( !inputFile.is_open() ) { // Exit on failure
		fprintf(stderr, "Cannot open %s for reading.\n", in_file_name);
		exit(1);
	}
	
	// Read the header
	string h_scaf, h_site, h_ref_nuc;
	vector <string> pop_info; // Stores population-specific labels.
	getline(inputFile, line);
	istringstream ss(line);
	ss >> h_scaf >> h_site >> h_ref_nuc;
	string str; // Temporarily stores population-specific labels
	pop_info.clear();
	while (true) {
		ss >> str;
		pop_info.push_back(str);
		if ( ss.eof() ) {
			break;
		}
	}
	int num_pops = (int)pop_info.size()/9;
	printf("%d populations to be analyzed\n", num_pops);

	FILE *outstream;

	// Open the output file
	outstream = fopen(out_file_name, "w");
	if (outstream == NULL ) { // Exit on failure
		fprintf(stderr, "Cannot open %s for writing.\n", out_file_name); 
		exit(1);
	}
	
	// Print out the field names
	fprintf(outstream, "scaffold\tsite\tref_nuc\ttot_cov\tnumber_of_alleles\tne_pops\tF\ttheta\tf\tHs\tHt\tMAF\n");
	// printf("scaffold\tsite\tref_nuc\ttot_cov\tnumber_of_alleles\tne_pops\tF\ttheta\tf\tHs\tHt\tMAF\n");
	
	// Read the main data
	string scaffold, ref_nuc, n1[num_pops+1], n2[num_pops+1], s_Ni[num_pops+1], s_best_p[num_pops+1], s_best_q[num_pops+1], best_error, s_best_H[num_pops+1], s_pol_llstat[num_pops+1];
	double Ni[num_pops+1], best_x[num_pops+1], best_H[num_pops+1], pol_llstat[num_pops+1]; 
	int site, pop_cov[num_pops+1], pg, num_alleles, ag;
	int tot_cov;		// total coverage (sum of the coverage across the populations)
	int ne_pops;            // total number of populations with data
	vector <string> alleles;   // store allele identities
	double sum_Ni, sum_sq_Ni, w_sum_af, w_sum_H, w_mean_af, mean_Ni, ss_af, var_af, w_mean_H, n_c, S_1, S_2, S_3, F, theta, f, Hs, Ht;
	double ma_w_mean_af[5];		// Weighted frequency of each of the alleles at multiallelic sites 
	double sum_S_1, sum_S_2, sum_S_3;
	double maf;

	while ( getline(inputFile, line) ) {
		istringstream ss(line);
		alleles.clear();
		tot_cov = 0;
		ne_pops = 0;
		ss >> scaffold >> site >> ref_nuc;
		for (pg = 1; pg <= num_pops; pg++) {
			ss >> n1[pg] >> n2[pg] >> pop_cov[pg] >> s_Ni[pg] >> s_best_p[pg] >> s_best_q[pg] >> best_error >> s_best_H[pg] >> s_pol_llstat[pg];
			tot_cov = tot_cov + pop_cov[pg];
			// printf("site: %d\tpop: %d\tn1: %s\tn2: %s\n", site, pg, n1[pg].c_str(), n2[pg].c_str());
			// fprintf(outstream, "site: %d\tpop: %d\tn1: %s\tn2: %s\n", site, pg, n1[pg].c_str(), n2[pg].c_str());
			if (n1[pg] != "NA") {
				Ni[pg] = atof(s_Ni[pg].c_str());
				if (Ni[pg] >= min_Ni) {
					ne_pops = ne_pops + 1;
                                	if ( find(alleles.begin(), alleles.end(), n1[pg]) == alleles.end() ) {
						// printf("%s\n", n1[pg].c_str());
                                        	alleles.push_back(n1[pg]);
                                	}
                        		if (n2[pg] != "NA") {
						pol_llstat[pg] = atof(s_pol_llstat[pg].c_str());
						if (pol_llstat[pg] > cv) {
                                			if ( find(alleles.begin(), alleles.end(), n2[pg]) == alleles.end() ) {
								// printf("%s\n", n2[pg].c_str());
                                        			alleles.push_back(n2[pg]);
							}
                                		}
					}
				}
                        }
		}
		// Count the number of alleles segregating in the population sample
		num_alleles = alleles.size();
		/*
		printf("site: %d\tnum_alleles: %d\n", site, num_alleles);
		fprintf(outstream, "site: %d\tnum_alleles: %d\n", site, num_alleles);
		for (ag=0; ag<num_alleles; ag++) {
			printf("%s\n", alleles.at(ag).c_str());
		}
		*/
		if (num_alleles == 2 && ne_pops > 1) {	// Calculate the fixations indixes for bi-allelic polymorphism
			sum_Ni = 0.0;
			sum_sq_Ni = 0.0;
			w_sum_af = 0.0;
			w_sum_H = 0.0;
			for (pg = 1; pg <= num_pops; pg++) { // Calculate the fixation indices using frequencies of the first allele stored in the vector
				if (n1[pg] != "NA" && Ni[pg] >= min_Ni) {	// Examine the population only when there are ML estimates with Ni equal to or greater than the specified value at the site
					if ( n1[pg] == alleles.at(0) ) {
						sum_Ni = sum_Ni + Ni[pg];
						sum_sq_Ni = sum_sq_Ni + Ni[pg]*Ni[pg];
						best_x[pg] = atof(s_best_p[pg].c_str());
						w_sum_af = w_sum_af + Ni[pg]*best_x[pg];
						best_H[pg] = atof(s_best_H[pg].c_str());
						w_sum_H = w_sum_H + Ni[pg]*best_H[pg];
					} else if ( n2[pg] == alleles.at(0) ) {
                                                sum_Ni = sum_Ni + Ni[pg];
                                                sum_sq_Ni = sum_sq_Ni + Ni[pg]*Ni[pg];
                                                best_x[pg] = atof(s_best_q[pg].c_str());
                                                w_sum_af = w_sum_af + Ni[pg]*best_x[pg];
                                                best_H[pg] = atof(s_best_H[pg].c_str());
                                                w_sum_H = w_sum_H + Ni[pg]*best_H[pg];
					} else {
						sum_Ni = sum_Ni + Ni[pg];
						sum_sq_Ni = sum_sq_Ni + Ni[pg]*Ni[pg];
						best_x[pg] = 0.0;
						w_sum_af = w_sum_af + Ni[pg]*best_x[pg];
						w_sum_H = w_sum_H + Ni[pg]*0.0;
					}
				}
			}
			w_mean_af = w_sum_af/sum_Ni;
			Ht = 2.0*w_mean_af*(1.0-w_mean_af);
			if (w_mean_af <= (double)1.0-w_mean_af) {
				maf = w_mean_af;
			} else {
				maf = (double)1.0-w_mean_af;
			}
			mean_Ni = sum_Ni/ne_pops;
			ss_af = 0.0;
			Hs = 0.0;
			for (pg = 1; pg <= num_pops; pg++) {
				if (n1[pg] != "NA" && Ni[pg] >= min_Ni) {
					Hs = Hs + Ni[pg]*2.0*best_x[pg]*(1.0-best_x[pg]);
					ss_af = ss_af + Ni[pg]*pow( best_x[pg]-w_mean_af, 2.0);
				}
			}
			Hs = Hs/sum_Ni;
			var_af = ( 1.0/( (ne_pops -1)*mean_Ni ) )*ss_af;
			w_mean_H = w_sum_H/sum_Ni;
			n_c = ( 1.0/((double)ne_pops -1.0) )*(sum_Ni - sum_sq_Ni/sum_Ni);
			S_1 = var_af - ( 1.0/(mean_Ni-1.0) )*( w_mean_af*(1.0-w_mean_af) - ( (ne_pops-1.0)/(double)ne_pops )*var_af - w_mean_H/4.0 );
			S_2 = w_mean_af*(1.0-w_mean_af) - ( mean_Ni/( (double)ne_pops*(mean_Ni-1.0) ) )*( ( ( (double)ne_pops*(mean_Ni-n_c) )/mean_Ni )*w_mean_af*(1.0-w_mean_af) - (1.0/mean_Ni)*( mean_Ni-1.0 + ( (double)ne_pops -1.0 )*(mean_Ni-n_c) )*var_af - ( ( (double)ne_pops*(mean_Ni-n_c) )/( 4.0*mean_Ni*n_c ) )*w_mean_H );
			S_3 = ( n_c/(2.0*mean_Ni) )*w_mean_H;
			fprintf(outstream, "%s\t%d\t%s\t%d\t%d\t%d\t", scaffold.c_str(), site, ref_nuc.c_str(), tot_cov, num_alleles, ne_pops);
			// printf("%s\t%d\t%s\t%d\t%d\t%d\t", scaffold.c_str(), site, ref_nuc.c_str(), tot_cov, num_alleles, ne_pops);
			if (S_2 > 1.0e-8 || S_2 < -1.0e-8) {
				F = 1.0 - S_3/S_2;
				theta = S_1/S_2;
				fprintf(outstream, "%f\t%f\t", F, theta);
				// printf("%f\t%f\t", F, theta);
				if (theta < 1.0) {
					f = (F-theta)/(1.0-theta);
					fprintf(outstream, "%f\t", f);
						// printf("%f\t", f);
				} else {
					fprintf(outstream, "NA\t");
					// printf("NA\t");
				}
			} else {
				fprintf(outstream, "NA\tNA\tNA\t");
				// printf("NA\tNA\tNA\t");
			}
			fprintf(outstream, "%f\t%f\t%f\n", Hs, Ht, maf);
			// printf("%f\t%f\t%f\n", Hs, Ht, maf);			
		} else if (num_alleles > 2 && ne_pops > 1) { // Calculate the fixations indixes at sites with more than two alleles segregating
			sum_S_1 = 0.0;
			sum_S_2 = 0.0;
			sum_S_3 = 0.0;
			for (ag=1; ag<=num_alleles; ag++) {
				sum_Ni = 0.0;
				sum_sq_Ni = 0.0;
				w_sum_af = 0.0;
				w_sum_H = 0.0;
				for (pg = 1; pg <= num_pops; pg++) { // Calculate the fixation indices using frequencies of each of the alleles stored in the vector
					if (n1[pg] != "NA" && Ni[pg] >= min_Ni) {	// Examine the population only when there are ML estimates with Ni equal to or greater than the specified value at the site
						if ( n1[pg] == alleles.at(ag-1) ) {
							sum_Ni = sum_Ni + Ni[pg];
							sum_sq_Ni = sum_sq_Ni + Ni[pg]*Ni[pg];
							best_x[pg] = atof(s_best_p[pg].c_str());
							w_sum_af = w_sum_af + Ni[pg]*best_x[pg];
							best_H[pg] = atof(s_best_H[pg].c_str());
							w_sum_H = w_sum_H + Ni[pg]*best_H[pg];
						} else if ( n2[pg] == alleles.at(ag-1) ) {
							sum_Ni = sum_Ni + Ni[pg];
							sum_sq_Ni = sum_sq_Ni + Ni[pg]*Ni[pg];
							best_x[pg] = atof(s_best_q[pg].c_str());
							w_sum_af = w_sum_af + Ni[pg]*best_x[pg];
							best_H[pg] = atof(s_best_H[pg].c_str());
							w_sum_H = w_sum_H + Ni[pg]*best_H[pg];
						} else {
							sum_Ni = sum_Ni + Ni[pg];
							sum_sq_Ni = sum_sq_Ni + Ni[pg]*Ni[pg];
							best_x[pg] = 0.0;
							w_sum_af = w_sum_af + Ni[pg]*best_x[pg];
							w_sum_H = w_sum_H + Ni[pg]*0.0;
						}
					}
				}
				ma_w_mean_af[ag] = w_sum_af/sum_Ni;
				mean_Ni = sum_Ni/ne_pops;
				ss_af = 0.0;
				for (pg = 1; pg <= num_pops; pg++) {
					if (n1[pg] != "NA" && Ni[pg] >= min_Ni) {
						ss_af = ss_af + Ni[pg]*pow( best_x[pg]-ma_w_mean_af[ag], 2.0);
					}
				}					
				var_af = ( 1.0/( (ne_pops -1)*mean_Ni ) )*ss_af;
				w_mean_H = w_sum_H/sum_Ni;
				n_c = ( 1.0/((double)ne_pops -1.0) )*(sum_Ni - sum_sq_Ni/sum_Ni);
				S_1 = var_af - ( 1.0/(mean_Ni-1.0) )*( ma_w_mean_af[ag]*(1.0-ma_w_mean_af[ag]) - ( (ne_pops-1.0)/(double)ne_pops )*var_af - w_mean_H/4.0 );
				S_2 = ma_w_mean_af[ag]*(1.0-ma_w_mean_af[ag]) - ( mean_Ni/( (double)ne_pops*(mean_Ni-1.0) ) )*( ( ( (double)ne_pops*(mean_Ni-n_c) )/mean_Ni )*ma_w_mean_af[ag]*(1.0-ma_w_mean_af[ag]) - (1.0/mean_Ni)*( mean_Ni-1.0 + ( (double)ne_pops -1.0 )*(mean_Ni-n_c) )*var_af - ( ( (double)ne_pops*(mean_Ni-n_c) )/( 4.0*mean_Ni*n_c ) )*w_mean_H );
				S_3 = ( n_c/(2.0*mean_Ni) )*w_mean_H;
				sum_S_1 = sum_S_1 + S_1;
				sum_S_2 = sum_S_2 + S_2;
				sum_S_3 = sum_S_3 + S_3;
			}       // end of the loop over the alleles
			// Calculate Hs
			Hs = 0.0;
			for (pg = 1; pg <= num_pops; pg++) { // Calculate Hs
				if (n1[pg] != "NA" && Ni[pg] >= min_Ni) {       // Examine the population only when there are ML estimates with Ni equal to or greater than the specified value at the site
					Hs = Hs + Ni[pg]*2.0*best_x[pg]*(1.0-best_x[pg]);
				}
			}
			Hs = Hs/sum_Ni;
			// Calculate Ht and maf
			if (num_alleles == 3) {
				Ht = 1.0 - pow(ma_w_mean_af[1], 2.0) - pow(ma_w_mean_af[2], 2.0) - pow(ma_w_mean_af[3], 2.0);
				if (ma_w_mean_af[1] <= ma_w_mean_af[2]) {
					maf = ma_w_mean_af[1];
				} else {
					maf = ma_w_mean_af[2];
				}
				if (ma_w_mean_af[3] < maf) {
					maf = ma_w_mean_af[3];
				}
			} else if (num_alleles == 4) {
				Ht = 1.0 - pow(ma_w_mean_af[1], 2.0) - pow(ma_w_mean_af[2], 2.0) - pow(ma_w_mean_af[3], 2.0) - pow(ma_w_mean_af[4], 2.0);
				if (ma_w_mean_af[1] <= ma_w_mean_af[2]) {
                                        maf = ma_w_mean_af[1];
                                } else {
                                        maf = ma_w_mean_af[2];
                                }
                                if (ma_w_mean_af[3] < maf) {
                                        maf = ma_w_mean_af[3];
                                }
				if (ma_w_mean_af[4] < maf) {
                                        maf = ma_w_mean_af[4];
                                }
			}
			fprintf(outstream, "%s\t%d\t%s\t%d\t%d\t%d\t", scaffold.c_str(), site, ref_nuc.c_str(), tot_cov, num_alleles, ne_pops);
                        // printf("%s\t%d\t%s\t%d\t%d\t%d\t", scaffold.c_str(), site, ref_nuc.c_str(), tot_cov, num_alleles, ne_pops);
                        if (sum_S_2 > 1.0e-8 || sum_S_2 < -1.0e-8) {
                        	F = 1.0 - sum_S_3/sum_S_2;
                        	theta = sum_S_1/sum_S_2;
                        	fprintf(outstream, "%f\t%f\t", F, theta);
                        	// printf("%f\t%f\t", F, theta);
                        	if (theta < 1.0) {
                        		f = (F-theta)/(1.0-theta);
                                	fprintf(outstream, "%f\t", f);
                                	// printf("%f\t", f);
                        	} else {
                                	fprintf(outstream, "NA\t");
                                	// printf("NA\t");
                        	}
			} else {
                                fprintf(outstream, "NA\tNA\tNA\n");
                                // printf("NA\tNA\tNA\n");
                        }
			fprintf(outstream, "%f\t%f\t%f\n", Hs, Ht, maf);
			// printf("%f\t%f\t%f\n", Hs, Ht, maf);
		} else {		// Print out the information at other sites 
			fprintf(outstream, "%s\t%d\t%s\t%d\t%d\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), tot_cov, num_alleles, ne_pops);
			// printf("%s\t%d\t%s\t%d\t%d\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), tot_cov, num_alleles, ne_pops);
		}						
	} 		

	return 0;
}
