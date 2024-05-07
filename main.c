#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "kvec.h"
#include "sys.h"
#include "paf.h"
#include "sdict.h"
#include "miniasm.h"

#define ZA_VERSION "1.0"

int main(int argc, char *argv[])
{
	ma_opt_t opt;
	int i, c, stage = 100, no_first = 0, no_second = 0, bi_dir = 1, o_set = 0, no_cont = 0;
	sdict_t *d, *excl = 0;
	ma_sub_t *sub = 0;
	ma_hit_t *hit;
	size_t n_hits;
	float cov = 40.0;
	char *fn_reads = 0, *outfmt = "ug";

	ma_opt_init(&opt);
	while ((c = getopt(argc, argv, "n:m:s:c:S:i:d:g:o:h:I:r:f:e:p:12VBRbF:")) >= 0) {
		if (c == 'm') opt.min_match = atoi(optarg);
		else if (c == 'i') opt.min_iden = atof(optarg);
		else if (c == 's') opt.min_span = atoi(optarg);
		else if (c == 'c') opt.min_dp = atoi(optarg);
		else if (c == 'o') opt.min_ovlp = atoi(optarg), o_set = 1;
		else if (c == 'S') stage = atoi(optarg);
		else if (c == 'd') opt.bub_dist = atoi(optarg);
		else if (c == 'g') opt.gap_fuzz = atoi(optarg);
		else if (c == 'h') opt.max_hang = atoi(optarg);
		else if (c == 'I') opt.int_frac = atof(optarg);
		else if (c == 'e') opt.max_ext = atoi(optarg);
		else if (c == 'f') fn_reads = optarg;
		else if (c == 'p') outfmt = optarg;
		else if (c == '1') no_first = 1;
		else if (c == '2') no_second = 1;
		else if (c == 'n') opt.n_rounds = atoi(optarg) - 1;
		else if (c == 'B') bi_dir = 1;
		else if (c == 'b') bi_dir = 0;
		else if (c == 'R') no_cont = 1;
		else if (c == 'F') opt.final_ovlp_drop_ratio = atof(optarg);
		else if (c == 'V') {
			printf("%s\n", ZA_VERSION);
			return 0;
		} else if (c == 'r') {
			char *s;
			opt.max_ovlp_drop_ratio = strtod(optarg, &s);
			if (*s == ',') opt.min_ovlp_drop_ratio = strtod(s + 1, &s);
		}
	}
	if (o_set == 0) opt.min_ovlp = opt.min_span;
	if (argc == optind) {
		fprintf(stderr, "Usage: miniasm [options] <in.paf>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  Pre-selection:\n");
		fprintf(stderr, "    -R          prefilter clearly contained reads (2-pass required)\n");
		fprintf(stderr, "    -m INT      min match length [%d]\n", opt.min_match);
		fprintf(stderr, "    -i FLOAT    min identity [%.2g]\n", opt.min_iden);
		fprintf(stderr, "    -s INT      min span [%d]\n", opt.min_span);
		fprintf(stderr, "    -c INT      min coverage [%d]\n", opt.min_dp);
		fprintf(stderr, "  Overlap:\n");
		fprintf(stderr, "    -o INT      min overlap [same as -s]\n");
		fprintf(stderr, "    -h INT      max over hang length [%d]\n", opt.max_hang);
		fprintf(stderr, "    -I FLOAT    min end-to-end match ratio [%.2g]\n", opt.int_frac);
		fprintf(stderr, "  Layout:\n");
		fprintf(stderr, "    -g INT      max gap differences between reads for trans-reduction [%d]\n", opt.gap_fuzz);
		fprintf(stderr, "    -d INT      max distance for bubble popping [%d]\n", opt.bub_dist);
		fprintf(stderr, "    -e INT      small unitig threshold [%d]\n", opt.max_ext);
		fprintf(stderr, "    -f FILE     read sequences []\n");
		fprintf(stderr, "    -n INT      rounds of short overlap removal [%d]\n", opt.n_rounds + 1);
		fprintf(stderr, "    -r FLOAT[,FLOAT]\n");
		fprintf(stderr, "                max and min overlap drop ratio [%.2g,%.2g]\n", opt.max_ovlp_drop_ratio, opt.min_ovlp_drop_ratio);
		fprintf(stderr, "    -F FLOAT    aggressive overlap drop ratio in the end [%.2g]\n", opt.final_ovlp_drop_ratio);
		fprintf(stderr, "  Miscellaneous:\n");
		fprintf(stderr, "    -p STR      output information: bed, paf, sg or ug [%s]\n", outfmt);
//		fprintf(stderr, "    -B          only one direction of an arc is present in input PAF\n"); // deprecated; for backward compatibility
		fprintf(stderr, "    -b          both directions of an arc are present in input\n");
		fprintf(stderr, "    -1          skip 1-pass read selection\n");
		fprintf(stderr, "    -2          skip 2-pass read selection\n");
		fprintf(stderr, "    -V          print version number\n");
		fprintf(stderr, "\nSee miniasm.1 for detailed description of the command-line options.\n");
		return 1;
	}

	sys_init();

    