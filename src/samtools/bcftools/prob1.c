#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <assert.h>
#include "prob1.h"

#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 16384)

#define MC_MAX_EM_ITER 16
#define MC_EM_EPS 1e-4
#define MC_DEF_INDEL 0.15

unsigned char seq_nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4 /*'-'*/, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

struct __bcf_p1aux_t {
	int n, M, n1, is_indel;
	uint8_t *ploidy; // haploid or diploid ONLY
	double *q2p, *pdg; // pdg -> P(D|g)
	double *phi, *phi_indel;
	double *z, *zswap; // aux for afs
	double *z1, *z2, *phi1, *phi2; // only calculated when n1 is set
	double t, t1, t2;
	double *afs, *afs1; // afs: accumulative AFS; afs1: site posterior distribution
	const uint8_t *PL; // point to PL
	int PL_len;
};

void bcf_p1_indel_prior(bcf_p1aux_t *ma, double x)
{
	int i;
	for (i = 0; i < ma->M; ++i)
		ma->phi_indel[i] = ma->phi[i] * x;
	ma->phi_indel[ma->M] = 1. - ma->phi[ma->M] * x;
}

static void init_prior(int type, double theta, int M, double *phi)
{
	int i;
	if (type == MC_PTYPE_COND2) {
		for (i = 0; i <= M; ++i)
			phi[i] = 2. * (i + 1) / (M + 1) / (M + 2);
	} else if (type == MC_PTYPE_FLAT) {
		for (i = 0; i <= M; ++i)
			phi[i] = 1. / (M + 1);
	} else {
		double sum;
		for (i = 0, sum = 0.; i < M; ++i)
			sum += (phi[i] = theta / (M - i));
		phi[M] = 1. - sum;
	}
}

void bcf_p1_init_prior(bcf_p1aux_t *ma, int type, double theta)
{
	init_prior(type, theta, ma->M, ma->phi);
	bcf_p1_indel_prior(ma, MC_DEF_INDEL);
}

void bcf_p1_init_subprior(bcf_p1aux_t *ma, int type, double theta)
{
	if (ma->n1 <= 0 || ma->n1 >= ma->M) return;
	init_prior(type, theta, 2*ma->n1, ma->phi1);
	init_prior(type, theta, 2*(ma->n - ma->n1), ma->phi2);
}

int bcf_p1_read_prior(bcf_p1aux_t *ma, const char *fn)
{
	gzFile fp;
	kstring_t s;
	kstream_t *ks;
	long double sum;
	int dret, k;
	memset(&s, 0, sizeof(kstring_t));
	fp = strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	ks = ks_init(fp);
	memset(ma->phi, 0, sizeof(double) * (ma->M + 1));
	while (ks_getuntil(ks, '\n', &s, &dret) >= 0) {
		if (strstr(s.s, "[afs] ") == s.s) {
			char *p = s.s + 6;
			for (k = 0; k <= ma->M; ++k) {
				int x;
				double y;
				x = strtol(p, &p, 10);
				if (x != k && (errno == EINVAL || errno == ERANGE)) return -1;
				++p;
				y = strtod(p, &p);
				if (y == 0. && (errno == EINVAL || errno == ERANGE)) return -1;
				ma->phi[ma->M - k] += y;
			}
		}
	}
	ks_destroy(ks);
	gzclose(fp);
	free(s.s);
	for (sum = 0., k = 0; k <= ma->M; ++k) sum += ma->phi[k];
	fprintf(stderr, "[prior]");
	for (k = 0; k <= ma->M; ++k) ma->phi[k] /= sum;
	for (k = 0; k <= ma->M; ++k) fprintf(stderr, " %d:%.3lg", k, ma->phi[ma->M - k]);
	fputc('\n', stderr);
	for (sum = 0., k = 1; k < ma->M; ++k) sum += ma->phi[ma->M - k] * (2.* k * (ma->M - k) / ma->M / (ma->M - 1));
	fprintf(stderr, "[%s] heterozygosity=%lf, ", __func__, (double)sum);
	for (sum = 0., k = 1; k <= ma->M; ++k) sum += k * ma->phi[ma->M - k] / ma->M;
	fprintf(stderr, "theta=%lf\n", (double)sum);
	bcf_p1_indel_prior(ma, MC_DEF_INDEL);
	return 0;
}

bcf_p1aux_t *bcf_p1_init(int n, uint8_t *ploidy)
{
	bcf_p1aux_t *ma;
	int i;
	ma = calloc(1, sizeof(bcf_p1aux_t));
	ma->n1 = -1;
	ma->n = n; ma->M = 2 * n;
	if (ploidy) {
		ma->ploidy = malloc(n);
		memcpy(ma->ploidy, ploidy, n);
		for (i = 0, ma->M = 0; i < n; ++i) ma->M += ploidy[i];
		if (ma->M == 2 * n) {
			free(ma->ploidy);
			ma->ploidy = 0;
		}
	}
	ma->q2p = calloc(256, sizeof(double));
	ma->pdg = calloc(3 * ma->n, sizeof(double));
	ma->phi = calloc(ma->M + 1, sizeof(double));
	ma->phi_indel = calloc(ma->M + 1, sizeof(double));
	ma->phi1 = calloc(ma->M + 1, sizeof(double));
	ma->phi2 = calloc(ma->M + 1, sizeof(double));
	ma->z = calloc(ma->M + 1, sizeof(double));
	ma->zswap = calloc(ma->M + 1, sizeof(double));
	ma->z1 = calloc(ma->M + 1, sizeof(double)); // actually we do not need this large
	ma->z2 = calloc(ma->M + 1, sizeof(double));
	ma->afs = calloc(ma->M + 1, sizeof(double));
	ma->afs1 = calloc(ma->M + 1, sizeof(double));
	for (i = 0; i < 256; ++i)
		ma->q2p[i] = pow(10., -i / 10.);
	bcf_p1_init_prior(ma, MC_PTYPE_FULL, 1e-3); // the simplest prior
	return ma;
}

int bcf_p1_set_n1(bcf_p1aux_t *b, int n1)
{
	if (n1 == 0 || n1 >= b->n) return -1;
	if (b->M != b->n * 2) {
		fprintf(stderr, "[%s] unable to set `n1' when there are haploid samples.\n", __func__);
		return -1;
	}
	b->n1 = n1;
	return 0;
}

void bcf_p1_destroy(bcf_p1aux_t *ma)
{
	if (ma) {
		free(ma->ploidy); free(ma->q2p); free(ma->pdg);
		free(ma->phi); free(ma->phi_indel); free(ma->phi1); free(ma->phi2);
		free(ma->z); free(ma->zswap); free(ma->z1); free(ma->z2);
		free(ma->afs); free(ma->afs1);
		free(ma);
	}
}

static int cal_pdg(const bcf1_t *b, bcf_p1aux_t *ma)
{
	int i, j;
	long *p, tmp;
	p = alloca(b->n_alleles * sizeof(long));
	memset(p, 0, sizeof(long) * b->n_alleles);
	for (j = 0; j < ma->n; ++j) {
		const uint8_t *pi = ma->PL + j * ma->PL_len;
		double *pdg = ma->pdg + j * 3;
		pdg[0] = ma->q2p[pi[2]]; pdg[1] = ma->q2p[pi[1]]; pdg[2] = ma->q2p[pi[0]];
		for (i = 0; i < b->n_alleles; ++i)
			p[i] += (int)pi[(i+1)*(i+2)/2-1];
	}
	for (i = 0; i < b->n_alleles; ++i) p[i] = p[i]<<4 | i;
	for (i = 1; i < b->n_alleles; ++i) // insertion sort
		for (j = i; j > 0 && p[j] < p[j-1]; --j)
			tmp = p[j], p[j] = p[j-1], p[j-1] = tmp;
	for (i = b->n_alleles - 1; i >= 0; --i)
		if ((p[i]&0xf) == 0) break;
	return i;
}
// f0 is the reference allele frequency
static double mc_freq_iter(double f0, const bcf_p1aux_t *ma)
{
	double f, f3[3];
	int i;
	f3[0] = (1.-f0)*(1.-f0); f3[1] = 2.*f0*(1.-f0); f3[2] = f0*f0;
	for (i = 0, f = 0.; i < ma->n; ++i) {
		double *pdg;
		pdg = ma->pdg + i * 3;
		f += (pdg[1] * f3[1] + 2. * pdg[2] * f3[2])
			/ (pdg[0] * f3[0] + pdg[1] * f3[1] + pdg[2] * f3[2]);
	}
	f /= ma->n * 2.;
	return f;
}

int bcf_p1_call_gt(const bcf_p1aux_t *ma, double f0, int k)
{
	double sum, g[3];
	double max, f3[3], *pdg = ma->pdg + k * 3;
	int q, i, max_i, ploidy;
	ploidy = ma->ploidy? ma->ploidy[k] : 2;
	if (ploidy == 2) {
		f3[0] = (1.-f0)*(1.-f0); f3[1] = 2.*f0*(1.-f0); f3[2] = f0*f0;
	} else {
		f3[0] = 1. - f0; f3[1] = 0; f3[2] = f0;
	}
	for (i = 0, sum = 0.; i < 3; ++i)
		sum += (g[i] = pdg[i] * f3[i]);
	for (i = 0, max = -1., max_i = 0; i <= ploidy; ++i) {
		g[i] /= sum;
		if (g[i] > max) max = g[i], max_i = i;
	}
	max = 1. - max;
	if (max < 1e-308) max = 1e-308;
	q = (int)(-4.343 * log(max) + .499);
	if (q > 99) q = 99;
	return q<<2|max_i;
}

#define TINY 1e-20

static void mc_cal_y_core(bcf_p1aux_t *ma, int beg)
{
	double *z[2], *tmp, *pdg;
	int _j, last_min, last_max;
	assert(beg == 0 || ma->M == ma->n*2);
	z[0] = ma->z;
	z[1] = ma->zswap;
	pdg = ma->pdg;
	memset(z[0], 0, sizeof(double) * (ma->M + 1));
	memset(z[1], 0, sizeof(double) * (ma->M + 1));
	z[0][0] = 1.;
	last_min = last_max = 0;
	ma->t = 0.;
	if (ma->M == ma->n * 2) {
		for (_j = beg; _j < ma->n; ++_j) {
			int k, j = _j - beg, _min = last_min, _max = last_max;
			double p[3], sum;
			pdg = ma->pdg + _j * 3;
			p[0] = pdg[0]; p[1] = 2. * pdg[1]; p[2] = pdg[2];
			for (; _min < _max && z[0][_min] < TINY; ++_min) z[0][_min] = z[1][_min] = 0.;
			for (; _max > _min && z[0][_max] < TINY; --_max) z[0][_max] = z[1][_max] = 0.;
			_max += 2;
			if (_min == 0) 
				k = 0, z[1][k] = (2*j+2-k)*(2*j-k+1) * p[0] * z[0][k];
			if (_min <= 1)
				k = 1, z[1][k] = (2*j+2-k)*(2*j-k+1) * p[0] * z[0][k] + k*(2*j+2-k) * p[1] * z[0][k-1];
			for (k = _min < 2? 2 : _min; k <= _max; ++k)
				z[1][k] = (2*j+2-k)*(2*j-k+1) * p[0] * z[0][k]
					+ k*(2*j+2-k) * p[1] * z[0][k-1]
					+ k*(k-1)* p[2] * z[0][k-2];
			for (k = _min, sum = 0.; k <= _max; ++k) sum += z[1][k];
			ma->t += log(sum / ((2. * j + 2) * (2. * j + 1)));
			for (k = _min; k <= _max; ++k) z[1][k] /= sum;
			if (_min >= 1) z[1][_min-1] = 0.;
			if (_min >= 2) z[1][_min-2] = 0.;
			if (j < ma->n - 1) z[1][_max+1] = z[1][_max+2] = 0.;
			if (_j == ma->n1 - 1) { // set pop1; ma->n1==-1 when unset
				ma->t1 = ma->t;
				memcpy(ma->z1, z[1], sizeof(double) * (ma->n1 * 2 + 1));
			}
			tmp = z[0]; z[0] = z[1]; z[1] = tmp;
			last_min = _min; last_max = _max;
		}
	} else { // this block is very similar to the block above; these two might be merged in future
		int j, M = 0;
		for (j = 0; j < ma->n; ++j) {
			int k, M0, _min = last_min, _max = last_max;
			double p[3], sum;
			pdg = ma->pdg + j * 3;
			for (; _min < _max && z[0][_min] < TINY; ++_min) z[0][_min] = z[1][_min] = 0.;
			for (; _max > _min && z[0][_max] < TINY; --_max) z[0][_max] = z[1][_max] = 0.;
			M0 = M;
			M += ma->ploidy[j];
			if (ma->ploidy[j] == 1) {
				p[0] = pdg[0]; p[1] = pdg[2];
				_max++;
				if (_min == 0) k = 0, z[1][k] = (M0+1-k) * p[0] * z[0][k];
				for (k = _min < 1? 1 : _min; k <= _max; ++k)
					z[1][k] = (M0+1-k) * p[0] * z[0][k] + k * p[1] * z[0][k-1];
				for (k = _min, sum = 0.; k <= _max; ++k) sum += z[1][k];
				ma->t += log(sum / M);
				for (k = _min; k <= _max; ++k) z[1][k] /= sum;
				if (_min >= 1) z[1][_min-1] = 0.;
				if (j < ma->n - 1) z[1][_max+1] = 0.;
			} else if (ma->ploidy[j] == 2) {
				p[0] = pdg[0]; p[1] = 2 * pdg[1]; p[2] = pdg[2];
				_max += 2;
				if (_min == 0) k = 0, z[1][k] = (M0-k+1) * (M0-k+2) * p[0] * z[0][k];
				if (_min <= 1) k = 1, z[1][k] = (M0-k+1) * (M0-k+2) * p[0] * z[0][k] + k*(M0-k+2) * p[1] * z[0][k-1];
				for (k = _min < 2? 2 : _min; k <= _max; ++k)
					z[1][k] = (M0-k+1)*(M0-k+2) * p[0] * z[0][k] + k*(M0-k+2) * p[1] * z[0][k-1] + k*(k-1)* p[2] * z[0][k-2];
				for (k = _min, sum = 0.; k <= _max; ++k) sum += z[1][k];
				ma->t += log(sum / (M * (M - 1.)));
				for (k = _min; k <= _max; ++k) z[1][k] /= sum;
				if (_min >= 1) z[1][_min-1] = 0.;
				if (_min >= 2) z[1][_min-2] = 0.;
				if (j < ma->n - 1) z[1][_max+1] = z[1][_max+2] = 0.;
			}
			tmp = z[0]; z[0] = z[1]; z[1] = tmp;
			last_min = _min; last_max = _max;
		}
	}
	if (z[0] != ma->z) memcpy(ma->z, z[0], sizeof(double) * (ma->M + 1));
}

static void mc_cal_y(bcf_p1aux_t *ma)
{
	if (ma->n1 > 0 && ma->n1 < ma->n && ma->M == ma->n * 2) { // NB: ma->n1 is ineffective when there are haploid samples
		int k;
		long double x;
		memset(ma->z1, 0, sizeof(double) * (2 * ma->n1 + 1));
		memset(ma->z2, 0, sizeof(double) * (2 * (ma->n - ma->n1) + 1));
		ma->t1 = ma->t2 = 0.;
		mc_cal_y_core(ma, ma->n1);
		ma->t2 = ma->t;
		memcpy(ma->z2, ma->z, sizeof(double) * (2 * (ma->n - ma->n1) + 1));
		mc_cal_y_core(ma, 0);
		// rescale z
		x = expl(ma->t - (ma->t1 + ma->t2));
		for (k = 0; k <= ma->M; ++k) ma->z[k] *= x;
	} else mc_cal_y_core(ma, 0);
}

static void contrast(bcf_p1aux_t *ma, double pc[4]) // mc_cal_y() must be called before hand
{
	int k, n1 = ma->n1, n2 = ma->n - ma->n1;
	long double sum1, sum2;
	pc[0] = pc[1] = pc[2] = pc[3] = -1.;
	if (n1 <= 0 || n2 <= 0) return;
	for (k = 0, sum1 = 0.; k <= 2*n1; ++k) sum1 += ma->phi1[k] * ma->z1[k];
	for (k = 0, sum2 = 0.; k <= 2*n2; ++k) sum2 += ma->phi2[k] * ma->z2[k];
	pc[2] = ma->phi1[2*n1] * ma->z1[2*n1] / sum1;
	pc[3] = ma->phi2[2*n2] * ma->z2[2*n2] / sum2;
	for (k = 2; k < 4; ++k) {
		pc[k] = pc[k] > .5? -(-4.343 * log(1. - pc[k] + TINY) + .499) : -4.343 * log(pc[k] + TINY) + .499;
		pc[k] = (int)pc[k];
		if (pc[k] > 99) pc[k] = 99;
		if (pc[k] < -99) pc[k] = -99;
	}
	pc[0] = ma->phi2[2*n2] * ma->z2[2*n2] / sum2 * (1. - ma->phi1[2*n1] * ma->z1[2*n1] / sum1);
	pc[1] = ma->phi1[2*n1] * ma->z1[2*n1] / sum1 * (1. - ma->phi2[2*n2] * ma->z2[2*n2] / sum2);
	pc[0] = pc[0] == 1.? 99 : (int)(-4.343 * log(1. - pc[0]) + .499);
	pc[1] = pc[1] == 1.? 99 : (int)(-4.343 * log(1. - pc[1]) + .499);
}

static double mc_cal_afs(bcf_p1aux_t *ma, double *p_ref_folded, double *p_var_folded)
{
	int k;
	long double sum = 0., sum2;
	double *phi = ma->is_indel? ma->phi_indel : ma->phi;
	memset(ma->afs1, 0, sizeof(double) * (ma->M + 1));
	mc_cal_y(ma);
	// compute AFS
	for (k = 0, sum = 0.; k <= ma->M; ++k)
		sum += (long double)phi[k] * ma->z[k];
	for (k = 0; k <= ma->M; ++k) {
		ma->afs1[k] = phi[k] * ma->z[k] / sum;
		if (isnan(ma->afs1[k]) || isinf(ma->afs1[k])) return -1.;
	}
	// compute folded variant probability
	for (k = 0, sum = 0.; k <= ma->M; ++k)
		sum += (long double)(phi[k] + phi[ma->M - k]) / 2. * ma->z[k];
	for (k = 1, sum2 = 0.; k < ma->M; ++k)
		sum2 += (long double)(phi[k] + phi[ma->M - k]) / 2. * ma->z[k];
	*p_var_folded = sum2 / sum;
	*p_ref_folded = (phi[k] + phi[ma->M - k]) / 2. * (ma->z[ma->M] + ma->z[0]) / sum;
	// the expected frequency
	for (k = 0, sum = 0.; k <= ma->M; ++k) {
		ma->afs[k] += ma->afs1[k];
		sum += k * ma->afs1[k];
	}
	return sum / ma->M;
}

int bcf_p1_cal(const bcf1_t *b, bcf_p1aux_t *ma, bcf_p1rst_t *rst)
{
	int i, k;
	long double sum = 0.;
	ma->is_indel = bcf_is_indel(b);
	// set PL and PL_len
	for (i = 0; i < b->n_gi; ++i) {
		if (b->gi[i].fmt == bcf_str2int("PL", 2)) {
			ma->PL = (uint8_t*)b->gi[i].data;
			ma->PL_len = b->gi[i].len;
			break;
		}
	}
	if (b->n_alleles < 2) return -1; // FIXME: find a better solution
	// 
	rst->rank0 = cal_pdg(b, ma);
	rst->f_exp = mc_cal_afs(ma, &rst->p_ref_folded, &rst->p_var_folded);
	rst->p_ref = ma->afs1[ma->M];
	for (k = 0, sum = 0.; k < ma->M; ++k)
		sum += ma->afs1[k];
	rst->p_var = (double)sum;
	// calculate f_flat and f_em
	for (k = 0, sum = 0.; k <= ma->M; ++k)
		sum += (long double)ma->z[k];
	rst->f_flat = 0.;
	for (k = 0; k <= ma->M; ++k) {
		double p = ma->z[k] / sum;
		rst->f_flat += k * p;
	}
	rst->f_flat /= ma->M;
	{ // calculate f_em
		double flast = rst->f_flat;
		for (i = 0; i < MC_MAX_EM_ITER; ++i) {
			rst->f_em = mc_freq_iter(flast, ma);
			if (fabs(rst->f_em - flast) < MC_EM_EPS) break;
			flast = rst->f_em;
		}
	}
	{ // estimate equal-tail credible interval (95% level)
		int l, h;
		double p;
		for (i = 0, p = 0.; i < ma->M; ++i)
			if (p + ma->afs1[i] > 0.025) break;
			else p += ma->afs1[i];
		l = i;
		for (i = ma->M-1, p = 0.; i >= 0; --i)
			if (p + ma->afs1[i] > 0.025) break;
			else p += ma->afs1[i];
		h = i;
		rst->cil = (double)(ma->M - h) / ma->M; rst->cih = (double)(ma->M - l) / ma->M;
	}
	rst->g[0] = rst->g[1] = rst->g[2] = -1.;
	contrast(ma, rst->pc);
	return 0;
}

void bcf_p1_dump_afs(bcf_p1aux_t *ma)
{
	int k;
	fprintf(stderr, "[afs]");
	for (k = 0; k <= ma->M; ++k)
		fprintf(stderr, " %d:%.3lf", k, ma->afs[ma->M - k]);
	fprintf(stderr, "\n");
	memset(ma->afs, 0, sizeof(double) * (ma->M + 1));
}
