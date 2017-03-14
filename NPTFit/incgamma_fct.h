#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_errno.h>

void inc_gamma_error_msg(const char* fn_name, const char* reason, const char* file, int line, int gsl_errno) {
	fprintf(stderr, "Error in GSL %s\n", fn_name);
	fprintf(stderr, "GSL Error %i: %s.\n", gsl_errno, gsl_strerror(gsl_errno));
	fprintf(stderr, "%s in %s:%i.\n", reason, file, line);
}

void gamma_inc_lower_error_handler(const char* reason, const char* file, int line, int gsl_errno) {
	if (gsl_errno != GSL_EUNDRFLW) {
		inc_gamma_error_msg("gamma_inc_P while computing lower incomplete gamma ratio.", reason, file, line, gsl_errno);
		abort();
	}
}

void gamma_inc_upper_error_handler(const char* reason, const char* file, int line, int gsl_errno) {
	if (gsl_errno != GSL_EUNDRFLW) {
		inc_gamma_error_msg("gamma_inc while computing upper incomplete gamma ratio.", reason, file, line, gsl_errno);
		abort();
	}
}
	
#define INC_GAMMA_NO_UNDERFLOW(fn, handler) do { \
		gsl_error_handler_t* old_handler = gsl_set_error_handler(& handler ); \
		fn; \
		gsl_set_error_handler(old_handler); \
	} while (0)

void inc_gamma_check_value(long double val, const char* fn_name) {
	if (isnan(val) || !isfinite(val)) {
		fprintf(stderr, "Error while computing %s.\n", fn_name);
		if (isnan(val)) {
			fprintf(stderr, "NaN detected. ");
		} else {
			fprintf(stderr, "Infinity detected. ");
		}
		fprintf(stderr, "This is likely the result of an overflow.\n");
		abort();
	}
}

/*
	Computes the upper Gamma(x+y,a)/Gamma(x) for x = [1 ... x_max]
*/
void inc_gamma_upper_ratio_ary(int x_max, double y, double a, double* result) {
	
	// x = 1 case
	size_t idx = 0;
	int x = 1;
	long double ln_gamma_xp1 = 0;
	// With fn(x) = Gamma(x+y, a)/Gamma(x)
	long double fn_x;
	INC_GAMMA_NO_UNDERFLOW(fn_x = gsl_sf_gamma_inc(x+y, a), gamma_inc_upper_error_handler);
	result[idx] = fn_x;
	
	for (;x < x_max;) {
		ln_gamma_xp1 += logl((long double) x);
		// Using the recurrence relation Gamma(s+1,a) = s Gamma(s,a) + a^s exp(-a)
		// We have Gamma((x+1)+y)/Gamma(x+1) = (x+y) ( Gamma(x+y)/Gamma(x) )/x + a^(x+y) exp(-a)/Gamma(x+1)
		// -> fn(x+1) = (x+y) fn(x)/x + exp((x+y) ln a - a - ln Gamma(x+1))
		long double fn_xp1 = (x+y) * fn_x/x + expl((x+y)*logl(a) - a - ln_gamma_xp1);
		result[idx+1] = (double)fn_xp1;
		fn_x = fn_xp1;
		idx += 1;
		x += 1;
	}
	// Make use of the fact that NaNs and infinities will propagate to the last value computed.
	inc_gamma_check_value(result[x_max-1], "upper incomplete gamma ratio");
}

/*
	Computes the lower gamma(x+y,a)/Gamma(x) for x = [1 ... x_max]
*/
void inc_gamma_lower_ratio_ary(int x_max, double y, double a, double* result) {
	
	// x = x_max case
	size_t idx = x_max - 1;
	int x = x_max;
	long double ln_gamma_x = gsl_sf_lngamma(x);
	// With fn(x) = gamma(x+y, a)/Gamma(x)
	long double fn_x;
	//INC_GAMMA_NO_UNDERFLOW(fn_x = expl((x+y)*logl(a) - a - logl(x+y) - ln_gamma_x)*gsl_sf_hyperg_1F1(1, x+y+1, a), hyperg_error_handler);
	INC_GAMMA_NO_UNDERFLOW(fn_x = gsl_sf_gamma_inc_P(x+y,a) * expl(gsl_sf_lngamma(x+y) - ln_gamma_x), gamma_inc_lower_error_handler);
	result[idx] = fn_x;
	
	for (;x>1;) {
		x -= 1;
		idx -= 1;
		ln_gamma_x -= logl((long double) x);
		long double fn_xp1 = fn_x;
		// We have gamma(s,a) = (gamma(s+1, a) + a^s exp(-a))/s
		// So gamma(x+y,a)/Gamma(x) = ( x gamma((x+1)+y, a)/Gamma(x+1) + a^(x+y) exp(-a)/Gamma(x) ) / s
		// -> fn(x) = ( x fn(x+1) + exp( (x+y) ln a - a - ln Gamma(x)! ) ) / (x+y)
		fn_x = ( x*fn_xp1 + expl( (x+y)*logl(a) - a - ln_gamma_x ) ) / (x+y);
		result[idx] = (double)fn_x;
	}
	// Make use of the fact that NaNs and infinities will propagate to the last value computed.
	inc_gamma_check_value(result[0], "lower incomplete gamma ratio");
}

