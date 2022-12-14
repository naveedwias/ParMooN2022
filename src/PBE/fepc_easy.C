#include "fepc_easy.h"
#include "basictools.h"

fepc_real_t legendre(int degree, fepc_real_t x) {
    if (degree == 0) {
        return 1;
    } else if (degree == 1) {
        return x;
    } else {
        return ((2*degree-1)*x*legendre(degree-1, x) - (degree-1)*legendre(degree-2, x)) / degree;
    }
}
 
fepc_real_t phi_l(int& step, int& v, int p, fepc_real_t& x, fepc_real_t& stepping) {
    fepc_real_t h_l, result;
    
    h_l = get_h_l(step, stepping);
    result = 1.;
    
	if (v*h_l > x || (v+1)*h_l < x) {
            return 0;
        }
    result *= std::sqrt((2*p + 1)/h_l)*legendre(p, 2.0*x/h_l - 2.0*v - 1.0); 

    return result;
}

fepc_real_t get_value_at_step(func_p0_p function, fepc_real_t& x, int& step, fepc_real_t& stepping, int* count) {
    ASSERT(function != nullptr && step <= function->maxlevel);
    
    int n, r;
	int lang = function->hierarchie[step]->vektor->lang;
    
    fepc_real_t result = 0.;
    
	if(step == function->maxlevel) {
		for (n = 0; n < count[step]; n++) {
			r = n % lang;

			result += function->hierarchie[step]->vektor->glied[r] * phi_l(step, r, 0, x, stepping);
		}
	} else {
		for (n = 0; n < count[step]; n++) {
			r = n % lang;

			if (!is_in_latter_interval(r, function->hierarchie[step+1]->vektor)) {
				result += function->hierarchie[step]->vektor->glied[r] * phi_l(step, r, 0, x, stepping);
			} 
		}

	}

    return result;
}

void set_gridstructure(func_p0_p function, interval_p* intervals, fepc_real_t stepping) {
    int start, lang;
   
    int n, steps;

    steps = function->maxlevel+1;
    
    fepc_real_t pow, h_l;
    
    pow = 2.;    
    
    for (n = 0; n < steps; n++) {
        pow *= 0.5;
        h_l = pow*stepping;
            
		if (intervals[n] != nullptr) {
            start = iround(intervals[n]->start / h_l);
            lang = iround(intervals[n]->end / h_l) - start;  // lang = length (integer), see p. 68; round() is needed to correct error
		}

        folge_del(function->hierarchie[n]->vektor);
		function->hierarchie[n]->vektor = folge_new(lang);
    }
}

fepc_real_t get_value_at_step(func_p1_p function, fepc_real_t& x, int& step, fepc_real_t& stepping, int* count) {
    ASSERT(function != nullptr && step <= function->maxlevel);
    
    int n, r;
	int lang = function->hierarchie[step]->vektor0->lang;
    
    fepc_real_t result = 0.;
    
	if(step == function->maxlevel) {
		for (n = 0; n < count[step]; n++) {
			r = n % lang;

			result += function->hierarchie[step]->vektor0->glied[r] * phi_l(step, r, 0, x, stepping);
			result += function->hierarchie[step]->vektor1->glied[r] * phi_l(step, r, 1, x, stepping);
		}
	} else {
		for (n = 0; n < count[step]; n++) {
			r = n % lang;

			if (!is_in_latter_interval(r, function->hierarchie[step+1]->vektor0)) {
				result += function->hierarchie[step]->vektor0->glied[r] * phi_l(step, r, 0, x, stepping);
				result += function->hierarchie[step]->vektor1->glied[r] * phi_l(step, r, 1, x, stepping);
			} 
		}

	}

    return result;
}

void set_gridstructure(func_p1_p function, interval_p* intervals, fepc_real_t stepping) {
    int start, lang;
   
    int n, steps;

    steps = function->maxlevel+1;
    
    fepc_real_t pow, h_l;
    
    pow = 2.;    
    
    for (n = 0; n < steps; n++) {
        pow *= 0.5;
        h_l = pow*stepping;
        
		if (intervals[n] != nullptr) {
            start = iround(intervals[n]->start / h_l);
            lang = iround(intervals[n]->end / h_l) - start;  // lang = length (integer), see p. 68; round() is needed to correct error
		}

        folge_p1_del(function->hierarchie[n]->vektor0);
		function->hierarchie[n]->vektor0 = folge_p1_new(lang);

        folge_p1_del(function->hierarchie[n]->vektor1);
		function->hierarchie[n]->vektor1 = folge_p1_new(lang);
    }
}