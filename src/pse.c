/*
 * National Center for Supercomputing Applications
 * University of Illinois at Urbana-Champaign
 *
 * Large-Scale Agent-Based Social Simulation
 * Les Gasser, NCSA Fellow
 *
 * Author: Santiago Nunez-Corrales
 */
#include <ranlib.h>
#include <rnglib.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <pse.h>

/*
 * Declaration of private functions
 *
 * Notes:
 * ------
 * 1. We define separate versions of randomization for int and unsigned int.
 *    This is useful for dealing separately with data structures. Also, we
 *    use separate definitions of randomization for integer values within some
 *    particular limit.
 * 2. Additionaly, we provide randomization alternatives for compute an update
 * 	  or compute and discard.
 */
int pse_sizeof(pse_variable *);
unsigned int pse_is_world_var(pse_variable *);
int pse_sample_int_distribution(int, double *, pse_distribution_type);
unsigned int pse_is_int_distribution(pse_distribution_type);
double pse_sample_double_distribution(double, double *, pse_distribution_type);
void pse_randomize(pse_variable *, pse_variable *, unsigned int location);
void pse_randomize_and_alter(pse_variable *, pse_variable *, unsigned int location, pse_error *error);

/*
 * Calculate the size of registered content
 */
int pse_sizeof(pse_variable *var) {
	if (var->array == PSE_SCALAR) {
		switch(var->storage) {
		case PSE_VAR_INT:
			return sizeof(int);
		case PSE_VAR_DOUBLE:
			return sizeof(double);
		case PSE_VAR_STRING:
			return PSE_MAX_STRLEN;
		case PSE_VAR_TIME:
			return sizeof(pse_time);
		default:
			return 0;
		}
	} else if (var->array == PSE_ARRAY) {
		switch(var->storage) {
		case PSE_VAR_INT:
			return sizeof(int)*var->size;
		case PSE_VAR_DOUBLE:
			return sizeof(double)*var->size;
		case PSE_VAR_STRING:
			return PSE_MAX_STRLEN*var->size;
		case PSE_VAR_TIME:
			return sizeof(pse_time)*var->size;
		default:
			return 0;
		}
	} else {
		return 0;
	}
}

unsigned int pse_is_world_var(pse_variable *var) {
	if (var->locality == PSE_WORLD)
		return PSE_TRUE;
	else
		return PSE_FALSE;
}

/*
 * Obtain a random number from one amongst many integer distributions
 */
int pse_sample_int_distribution(int value, double *pars,
										pse_distribution_type distribution) {
	int min;
	int max;
	double p;
	double mu;

	switch(distribution) {
	case PSE_DIST_UNIFORM_INT_SELF:
		min = 0;
		return ignuin(min, value);
		break;
	case PSE_DIST_UNIFORM_INT_BOUNDED:
		min = round(pars[0]);
		max = round(pars[1]);
		return ignuin(min, max);
	case PSE_DIST_BERNOULLI:
		p = pars[0];
		return genunf(0,1) > p ? PSE_HEADS : PSE_TAILS;
	case PSE_DIST_BINOMIAL:
		max = round(pars[0]);
		p = pars[1];
		return ignbin(max,p);
	case PSE_DIST_BINOMIAL_SELF:
		p = pars[0];
		max = value;
		return ignbin(max,p);
	case PSE_DIST_NEG_BINOMIAL:
		p = pars[0];
		max = round(pars[1]);
		return ignnbn(max,p);
	case PSE_DIST_NEG_BINOMIAL_SELF:
		p = pars[0];
		max = value;
		return ignnbn(max,p);
	case PSE_DIST_POISSON:
		mu = pars[0];
		return ignpoi(mu);
	case PSE_DIST_POISSON_SELF:
		mu = value;
		return ignpoi(mu);
	case PSE_DIST_NONE:
		return value;
	default:
		return 0;
	}
}

/*
 * Obtain a random number from one amongst many integer distributions
 */
double pse_sample_double_distribution(double value, double *pars,
										pse_distribution_type distribution) {
	/*
	 * Parameter names are preserved for clarity
	 */
	double min;
	double max;
	double alpha;
	double beta;
	double mu;
	double sigma;
	double dfn;
	double dfd;
	double df;

	switch(distribution) {
	case PSE_DIST_UNIFORM_DOUBLE_SELF:
		min = 0;
		return genunf(min, value);
		break;
	case PSE_DIST_UNIFORM_DOUBLE_BOUNDED:
		min = pars[0];
		max = pars[1];
		return genunf(min, max);
	case PSE_DIST_NORMAL:
		mu = pars[0];
		sigma = pars[1];
		return gennor(mu,sigma);
	case PSE_DIST_NORMAL_SELF:
		mu = value;
		sigma = pars[0];
		return gennor(mu,sigma);
	case PSE_DIST_EXPONENTIAL:
		mu = pars[0];
		return genexp(mu);
	case PSE_DIST_EXPONENTIAL_SELF:
		mu = value;
		return genexp(mu);
	case PSE_DIST_GAMMA:
		/*
		 * alpha: shape constant
		 * beta: rate constant
		 *
		 * Note: parameter order in ranlib is counter-intuitive
		 */
		alpha = pars[1];
		beta = pars[0];
		return gengam(beta, alpha);
	case PSE_DIST_GAMMA_SELF:
		alpha = pars[1];
		beta = value;
		return gengam(beta, alpha);
	case PSE_DIST_F:
		/*
		 * F statistics are independent of value. They represent a proportion
		 * of the ratio of variations between sample and population variance
		 * for two populations.
		 */
		dfn = pars[0];
		dfd = pars[1];
		return genf(dfn, dfd);
	case PSE_DIST_BETA:
		alpha = pars[0];
		beta = pars[1];
		return genbet(alpha, beta);
	case PSE_DIST_CHISQ:
		df = pars[0];
		return genchi(df);
	case PSE_DIST_CHISQ_SELF:
		df = value;
		return genchi(df);
	case PSE_DIST_NONE:
		return value;
	case PSE_DIST_FOKKER_PLANCK:
		/*
		 * Placeholder for future implemenentation
		 */
		return value;
	case PSE_DIST_CUSTOM:
		/*
		 * Placeholder for future implementation
		 */
		return value;
	default:
		return 0;
	}
}

unsigned int pse_is_int_distribution(pse_distribution_type distribution) {
	switch(distribution) {
	case PSE_DIST_UNIFORM_INT_SELF:
	case PSE_DIST_UNIFORM_INT_BOUNDED:
	case PSE_DIST_BERNOULLI:
	case PSE_DIST_BINOMIAL:
	case PSE_DIST_BINOMIAL_SELF:
	case PSE_DIST_NEG_BINOMIAL:
	case PSE_DIST_NEG_BINOMIAL_SELF:
	case PSE_DIST_POISSON:
	case PSE_DIST_POISSON_SELF:
		return PSE_TRUE;
	default:
		return PSE_FALSE;
	}
}

unsigned int pse_is_self_distribution(pse_distribution_type distribution) {
	switch(distribution) {
	case PSE_DIST_UNIFORM_INT_SELF:
	case PSE_DIST_BINOMIAL_SELF:
	case PSE_DIST_NEG_BINOMIAL_SELF:
	case PSE_DIST_POISSON_SELF:
	case PSE_DIST_UNIFORM_DOUBLE_SELF:
	case PSE_DIST_NORMAL_SELF:
	case PSE_DIST_EXPONENTIAL_SELF:
	case PSE_DIST_GAMMA_SELF:
	case PSE_DIST_CHISQ_SELF:
		return PSE_TRUE;
	default:
		return PSE_FALSE;
	}
}

/*
 * Randomize provides stochasticity into agent models.
 *
 * In the case of strings, two steps are required:
 * 1. Find a randomization value corresponding to a valid location in the array.
 * 2. Find a randomization value adequate for ASCII text. For the moment, we do
 *    not concern ourselves with unicode.
 *
 * In any case, we assume that the probability distributions have the adequate
 * parameters to generate the values. The responsibility is in the hands of model
 * developers to understand the statistics behind ay phenomenology being portrayed.
 *
 * TODO: this method is ugly for strings. A common factorization should be possible.
 */
void pse_randomize(pse_variable *ptr_out, pse_variable *var, unsigned int location) {
	unsigned int array_location;
	unsigned int str_len;
	char char_median = (char)127;
	char char_max = (char)255;

	if (var->array == PSE_SCALAR) {
		switch(var->storage) {
		case PSE_VAR_INT:
			ptr_out->content.cint = pse_sample_int_distribution(var->content.cint,
										var->point_parameters, var->point_distribution);
			break;
		case PSE_VAR_DOUBLE:
			ptr_out->content.cdouble = pse_sample_double_distribution(var->content.cdouble,
										var->point_parameters, var->point_distribution);
			break;
		case PSE_VAR_STRING:
			/*
			 * Copy the contents and the parameters
			 */
			strcpy(ptr_out->content.cstring, var->content.cstring);
			memcpy(ptr_out->array_parameters, ptr_out->array_parameters,
										sizeof(double)*PSE_MAX_DIST_PARAMS);
			str_len = strlen(ptr_out->content.cstring);
			/*
			 * First, we start by determining the location to be altered in the string
			 * based on its size and parameters.
			 */
			if (pse_is_int_distribution(var->array_distribution) == PSE_TRUE) {
				/*
				 * Handle integer distributions
				 */
				if (pse_is_self_distribution(var->array_distribution) == PSE_TRUE) {
					do {
						array_location = pse_sample_int_distribution(str_len/2,
								ptr_out->array_parameters, var->array_distribution);
					} while(array_location >= str_len);
				} else {
					ptr_out->array_parameters[0] = str_len/2;

					do {
						array_location = pse_sample_int_distribution(str_len/2,
								ptr_out->array_parameters, var->array_distribution);
					} while(array_location >= str_len);
				}
			} else {
				/*
				 * Handle floating point distributions
				 */
				if (pse_is_self_distribution(var->array_distribution) == PSE_TRUE) {
					do {
						array_location = (unsigned int)pse_sample_double_distribution(str_len/2,
								ptr_out->array_parameters, var->array_distribution);
					} while(array_location >= str_len);
				} else {
					ptr_out->array_parameters[0] = str_len/2;

					do {
						array_location = (unsigned int)pse_sample_int_distribution(str_len/2,
								ptr_out->array_parameters, var->array_distribution);
					} while(array_location >= str_len);
				}
			}

			/*
			 * Once we have a proper location, alter the char data.
			 */
			if (pse_is_int_distribution(var->point_distribution) == PSE_TRUE) {
				/*
				 * Handle integer distributions
				 */
				if (pse_is_self_distribution(var->point_distribution) == PSE_TRUE) {
					do {
						ptr_out->content.cstring[array_location] =
								(char)pse_sample_int_distribution(char_median,
								ptr_out->point_parameters, var->point_distribution);
					} while(ptr_out->content.cstring[array_location] > char_max);
				} else {
					ptr_out->array_parameters[0] = str_len/2;

					do {
						ptr_out->content.cstring[array_location] =
								(char)pse_sample_int_distribution(char_median,
								ptr_out->point_parameters, var->point_distribution);
					} while(ptr_out->content.cstring[array_location] > char_max);
				}
			} else {
				/*
				 * Handle floating point distributions
				 */
				if (pse_is_self_distribution(var->array_distribution) == PSE_TRUE) {
					do {
						ptr_out->content.cstring[array_location] =
								(char)pse_sample_double_distribution(char_median,
								ptr_out->point_parameters, var->point_distribution);
					} while(ptr_out->content.cstring[array_location] > char_max);
				} else {
					ptr_out->array_parameters[0] = char_median;

					do {
						ptr_out->content.cstring[array_location] =
								(char)pse_sample_double_distribution(char_median,
								ptr_out->point_parameters, var->point_distribution);
					} while(ptr_out->content.cstring[array_location] > char_max);
				}
			}

			break;
		case PSE_VAR_TIME:
			ptr_out->content.ctime = pse_sample_double_distribution(var->content.ctime,
										var->point_parameters, var->point_distribution);
			break;
		default:
			break;
		}
	} else {
		switch(var->storage) {
		case PSE_VAR_INT:
			ptr_out->content.cint_a[location] =
					pse_sample_int_distribution(var->content.cint_a[location],
										var->point_parameters, var->point_distribution);
			break;
		case PSE_VAR_DOUBLE:
			ptr_out->content.cdouble_a[location] =
					pse_sample_double_distribution(var->content.cdouble_a[location],
										var->point_parameters, var->point_distribution);
			break;
		case PSE_VAR_STRING:
			/*
			 * Copy the contents and the parameters
			 */
			strcpy(ptr_out->content.cstring_a[location], var->content.cstring_a[location]);
			memcpy(ptr_out->array_parameters, ptr_out->array_parameters,
										sizeof(double)*PSE_MAX_DIST_PARAMS);
			str_len = strlen(ptr_out->content.cstring_a[location]);
			/*
			 * First, we start by determining the location to be altered in the string
			 * based on its size and parameters.
			 */
			if (pse_is_int_distribution(var->array_distribution) == PSE_TRUE) {
				/*
				 * Handle integer distributions
				 */
				if (pse_is_self_distribution(var->array_distribution) == PSE_TRUE) {
					do {
						array_location = pse_sample_int_distribution(str_len/2,
								ptr_out->array_parameters, var->array_distribution);
					} while(array_location >= str_len);
				} else {
					ptr_out->array_parameters[0] = str_len/2;

					do {
						array_location = pse_sample_int_distribution(str_len/2,
								ptr_out->array_parameters, var->array_distribution);
					} while(array_location >= str_len);
				}
			} else {
				/*
				 * Handle floating point distributions
				 */
				if (pse_is_self_distribution(var->array_distribution) == PSE_TRUE) {
					do {
						array_location = (unsigned int)pse_sample_double_distribution(str_len/2,
								ptr_out->array_parameters, var->array_distribution);
					} while(array_location >= str_len);
				} else {
					ptr_out->array_parameters[0] = str_len/2;

					do {
						array_location = (unsigned int)pse_sample_int_distribution(str_len/2,
								ptr_out->array_parameters, var->array_distribution);
					} while(array_location >= str_len);
				}
			}

			/*
			 * Once we have a proper location, alter the char data.
			 */
			if (pse_is_int_distribution(var->point_distribution) == PSE_TRUE) {
				/*
				 * Handle integer distributions
				 */
				if (pse_is_self_distribution(var->point_distribution) == PSE_TRUE) {
					do {
						ptr_out->content.cstring_a[location][array_location] =
								(char)pse_sample_int_distribution(char_median,
								ptr_out->point_parameters, var->point_distribution);
					} while(ptr_out->content.cstring_a[location][array_location] > char_max);
				} else {
					ptr_out->array_parameters[0] = char_median;

					do {
						ptr_out->content.cstring[array_location] = (char)pse_sample_int_distribution(char_median,
								ptr_out->array_parameters, var->array_distribution);
					} while(ptr_out->content.cstring_a[location][array_location] > char_max);
				}
			} else {
				/*
				 * Handle floating point distributions
				 */
				if (pse_is_self_distribution(var->array_distribution) == PSE_TRUE) {
					do {
						ptr_out->content.cstring_a[location][array_location] =
								(char)pse_sample_double_distribution(str_len/2,
								ptr_out->array_parameters, var->array_distribution);
					} while(ptr_out->content.cstring[array_location] > char_max);
				} else {
					ptr_out->array_parameters[0] = str_len/2;

					do {
						ptr_out->content.cstring_a[location][array_location] =
								(char)pse_sample_double_distribution(char_median,
								ptr_out->array_parameters, var->array_distribution);
					} while(ptr_out->content.cstring_a[location][array_location] > char_max);
				}
			}
			break;
		case PSE_VAR_TIME:
			ptr_out->content.ctime_a[location] =
					pse_sample_double_distribution(var->content.ctime_a[location],
										var->point_parameters, var->point_distribution);
			break;
		default:
			break;
		}
	}

	return;
}

/*
 * Randomize and alter, used for replacing values and associated more closely
 * with SELF distributions.
 */
void pse_randomize_and_alter(pse_variable *ptr_out, pse_variable *var, unsigned int location, pse_error *error) {
	if (var->read_and_alter == PSE_FALSE) {
		*error = PSE_ERROR_VARIABLE_IS_IMMUTABLE;
		return;
	}

	pse_randomize(ptr_out, var, location);

	/*
	 * Update contents of the original variable
	 */
	var->content = ptr_out->content;
	memcpy(var->point_parameters, ptr_out->point_parameters, PSE_MAX_DIST_PARAMS*sizeof(double));
	memcpy(var->array_parameters, ptr_out->array_parameters, PSE_MAX_DIST_PARAMS*sizeof(double));

	return;
}

/*
 * PSE initialization.
 *
 * PSE stubs cannot be reused.
 */
pse_error pse_init(pse_agent_stub *pse) {
	int i;

	if (pse->state == INITIALIZED)
		return PSE_ERROR_ALREADY_INITIALIZED;

	if (pse->state == STARTED)
			return PSE_ERROR_ALREADY_STARTED;

	if (pse->state == FINALIZED)
			return PSE_ERROR_ALREADY_FINALIZED;

	for (i = 0; i < PSE_MAX_VARIABLES; i++) {
		pse->variables[i] = NULL;
		pse->dependencies[i] = NULL;
	}

	pse->var_count = 0;
	pse->var_limit = 0;
	pse->state = INITIALIZED;

	return PSE_ERROR_OK;
}

/*
 * PSE start
 *
 * The start state marks when basic initialization has occurred and variable
 * registration has occurred. This code is likely to vary in the future for
 * dynamic cognitive models. The apparent inefficiency in case handling
 * (explicitness) is required to remind future implementors of this. For the
 * current time, it is only a flag.
 *
 * When the machine is started, the random number generators are initialized
 * in each agent with a provided seed.
 */
pse_error pse_start(pse_agent_stub *pse, int seed_1, int seed_2) {
	if (pse->state == CREATED)
			return PSE_ERROR_NOT_INITIALIZED;

	if (pse->state == STARTED)
			return PSE_ERROR_ALREADY_STARTED;

	if (pse->state == FINALIZED)
			return PSE_ERROR_ALREADY_FINALIZED;

	/*
	 * Initialize the random number generators and set both seeds.
	 */
	initialize();
	set_seed(seed_1, seed_2);

	pse->state = STARTED;

	return PSE_ERROR_OK;
}

/*
 * PSE finalization
 */
pse_error pse_finalize(pse_agent_stub *pse) {
	int i;

	if (pse->state == FINALIZED)
		return PSE_ERROR_ALREADY_FINALIZED;

	if (pse->state == INITIALIZED)
		return PSE_ERROR_NOT_STARTED;

	for (i = 0; i < PSE_MAX_VARIABLES; i++) {
		if (pse->variables[i] != NULL) {
			if (pse->variables[i]->array == PSE_ARRAY) {
				switch(pse->variables[i]->storage) {
				case PSE_VAR_INT:
					free(pse->variables[i]->content.cint_a);
					break;
				case PSE_VAR_DOUBLE:
					free(pse->variables[i]->content.cdouble_a);
					break;
				case PSE_VAR_STRING:
					free(pse->variables[i]->content.cstring_a);
					break;
				case PSE_VAR_TIME:
					free(pse->variables[i]->content.ctime_a);
					break;
				default:
					return PSE_ERROR_TYPE_UNKNOWN;
				}
			}

			free(pse->variables[i]);
		}

		if (pse->dependencies[i] != NULL)
			free(pse->dependencies[i]);
	}

	pse->var_count = 0;
	pse->var_limit = 0;
	pse->state = FINALIZED;

	return PSE_ERROR_OK;
}

/*
 * PSE variable registration
 *
 * Variables are registered based on the last step prior to serving.
 */
pse_varid pse_register(pse_agent_stub *pse, pse_storage_type storage,
						pse_model_type model, pse_locality_type locality,
						pse_distribution_type point_distribution,
						double *point_parameters, pse_array_type array,
						unsigned int size, unsigned int read_and_alter,
						pse_distribution_type array_distribution,
						double *array_parameters, char *name) {
	int i;

	pse_varid next_available_varid = -1;
	pse_variable *p_to_var;

	if (pse->state == CREATED)
		return PSE_ERROR_NOT_INITIALIZED;

	if (pse->state == STARTED)
		return PSE_ERROR_ALREADY_STARTED;

	if (pse->state == FINALIZED)
		return PSE_ERROR_ALREADY_FINALIZED;

	if (pse->var_count == PSE_MAX_VARIABLES)
		return PSE_ERROR_TOO_MANY_VARIABLES;

	/*
	 * Important: variable count is not the basis for this value. Future
	 * careless use of it may lead to severe memory fragmentation in agents
	 * that change their internal representation. Conclusion: information
	 * processing classifiable as complex dissipates much more energy.
	 */
	next_available_varid = pse->var_limit;

	if (pse->variables[next_available_varid] != NULL)
		return PSE_ERROR_VARIABLE_ALREADY_REGISTERED;

	pse->variables[next_available_varid] = (pse_variable *) malloc(sizeof(pse_variable));
	p_to_var = pse->variables[next_available_varid];
	p_to_var->storage = storage;
	p_to_var->model = model;
	p_to_var->locality = locality;
	p_to_var->point_distribution = point_distribution;
	p_to_var->array = array;
	p_to_var->size = size;

	memcpy(p_to_var->point_parameters, point_parameters, PSE_MAX_DIST_PARAMS*sizeof(double));
	p_to_var->has_dependencies = PSE_FALSE;
	p_to_var->read_and_alter = read_and_alter;
	p_to_var->array_distribution = array_distribution;
	memcpy(p_to_var->array_parameters, array_parameters, PSE_MAX_DIST_PARAMS*sizeof(double));
	strcpy(p_to_var->name, name);


	/*
	 * We process registration based on content type
	 */
	if(p_to_var->array == PSE_ARRAY) {
		switch(p_to_var->storage) {
		case PSE_VAR_INT:
			p_to_var->content.cint_a = (int *)malloc(sizeof(int)*p_to_var->size);
			break;
		case PSE_VAR_DOUBLE:
			p_to_var->content.cdouble_a = (double *)malloc(sizeof(double)*p_to_var->size);
			break;
		case PSE_VAR_TIME:
			p_to_var->content.ctime_a = (pse_time *)malloc(sizeof(pse_time)*p_to_var->size);
			break;
		case PSE_VAR_STRING:
			/*
			 * Initialize the structure and then all the strings associated to it.
			 * This allows constructing models that contain a variable number of
			 * strings of at most length 1000.
			 */
			p_to_var->content.cstring_a = (char **)malloc(sizeof(char *)*p_to_var->size);
			for (i = 0; i < p_to_var->size; i++)
				p_to_var->content.cstring_a[i] = (char *)malloc(sizeof(char)*PSE_MAX_STRLEN);
			break;
		default:
			return PSE_ERROR_TYPE_UNKNOWN;
		}

		p_to_var->array_distribution = array_distribution;
	} else {
		/*
		 * Strings require initialization
		 */
		switch(p_to_var->storage) {
		case PSE_VAR_STRING:
			p_to_var->content.cstring = (char *)malloc(sizeof(char)*PSE_MAX_STRLEN);
			break;
		default:
			break;
		}

		p_to_var->array_distribution = PSE_DIST_NONE;
	}

	pse->var_count++;
	pse->var_limit++;

	return next_available_varid;
}

/*
 * PSE variable deregistration
 *
 * This function is reseved for future use, but for the moment mirrors the
 * role of register in terms of the underlying state machine.
 */
pse_error pse_deregister(pse_agent_stub *pse, pse_varid varid) {
	int i;

	if (pse->state == CREATED)
		return PSE_ERROR_NOT_INITIALIZED;

	if (pse->state == STARTED)
		return PSE_ERROR_ALREADY_STARTED;

	if (pse->state == FINALIZED)
		return PSE_ERROR_ALREADY_FINALIZED;

	if (pse->var_count == 0)
		return PSE_ERROR_TOO_MANY_VARIABLES;

	if (pse->variables[varid] == NULL)
		return PSE_ERROR_VARIABLE_UNKNOWN;

	/*
	 * Similarly, de-registration involves removing memory assignments
	 * of array structures.
	 */
	if(pse->variables[varid]->array == PSE_ARRAY) {
		switch(pse->variables[varid]->storage) {
		case PSE_VAR_INT:
			free(pse->variables[varid]->content.cint_a);
			pse->variables[varid]->content.cint_a = NULL;
			break;
		case PSE_VAR_DOUBLE:
			free(pse->variables[varid]->content.cdouble_a);
			pse->variables[varid]->content.cdouble_a = NULL;
			break;
		case PSE_VAR_TIME:
			free(pse->variables[varid]->content.ctime_a);
			pse->variables[varid]->content.ctime_a = NULL;
			break;
		case PSE_VAR_STRING:
			/*
			 * Initialize the structure and then all the strings associated to it.
			 * This allows constructing models that contain a variable number of
			 * strings of at most length 1000.
			 */
			for (i = 0; i < pse->variables[varid]->size; i++)
				free(pse->variables[varid]->content.cstring_a[i]);
			break;

			free(pse->variables[varid]->content.cstring_a);
			pse->variables[varid]->content.cstring_a = NULL;
		default:
			return PSE_ERROR_TYPE_UNKNOWN;
		}
	} else {
		/*
		 * Strings require initialization
		 */
		switch(pse->variables[varid]->storage) {
		case PSE_VAR_STRING:
			free(pse->variables[varid]->content.cstring);
			pse->variables[varid]->content.cstring = NULL;
			break;
		default:
			return PSE_ERROR_TYPE_UNKNOWN;
		}
	}

	free(pse->variables[varid]);
	pse->variables[varid] = NULL;
	pse->var_count--;

	return PSE_ERROR_OK;
}

/*
 * Add Bayesian dependencies (priors)
 */
pse_error pse_add_dependencies(pse_agent_stub *pse, pse_varid varid, int *conditionals,
																unsigned int count) {
	if (pse->state == CREATED)
			return PSE_ERROR_NOT_INITIALIZED;

	if (pse->state == STARTED)
		return PSE_ERROR_ALREADY_STARTED;

	if (pse->state == FINALIZED)
		return PSE_ERROR_ALREADY_FINALIZED;

	if (pse->variables[varid] == NULL)
		return PSE_ERROR_VARIABLE_UNKNOWN;

	if (pse->dependencies[varid] != NULL)
		return PSE_ERROR_DEPENDENCY_ALREADY_EXISTS;

	/*
	 * We only allow dependencies to variables from the world model. It is
	 * significant insofar joint dependencies imply simultaneity, which is not
	 * strictly defined in a purely relativistic universe. We do not rule out
	 * however an interpretation of belief in simultaneity, but that is an
	 * information-type event inside the agent, not an event in the world.
	 */
	if (pse_is_world_var(pse->variables[varid]) == PSE_FALSE)
		return PSE_ERROR_DEPENDENCY_NOT_WORLD;

	pse->dependencies[varid] = (pse_dependency *)malloc(sizeof(pse_dependency));
	pse->dependencies[varid]->count = count;
	memcpy(pse->dependencies[varid]->conditionals, conditionals, count*sizeof(int));
	pse->dependencies[varid]->priors = NULL;

	pse->variables[varid]->has_dependencies = PSE_TRUE;

	return PSE_ERROR_OK;
}

/*
 * Remove Bayesian depedencies.
 *
 * This function provides hooks for future implementations of advanced cognitive
 * models that change their internal representation by adding or removing
 * dependencies. For instance, if a certain dependency leads to small likelihood
 * of certain events, a bootstrapping procedure may help agents decide which
 * dependencies are artificial.
 */
pse_error pse_rm_dependencies(pse_agent_stub *pse, pse_varid varid) {
	if (pse->state == CREATED)
		return PSE_ERROR_NOT_INITIALIZED;

	if (pse->state == STARTED)
		return PSE_ERROR_ALREADY_STARTED;

	if (pse->state == FINALIZED)
		return PSE_ERROR_ALREADY_FINALIZED;

	if (pse->variables[varid] == NULL)
		return PSE_ERROR_VARIABLE_UNKNOWN;

	if (pse->dependencies[varid] == NULL)
		return PSE_ERROR_DEPENDENCY_UNKNOWN;

	free(pse->dependencies[varid]);

	pse->variables[varid]->has_dependencies = PSE_FALSE;

	return PSE_ERROR_OK;
}

/*
 * Supply the model for prior probabilities when dependencies exist. The
 * function can be constructed in terms of existing probability density
 * functions or can be supplied. In the future, a test for a function to be
 * an actual PDF should be required.
 *
 * We separate this part of the declaration of dependencies on grounds of
 * the intellectual complexity it may involve. As an example, this function
 * may serve as the connection to machine learning, classification or
 * a full-blown cognitive model.
 */
pse_error pse_supply_prior(pse_agent_stub * pse, pse_varid varid,
							double (*priors)(unsigned int, unsigned int *)) {
	if (pse->state == CREATED)
		return PSE_ERROR_NOT_INITIALIZED;

	if (pse->state == STARTED)
		return PSE_ERROR_ALREADY_STARTED;

	if (pse->state == FINALIZED)
		return PSE_ERROR_ALREADY_FINALIZED;

	if (pse->variables[varid] == NULL)
		return PSE_ERROR_VARIABLE_UNKNOWN;

	if (pse->dependencies[varid] == NULL)
		return PSE_ERROR_DEPENDENCY_UNKNOWN;

	pse->dependencies[varid]->priors = priors;

	return PSE_ERROR_OK;
}

/*
 * Use a variable as a template for another one. This is equivalent to the
 * second assignment step.
 */
pse_variable * pse_template(pse_variable *ptr_out, pse_variable *var) {
	int i;

	ptr_out = (pse_variable *) malloc(sizeof(pse_variable));
	ptr_out->storage = var->storage;
	ptr_out->model = var->model;
	ptr_out->locality = var->locality;
	ptr_out->point_distribution = var->point_distribution;
	ptr_out->array = var->array;
	ptr_out->size = var->size;
	memcpy(ptr_out->point_parameters, var->point_parameters, PSE_MAX_DIST_PARAMS*sizeof(double));
	ptr_out->has_dependencies = var->has_dependencies;
	ptr_out->read_and_alter = var->read_and_alter;
	ptr_out->array_distribution = var->array_distribution;
	memcpy(ptr_out->array_parameters, var->array_parameters, PSE_MAX_DIST_PARAMS*sizeof(double));
	strcpy(ptr_out->name, var->name);

	/*
	 * We process registration based on content type
	 */
	if(ptr_out->array == PSE_ARRAY) {
		switch(ptr_out->storage) {
		case PSE_VAR_INT:
			ptr_out->content.cint_a = (int *)malloc(sizeof(int)*ptr_out->size);
			break;
		case PSE_VAR_DOUBLE:
			ptr_out->content.cdouble_a = (double *)malloc(sizeof(double)*ptr_out->size);
			break;
		case PSE_VAR_TIME:
			ptr_out->content.ctime_a = (pse_time *)malloc(sizeof(pse_time)*ptr_out->size);
			break;
		case PSE_VAR_STRING:
			/*
			 * Initialize the structure and then all the strings associated to it.
			 * This allows constructing models that contain a variable number of
			 * strings of at most length 1000.
			 */
			ptr_out->content.cstring_a = (char **)malloc(sizeof(char *)*ptr_out->size);
			for (i = 0; i < ptr_out->size; i++)
				ptr_out->content.cstring_a[i] = (char *)malloc(sizeof(char)*PSE_MAX_STRLEN);
			break;
		default:
			break;
		}
	} else {
		/*
		 * Strings require initialization
		 */
		switch(ptr_out->storage) {
		case PSE_VAR_STRING:
			ptr_out->content.cstring = (char *)malloc(sizeof(char)*PSE_MAX_STRLEN);
			break;
		default:
			break;
		}
	}

	return ptr_out;
}

void pse_scratch(pse_variable *ptr_out) {
	if (ptr_out->array == PSE_ARRAY) {
		switch(ptr_out->storage) {
		case PSE_VAR_INT:
			free(ptr_out->content.cint_a);
			break;
		case PSE_VAR_DOUBLE:
			free(ptr_out->content.cdouble_a);
			break;
		case PSE_VAR_STRING:
			free(ptr_out->content.cstring_a);
			break;
		case PSE_VAR_TIME:
			free(ptr_out->content.ctime_a);
			break;
		default:
			return;
		}
	}

	free(ptr_out);

	return;
}

/*
 * Prepare the state of a variable
 *
 * A location is necessary when the data type is an array.
 * Input content types can only be scalars. This is a mechanism to avoid
 * 'cheating' when developing models that update large states. Therefore,
 * agents must update large states incrementally based on information they
 * are able to gather through communication.
 */
void pse_prepare(pse_agent_stub *pse, pse_varid varid, pse_content content,
				unsigned int location, pse_storage_type storage, pse_error *error) {
	pse_variable *p_to_var;

	if (pse->state == CREATED) {
		*error = PSE_ERROR_NOT_INITIALIZED;
		return;
	}

	if (pse->state == INITIALIZED) {
		*error = PSE_ERROR_NOT_INITIALIZED;
		return;
	}

	if (pse->state == FINALIZED) {
		*error = PSE_ERROR_ALREADY_FINALIZED;
		return;
	}

	if (pse->variables[varid] == NULL){
		*error = PSE_ERROR_VARIABLE_UNKNOWN;
		return;
	}

	p_to_var = pse->variables[varid];

	if (p_to_var->storage != storage) {
		*error = PSE_ERROR_TYPE_MISMATCH;
		return;
	}

	if (p_to_var->array == PSE_SCALAR) {
		switch(p_to_var->storage) {
		case PSE_VAR_INT:
			p_to_var->content.cint = content.cint;
			*error = PSE_ERROR_OK;
			break;
		case PSE_VAR_DOUBLE:
			p_to_var->content.cdouble = content.cdouble;
			*error = PSE_ERROR_OK;
			break;
		case PSE_VAR_STRING:
			strcpy(p_to_var->content.cstring, content.cstring);
			*error = PSE_ERROR_OK;
			break;
		case PSE_VAR_TIME:
			p_to_var->content.ctime = content.ctime;
			*error = PSE_ERROR_OK;
			break;
		default:
			*error = PSE_ERROR_TYPE_UNKNOWN;
			break;
		}
		return;
	} else {
		if (location >= p_to_var->size) {
			*error = PSE_ERROR_ARRAY_OUTOFBOUNDS;
			*error = PSE_ERROR_OK;
			return;
		}

		switch(p_to_var->storage) {
		case PSE_VAR_INT:
			p_to_var->content.cint_a[location] = content.cint;
			*error = PSE_ERROR_OK;
			break;
		case PSE_VAR_DOUBLE:
			p_to_var->content.cdouble_a[location] = content.cdouble;
			*error = PSE_ERROR_OK;
			break;
		case PSE_VAR_STRING:
			strcpy(p_to_var->content.cstring_a[location], content.cstring);
			*error = PSE_ERROR_OK;
			break;
		case PSE_VAR_TIME:
			p_to_var->content.cdouble_a[location] = content.cdouble;
			*error = PSE_ERROR_OK;
			break;
		default:
			*error = PSE_ERROR_TYPE_UNKNOWN;
			return;
		}
	}

	/*
	 * If a variable is stochastic, variation must be ensured.
	 */
	if (p_to_var->model == PSE_VAR_STOCHASTIC) {
		/*
		 * Separate by models that have dependencies.
		 */
		if (p_to_var->has_dependencies == PSE_FALSE) {
			/*
			 * We use a helper function to prepare the state with the desired
			 * distribution function(s).
			 */
			pse_randomize_and_alter(p_to_var, p_to_var, location, error);
			*error = PSE_ERROR_OK;
		} else {
			/*
			 * TODO: this needs to be implemented as a Bayesian distribution computation.
			 */
			*error = PSE_ERROR_OK;
		}
	} else {
		*error = PSE_ERROR_OK;
		return;
	}
}

/*
 * Observe function
 *
 * In the model, the equivalent of a read operation is an observe statement.For
 * deterministic variables, observations do not alter the state of the
 */
void pse_observe(pse_agent_stub *pse, pse_varid varid,
						unsigned int location, pse_variable *ptr_out, pse_error *error) {
	int csize;
	pse_variable *p_to_var;

	if (pse->state == CREATED) {
		*error = PSE_ERROR_NOT_INITIALIZED;
		return;
	}

	if (pse->state == INITIALIZED) {
		*error = PSE_ERROR_NOT_INITIALIZED;
		return;
	}

	if (pse->state == FINALIZED) {
		*error = PSE_ERROR_ALREADY_FINALIZED;
		return;
	}

	if (pse->variables[varid] == NULL){
		*error = PSE_ERROR_VARIABLE_UNKNOWN;
		return;
	}

	p_to_var = pse->variables[varid];

	if (p_to_var->model == PSE_VAR_DETERMINISTIC) {
		csize = pse_sizeof(p_to_var);

		if (csize != 0) {
			memcpy(ptr_out, &(p_to_var->content), csize);
			*error = PSE_ERROR_OK;
		} else {
			*error = PSE_ERROR_TYPE_UNKNOWN;
		}

		return;
	} else if (p_to_var->model == PSE_VAR_STOCHASTIC) {
		/*
		 * Separate by models that have dependencies.
		 */
		if (p_to_var->has_dependencies == PSE_FALSE) {
			if (p_to_var->read_and_alter == PSE_TRUE)
				pse_randomize_and_alter(ptr_out, p_to_var, location, error);
			else
				pse_randomize(ptr_out, p_to_var, location);

			*error = PSE_ERROR_OK;
		} else {
			/*
			 * TODO: this needs to be implemented as a Bayesian distribution computation.
			 */
			*error = PSE_ERROR_OK;

			return;
		}
	} else {
		*error = PSE_ERROR_OK;
		return;
	}
}

/*
 * Message to error logs depending on error type.
 */

void pse_error_log(pse_error error, char *buffer, char *arg) {
	char *final_arg =  (arg == NULL) ? "none" :  arg;
	switch(error) {
	/*
	 * No error should be printed.
	 * TODO: the PSE needs to have debug levels.
	 */
	case PSE_ERROR_OK:
		sprintf(buffer, PSE_ERROR_FMT, "Operation successful", final_arg);
		break;
	case PSE_ERROR_ALREADY_INITIALIZED:
		sprintf(buffer, PSE_ERROR_FMT, "The PSE has already been initialized", final_arg);
		break;
	case PSE_ERROR_ALREADY_FINALIZED:
		sprintf(buffer, PSE_ERROR_FMT, "The PSE has already been finalized", final_arg);
		break;
	case PSE_ERROR_ALREADY_STARTED:
		sprintf(buffer, PSE_ERROR_FMT, "The PSE has already been started", final_arg);
		break;
	case PSE_ERROR_NOT_INITIALIZED:
		sprintf(buffer, PSE_ERROR_FMT, "The PSE has not yet been initialized", final_arg);
		break;
	case PSE_ERROR_NOT_STARTED:
		sprintf(buffer, PSE_ERROR_FMT, "The PSE has not yet been started", final_arg);
		break;
	case PSE_ERROR_TOO_MANY_VARIABLES:
		sprintf(buffer, PSE_ERROR_FMT, "The PSE has too many variables registered in this agent", final_arg);
		break;
	case PSE_ERROR_VARIABLE_ALREADY_REGISTERED:
		sprintf(buffer, PSE_ERROR_FMT, "The PSE already contains this variable", final_arg);
		break;
	case PSE_ERROR_VARIABLE_UNKNOWN:
		sprintf(buffer, PSE_ERROR_FMT, "The PSE contains no such variable", final_arg);
		break;
	case PSE_ERROR_DEPENDENCY_ALREADY_EXISTS:
		sprintf(buffer, PSE_ERROR_FMT, "The PSE has contains this dependency", final_arg);
		break;
	case PSE_ERROR_DEPENDENCY_UNKNOWN:
		sprintf(buffer, PSE_ERROR_FMT, "The PSE has contains no such dependency", final_arg);
		break;
	case PSE_ERROR_DEPENDENCY_NOT_WORLD:
		sprintf(buffer, PSE_ERROR_FMT, "The dependency refers to an entity not in the world model", final_arg);
		break;
	case PSE_ERROR_TYPE_UNKNOWN:
		sprintf(buffer, PSE_ERROR_FMT, "The PSE recognizes no such type", final_arg);
		break;
	case PSE_ERROR_TYPE_MISMATCH:
		sprintf(buffer, PSE_ERROR_FMT, "Type mismatch for variable", final_arg);
		break;
	case PSE_ERROR_ARRAY_OUTOFBOUNDS:
		sprintf(buffer, PSE_ERROR_FMT, "Illegal out-of-bounds access of array attempted", final_arg);
		break;
	case PSE_ERROR_VARIABLE_IS_IMMUTABLE:
		sprintf(buffer, PSE_ERROR_FMT, "Illegal attempt to change immutable variable", final_arg);
		break;
	default:
		sprintf(buffer, PSE_ERROR_FMT, "Operation successful", final_arg);
		break;
	}

	return;
}

int pse_read_int(pse_agent_stub *pse, pse_varid varid) {
	return pse->variables[varid]->content.cint;
}

double pse_read_double(pse_agent_stub *pse, pse_varid varid) {
	return pse->variables[varid]->content.cdouble;
}

char * pse_read_string(pse_agent_stub *pse, pse_varid varid) {
	return pse->variables[varid]->content.cstring;
}

pse_time pse_read_time(pse_agent_stub *pse, pse_varid varid) {
	return pse->variables[varid]->content.ctime;
}
