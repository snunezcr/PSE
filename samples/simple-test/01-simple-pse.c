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
#include <stdio.h>
#include <pse.h>

#define ERROR_BUFF_SIZE 200
#define TEST_CYCLES	30

/*
 * Purpose of the test:
 * --------------------
 *
 * Show basic functionality of the PSE with variable
 * registration and manipulation.
 */
int main(int argc, char **argv) {
	/*
	 * Temporary reference data structures
	 */
	pse_agent_stub test_pse;
	pse_variable *temp_var = NULL;
	pse_content temp_content;

	/*
	 * Distribution parameters
	 */
	double point_params[PSE_MAX_DIST_PARAMS] = {0.0,0.0,0.0,0.0,0.0};
	double array_params[PSE_MAX_DIST_PARAMS] = {0.0,0.0,0.0,0.0,0.0};

	/*
	 * Variable id pointers. This needs to be replaced with a translation map
	 * for TinySocial where variable names are the keys.
	 */
	int varid_double;
	int varid_int;
	int varid_time;

	int i;

	/*
	 * Error handling
	 */
	pse_error errno;
	char errmsg[ERROR_BUFF_SIZE];

	/*
	 * Test 1: spurious start
	 */
//	errno = pse_start(&test_pse,103,29);
//	pse_error_log(errno, errmsg, NULL);
//	fprintf(stderr, "%s", errmsg);

	/*
	 * Initialize the PSE
	 */
	errno = pse_init(&test_pse);
	pse_error_log(errno, errmsg, "init");
	fprintf(stderr, "%s", errmsg);

	/*
	 * Register test variables
	 */

	point_params[0] = 50.0;
	point_params[1] = 2.3;

	array_params[0] = 0.0;

	varid_double = pse_register(&test_pse, PSE_VAR_DOUBLE, PSE_VAR_STOCHASTIC,
								PSE_AGENT, PSE_DIST_NORMAL_SELF, point_params,
								PSE_SCALAR, 1, PSE_TRUE, PSE_DIST_NONE,
								array_params, "distance");
	pse_error_log(varid_double, errmsg, "distance");
	fprintf(stderr, "%s", errmsg);

	point_params[0] = 100.0;
	point_params[1] = 0.5;
	array_params[0] = 0.0;

	varid_int = pse_register(&test_pse, PSE_VAR_INT, PSE_VAR_STOCHASTIC,
							PSE_AGENT, PSE_DIST_BINOMIAL, point_params,
							PSE_SCALAR, 1, PSE_FALSE, PSE_DIST_NONE,
							array_params, "hopping_steps");
	pse_error_log(varid_int, errmsg, "hopping_steps");
	fprintf(stderr, "%s", errmsg);

	point_params[0] = 0.0;
	array_params[0] = 0.0;

	varid_time = pse_register(&test_pse, PSE_VAR_TIME, PSE_VAR_DETERMINISTIC,
							 PSE_AGENT, PSE_DIST_NONE, point_params,
							 PSE_SCALAR, 1, PSE_FALSE, PSE_DIST_NONE,
							 array_params, "deterministic_time");
	pse_error_log(varid_time, errmsg, "deterministic_time");
	fprintf(stderr, "%s", errmsg);

	/*
	 * Non spurious start
	 */
	errno = pse_start(&test_pse,103,29);
	pse_error_log(errno, errmsg, "start");
	fprintf(stderr, "%s", errmsg);

	/*
	 * Prepare variables
	 */
	temp_content.cdouble = 12.4;
	pse_prepare(&test_pse, varid_double, temp_content, 0, PSE_VAR_DOUBLE, &errno);
	pse_error_log(errno, errmsg, "prepare distance");
	fprintf(stderr, "%s", errmsg);

	temp_content.cint = 130;
	pse_prepare(&test_pse, varid_int, temp_content, 0, PSE_VAR_INT, &errno);
	pse_error_log(errno, errmsg, "hopping_steps");
	fprintf(stderr, "%s", errmsg);

	temp_content.ctime = 54.5;
	pse_prepare(&test_pse, varid_time, temp_content, 0, PSE_VAR_TIME, &errno);
	pse_error_log(errno, errmsg, "prepare deterministic_time");
	fprintf(stderr, "%s", errmsg);

	/*
	 * Observe double variables iteratively
	 */
	for (i = 0; i < TEST_CYCLES; i++) {
		temp_var = pse_template(temp_var,test_pse.variables[varid_double]);

		pse_observe(&test_pse, varid_double, 0, temp_var, &errno);
		pse_error_log(errno, errmsg, "observe distance");
		fprintf(stderr, "%s", errmsg);

		if (errno == PSE_ERROR_OK) {
			printf("[PSE Runtime] Iteration: %d\tRead double (return) value: %lf\n", i,
					temp_var->content.cdouble);
			printf("[PSE Runtime] Iteration: %d\tRead double (stored) value: %lf\n", i,
					pse_read_double(&test_pse, varid_double));
		}

		pse_scratch(temp_var);
	}

	/*
	 * Observe int variables iteratively
	 */
	for (i = 0; i < TEST_CYCLES; i++) {
		temp_var = pse_template(temp_var,test_pse.variables[varid_int]);

		pse_observe(&test_pse, varid_int, 0, temp_var, &errno);
		pse_error_log(errno, errmsg, "observe hopping_step");
		fprintf(stderr, "%s", errmsg);

		if (errno == PSE_ERROR_OK) {
			printf("[PSE Runtime] Iteration: %d\tRead int (return) value: %d\n", i,
					temp_var->content.cint);
			printf("[PSE Runtime] Iteration: %d\tRead int (stored) value: %d\n", i,
					pse_read_int(&test_pse, varid_int));
		}

		pse_scratch(temp_var);
	}

	/*
	 * Observe time variables iteratively
	 */
	for (i = 0; i < TEST_CYCLES; i++) {
		temp_var = pse_template(temp_var,test_pse.variables[varid_time]);

		pse_observe(&test_pse, varid_time, 0, temp_var, &errno);
		pse_error_log(errno, errmsg, "observe deterministic_time");
		fprintf(stderr, "%s", errmsg);

		if (errno == PSE_ERROR_OK)
			printf("[PSE Runtime] Iteration: %d\tRead time value: %lf\n", i,
					pse_read_time(&test_pse, varid_time));

		pse_scratch(temp_var);
	}

	/*
	 * Finalize the PSE
	 */
	errno = pse_finalize(&test_pse);
	pse_error_log(errno, errmsg, "finalize");
	fprintf(stderr, "%s", errmsg);

	return PSE_ERROR_OK;
}
