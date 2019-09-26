/*
 * National Center for Supercomputing Applications
 * University of Illinois at Urbana-Champaign
 *
 * Large-Scale Agent-Based Social Simulation
 * Les Gasser, NCSA Fellow
 *
 * Author: Santiago Nunez-Corrales
 */

#define PSE_MAX_VARIABLES 	2000
#define PSE_VARNAME_SIZE 	50
#define PSE_MAX_STRLEN 		1000
#define PSE_MAX_DIST_PARAMS	5
#define PSE_TRUE 			1
#define PSE_FALSE			0
#define PSE_HEADS 			1
#define PSE_TAILS			0

/*
 * Fixed strings
 */
#define PSE_ERROR_FMT	"[PSE Runtime] %s. Argument: %s.\n"

typedef enum pse_storage_type {
	PSE_VAR_INT,
	PSE_VAR_DOUBLE,
	PSE_VAR_STRING,
	PSE_VAR_TIME
} pse_storage_type;

typedef enum pse_array_type {
	PSE_SCALAR,
	PSE_ARRAY
} pse_array_type;

typedef enum pse_model_type {
	PSE_VAR_STOCHASTIC,
	PSE_VAR_DETERMINISTIC
} pse_model_type;

/*
 * Distributions come in two flavors:
 * - Those that do not use the current value as input
 * - Those that do (SELF).
 */

typedef enum pse_distribution_type {
	PSE_DIST_UNIFORM_INT_SELF,
	PSE_DIST_UNIFORM_INT_BOUNDED,
	PSE_DIST_BERNOULLI,
	PSE_DIST_BINOMIAL,
	PSE_DIST_BINOMIAL_SELF,
	PSE_DIST_NEG_BINOMIAL,
	PSE_DIST_NEG_BINOMIAL_SELF,
	PSE_DIST_POISSON,
	PSE_DIST_POISSON_SELF,
	PSE_DIST_UNIFORM_DOUBLE_SELF,
	PSE_DIST_UNIFORM_DOUBLE_BOUNDED,
	PSE_DIST_NORMAL,
	PSE_DIST_NORMAL_SELF,
	PSE_DIST_EXPONENTIAL,
	PSE_DIST_EXPONENTIAL_SELF,
	PSE_DIST_GAMMA,
	PSE_DIST_GAMMA_SELF,
	PSE_DIST_CHISQ,
	PSE_DIST_CHISQ_SELF,
	PSE_DIST_F,
	PSE_DIST_BETA,
	PSE_DIST_FOKKER_PLANCK,
	PSE_DIST_CUSTOM,
	PSE_DIST_NONE
} pse_distribution_type;

/*
 * Should there be a MESSAGE locality type?
 *
 * Message internals would be accesible outside some way.
 */
typedef enum pse_locality_type {
	PSE_AGENT,
	PSE_WORLD
} pse_locality_type;

/*
 * In this models, variables are initialized as needed. No relocation is made
 * when a variable is deleted. In general, preliminary models will not require
 * de-registering variables. More importantly, there is no need for concurrency
 * since agents have a local PSE stub in their Charm++/ROSS definition.
 */

typedef int pse_varid;
typedef int pse_depid;

/* We assume that Bayes rule applies in the form of
 *
 * P(y | x1, x2, x3, ... , xn) = .
 *
 * That is, computing Bayesian simultaneity of a posteriori events is invalid.
 *
 * Also, computing is performed in double floating-point arithmetic and recast
 * to integer when needed. The prior probability function is needed when a
 * dependency is indicated.
 */

typedef struct pse_dependency {
	unsigned int count;
	pse_depid conditionals[PSE_MAX_VARIABLES];
	double (*priors)(unsigned int, unsigned int *);
} pse_dependency;


/*
 * At present the time representation within Charm++/ROSS is equivalent to a double.
 * This is not ideal since other frameworks may differ. The code below is the point
 * of contact for establishing a time representation.
 */
typedef double pse_time;

/*
 * A PSE variable is an object that can be measured with respect to a prior
 * observed value and a set of dependencies.
 */

typedef union pse_content {
	int cint;
	double cdouble;
	pse_time ctime;
	char *cstring;
	int *cint_a;
	double *cdouble_a;
	pse_time *ctime_a;
	char **cstring_a;
} pse_content;

/*
 * Definition of variables
 *
 * There are two distributions per variable depending on the array type.
 * The first distribution is concerned with how data vary at a point location
 * in scalar fashion. The second distribution, if the array is present, drives
 * the location in an array where variation occurs.
 *
 * An interesting flag is read_and_alter. This is a destructive operation in
 * the sense in which measurements modify the content of a variable. If active,
 * each observe call replaces the value with the most recent stochastic one.
 */
typedef struct pse_variable {
	pse_storage_type storage;
	pse_model_type model;
	pse_locality_type locality;
	pse_distribution_type point_distribution;
	double point_parameters[PSE_MAX_DIST_PARAMS];
	unsigned int has_dependencies;
	unsigned int read_and_alter;
	char name[PSE_VARNAME_SIZE];
	pse_content content;
	pse_array_type array;
	unsigned int size;
	pse_distribution_type array_distribution;
	double array_parameters[PSE_MAX_DIST_PARAMS];
} pse_variable;

typedef enum pse_state {
	CREATED,
	INITIALIZED,
	STARTED,
	FINALIZED
} pse_state;

/*
 * Var count and var limit differ in terms of what has been used in the array
 * and how many variables are used.
 */
typedef struct pse_agent_stub {
	pse_state state;
	unsigned int var_count;
	unsigned int var_limit;
	pse_variable *variables[PSE_MAX_VARIABLES];
	pse_dependency *dependencies[PSE_MAX_VARIABLES];
} pse_agent_stub;


typedef enum pse_error {
	PSE_ERROR_OK 							= 0,
	PSE_ERROR_ALREADY_INITIALIZED 			= -1,
	PSE_ERROR_ALREADY_FINALIZED 				= -3,
	PSE_ERROR_ALREADY_STARTED				= -4,
	PSE_ERROR_NOT_INITIALIZED 				= -5,
	PSE_ERROR_NOT_STARTED	 				= -6,
	PSE_ERROR_TOO_MANY_VARIABLES 			= -7,
	PSE_ERROR_VARIABLE_ALREADY_REGISTERED 	= -9,
	PSE_ERROR_VARIABLE_UNKNOWN				= -11,
	PSE_ERROR_DEPENDENCY_ALREADY_EXISTS		= -13,
	PSE_ERROR_DEPENDENCY_UNKNOWN				= -15,
	PSE_ERROR_DEPENDENCY_NOT_WORLD			= -16,
	PSE_ERROR_TYPE_UNKNOWN					= -17,
	PSE_ERROR_TYPE_MISMATCH					= -19,
	PSE_ERROR_ARRAY_OUTOFBOUNDS				= -21,
	PSE_ERROR_VARIABLE_IS_IMMUTABLE			= -23
} pse_error;

/*
 * Function prototypes. The operation of the machine is simplified through
 * a set of complementary functions. No arithmetics is performed inside the
 * PSE, which is consistent with the fact that agents make decisions based on
 * measurements and information towards having a cognitive architecture, and
 * not through the direct manipulation of the laws of physics of their
 * respective worlds. Templating is a second step in variable declaration
 * for our model. A variable is unusable if after its declaration it has no
 * structure (hence, all Charm++/ROSS variables will really be pointers).
 * When it is not necessary, needs to be discarded (freed and set to NULL).
 */

pse_error pse_init(pse_agent_stub *);
pse_error pse_start(pse_agent_stub *, int, int);
pse_error pse_finalize(pse_agent_stub *);

pse_varid pse_register(pse_agent_stub *, pse_storage_type, pse_model_type,
						pse_locality_type, pse_distribution_type, double *,
						pse_array_type, unsigned int, unsigned int,
						pse_distribution_type, double *, char *);
pse_error pse_deregister(pse_agent_stub *, pse_varid);

pse_error pse_add_dependencies(pse_agent_stub *, pse_varid, int *, unsigned int);
pse_error pse_rm_dependencies(pse_agent_stub *, pse_varid);

pse_error pse_supply_prior(pse_agent_stub *, pse_varid,
							double (*priors)(unsigned int, unsigned int *));


pse_variable * pse_template(pse_variable *, pse_variable *);
void pse_scratch(pse_variable *);

void pse_prepare(pse_agent_stub *, pse_varid, pse_content, unsigned int,
						pse_storage_type,pse_error *);
void pse_observe(pse_agent_stub *, pse_varid, unsigned int, pse_variable *, pse_error *);

void pse_error_log(pse_error, char *, char *);

int pse_read_int(pse_agent_stub *, pse_varid);
double pse_read_double(pse_agent_stub *, pse_varid);
char * pse_read_string(pse_agent_stub *, pse_varid);
pse_time pse_read_time(pse_agent_stub *, pse_varid);
