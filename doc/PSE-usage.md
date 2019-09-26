# Polarizing Stochastic Engine Developer Manual

The Polarizing Stochastic Engine is a standalone software library that can be
used to simulate stochastic variables. A PSE is similar in features to a 
probabilistic memory, but with fine-grained control. We briefly explain how to
build programs using the PSE in this document.

## Includes

Using the PSE library requires importing *pse.h* in the file where needed.

```c
#include <pse.h>
```

This file provides all required constructs for using the PSE, including error
reporting:

```c
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
```

All functions that alter the data structures in the PSE return an error code
representing the result of the operation. Some functions that manage data
contents also return error codes. In addition, a log function generates strings
accordingly.

```c
	pse_error_log(errno, errmsg, "init");
```

## PSE stubs

A PSE is a *stub* of memory with dynamic allocation properties, including the
ability to register and remove variables at runtime. They should be declared at
the beginning of the respective *.c file.

```c
	pse_agent_stub test_pse;
```

Stub contain all aspects of the state of the PSE and should be manipulated by
means of the functions provided in the API. Direct manipulation is discouraged
for reasons of consistency. PSE's can be:

- initialized

```c
	errno = pse_init(&test_pse);
```

- started

```c
	errno = pse_start(&test_pse,rnd1,rnd2);
```

where rnd1 and rnd2 are random seeds. These seeds may be utilized in 
distributed settings to control randomness across the execution basis.

- finalized

```c
	errno = pse_finalize(&test_pse);
```

## Data contents

Data access from and to the PSE requires a transfer data structure known as the
content handler. The number of content handlers depends on the use case.

```c
	pse_content temp_content;
```

Internally, the *pse_content* structure is a C union as follows:

```c
typedef union pse_content {
	int cint;
	double cdouble;
	time ctime;
	char *cstring;
	int *cint_a;
	double *cdouble_a;
	time *ctime_a;
	char **cstring_a;
} pse_content;
```

## Variables in the PSE

### Internal vs temporary variables
Variables belong to two different classes in the PSE:

* Variables that are manipulated indirectly inside the PSE
* Variables that are manipulated directly outside the PSE

Variables inside the PSE need to be registered and prepared, and may be 
observed to obtain a measurement of their value. Their measurement depends on
various factors:

* If the variable is deterministic, it will return the stored valued without
  alterations.
* If the variable is stochastic but not *self-updating*, the stored value will
  not change, but the results will. This variable is useful to represent
  identically distributed random variables with a given distribution.
* If the variable is stochastic and self-updating, the stored value will change
  as specififed by the probability distribution selected at its registration
  time.

### Registering variables

Variables in the PSE are not *declared*, but rather *registered*. The semantics
behind this process is a constructive one: the abstract model of the machine is
built for particular cases and, only when it is defined, can be used.

The result of registering a variable is a unique unsigned integer identificator.
Variables then need to be accessed using this quantity. It is advisable not to
handle variable identificators manually, but to use the dictionary library that
is provided with this code.

```c
	point_params[0] = 50.0;
	point_params[1] = 2.3;

	array_params[0] = 0.0;

	varid_double = pse_register(&test_pse, PSE_VAR_DOUBLE, PSE_VAR_STOCHASTIC,
								PSE_AGENT, PSE_DIST_NORMAL_SELF, point_params,
								PSE_SCALAR, 1, PSE_TRUE, PSE_DIST_NONE,
								array_params, "distance");
```
The code above registers a variable named *distance* that is stochastic, is
stored in the internal representation of an agent and is normally distributed.
Moreover, it is self-updating and of scalar type (hence of size 1). Point 
parameters are parameters for scalar distributions, and array parameters are
applied to select elements that need to be selected from an array before the
application of point paramters. This provides a level of fine control over what
changes in a stochastic program.

### Preparing variables

All variables in the PSE need to be prepared before usage. Accessing an
unprepared variable leads to undefined behavior by design. Data types need to
be encapsulated with a generic union type that simulates polymorphism.

```c
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
```

For the previous case, the prepare statement proceeds as follows:

```c
	pse_prepare(&test_pse, varid_double, temp_content, 0, PSE_VAR_DOUBLE, &errno);
```

### Observing variables

In order to interact with code outside the PSE, variables need to be observed.
The latter means that, depending on the properties specified during allocation
of the variable, it may be altered upon reading it. Similar to other virtual
memory devices, programs need to define temporary scratch space to handle data.

```c
	pse_variable *temp_var = NULL;
```

That is not sufficent, however, to copy the data back. Because of the simulated
polymorphism involved in data types that can hold scalar/array data types, and
the integer/floating/time/string contents, variables need to be templated for
particular uses, as manner of registers.

```c
	temp_var = pse_template(temp_var,test_pse.variables[varid_double]);
```

Once the variable has been templated, it is ready to hold data that will be 
converted into regular variables outside the PSE.

```c
	pse_observe(&test_pse, varid_double, 0, temp_var, &errno);
```

Notice that *varid_double* is of *pse_varid* type. For convenience, we provide
a dictionary library (implemented as a singly-linked list) that uses strings
representing variable names as keys and pse_varid's as values. The latter was
designed to fit the needs of the Social Theory Scaling Environment (STSE), as
keeping track of possibly many variable names in a simulation would be a 
dauting task only by using memory locations inside the PSE.

### Instrumentation

All data types have associated a *pse_read_X* function where X is the data 
type. These functions are to be used only for instrumentation purposes and do
not replace calls to *pse_observe*. 