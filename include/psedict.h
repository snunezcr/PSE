/*
 * National Center for Supercomputing Applications
 * University of Illinois at Urbana-Champaign
 *
 * Large-Scale Agent-Based Social Simulation
 * Les Gasser, NCSA Fellow
 *
 * Author: Santiago Nunez-Corrales
 */

#include <pse.h>

/**
 * Structures for a PSE dictionary, to be handled as a singly-linked list.
 */
typedef struct pse_entry {
	char varname[PSE_VARNAME_SIZE];
	pse_varid assigned;
	struct pse_entry *next;
} pse_entry;

typedef struct pse_dictionary {
	int length;
	pse_entry *head;
} pse_dictionary;

/**
 * Functions to search inside the PSE dictionary. This provides a symbolic
 * link between PSE registration return values and variable names. This handle
 * is used by the OCaml compiler to make code access simpler and avoid generating
 * one temporary variable per call.
 */
int pse_dict_init(pse_dictionary *);
void pse_dict_finalize(pse_dictionary *);
int pse_dict_add(pse_dictionary *, char *, pse_varid);
int pse_dict_remove(pse_dictionary *, char *);
int pse_dict_search(pse_dictionary *, char *);
