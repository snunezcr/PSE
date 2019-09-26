/*
 * National Center for Supercomputing Applications
 * University of Illinois at Urbana-Champaign
 *
 * Large-Scale Agent-Based Social Simulation
 * Les Gasser, NCSA Fellow
 *
 * Author: Santiago Nunez-Corrales
 */

#include <psedict.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

/**
 * Initialize a dictionary
 * @param A dictionary
 * @return true if not previously initialized, false otherwise
 */
int pse_dict_init(pse_dictionary *dict) {
	if (dict->head != NULL)
		return 0;
	else {
		dict->head = NULL;
		dict->length = 0;
		return 1;
	}
}

/**
 * Finalize a dictionary
 * @param A dictionary
 */
void pse_dict_finalize(pse_dictionary *dict) {
	if (dict->head == NULL)
		return;
	else {
		pse_entry *to_delete = dict->head;
		pse_entry *next = NULL;

		do {
			next = to_delete->next;
			free(to_delete);
		} while (next != NULL);
	}
}

/**
 * Add an entry to the dictionary
 * @param A dictionary
 * @param A variable name
 * @param An index from the PSE corresponding to that variable registration
 * @return 1 if the addition was successful, 0 if it already exists
 */
int pse_dict_add(pse_dictionary *dict, char *variable, pse_varid pseval) {
	if (pse_dict_search(dict, variable) != -1)
		return 0;
	else if (dict->head == NULL) {
		dict->head = (pse_entry *) malloc(sizeof(pse_entry));
		strcpy(dict->head->varname, variable);
		dict->head->assigned = pseval;
		dict->head->next = NULL;
		return 1;
	}
	else {
		pse_entry *iterate = dict->head;

		while (iterate->next != NULL)
			iterate = iterate->next;

		// Once found, fill data in
		iterate->next = (pse_entry *) malloc(sizeof(pse_entry));
		strcpy(iterate->varname, variable);
		iterate->next->assigned = pseval;
		iterate->next->next = NULL;

		return 1;
	}
}

/**
 *
 * @param A dictionary
 * @param A variable name to look for and remove
 * @return 1 if remove was successful, 0 otherwise
 */
int pse_dict_remove(pse_dictionary *dict, char *variable) {
	if (dict->head == NULL)
		return 0;
	else {
		pse_entry *prev;
		pse_entry *curr;

		// Case 1: it is the head value
		if (strcmp(dict->head->varname, variable) == 0) {
			curr = dict->head;
			dict->head = curr->next;
			free(curr);
			return 1;
		} else {
			prev = dict->head;
			curr = dict->head->next;

			while (curr != NULL) {
				if (strcmp(curr->varname, variable) != 0) {
					prev = curr;
					curr = curr->next;
				} else {
					prev->next = curr->next;
					free(curr);
					return 1;
				}
			}

			return 0;
		}
	}
}

/**
 * Search an entry in the dictionary
 * @param A dictionary
 * @param A variable
 * @return -1 if the variable was not found, a value
 */
pse_varid pse_dict_search(pse_dictionary *dict, char *variable) {
	pse_entry *curr = dict->head;

	while (curr != NULL) {
		if (strcmp(curr->varname, variable) == 0)
			return curr->assigned;
	}

	return -1;
}

