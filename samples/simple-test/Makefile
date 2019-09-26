# National Center for Supercomputing Applications
# University of Illinois at Urbana-Champaign
# 
# Large-Scale Agent-Based Social Simulation
# Les Gasser, NCSA Fellow
   
# Author: Santiago Nunez-Corrales
BASE_DIR=../../../src
RAND_DIR=$(BASE_DIR)/rand
PSE_DIR=$(BASE_DIR)
INCLUDE_DIR=$(BASE_DIR)/include
TEST_NAME=01-simple-pse

CFLAGS=-Wall
LDFLAGS=-I$(INCLUDE_DIR)

all:
	@echo "Building test application $(TEST_NAME)..."
	@gcc $(CFLAGS) $(TEST_NAME).c $(PSE_DIR)/pse.c $(RAND_DIR)/ranlib.c $(RAND_DIR)/rnglib.c $(LDFLAGS) -o $(TEST_NAME)
	@echo "Done."
	
clean:
	@echo "Cleaning build for $(TEST_NAME)..."
	@rm $(TEST_NAME)
	@echo "Done."
