CC := gcc
SRC := src
PROJECT_NAME := main


SRC_FILES := $(shell find $(SRC) -name '*.c')

inc := -I inc -lgsl  -lm 

CFLAGS= -Wall

all:
	$(CC) $(CFLAGS) $(inc) $(SRC_FILES)  -o $(PROJECT_NAME) 

.PHONY: clean
clean:
	rm -rf $(PROJECT_NAME)