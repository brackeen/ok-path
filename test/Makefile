CC = clang -Wall -g -std=c99
BIN = ./bin
INCLUDES = -I .. 
FILES = *.c ../*.c
OUT_NAME = ok_path_test
HAS_VALGRIND := $(shell type valgrind 1>/dev/null 2>&1 ; echo $$?)

test: build
ifeq ($(HAS_VALGRIND),0)
	valgrind --tool=memcheck --leak-check=full $(BIN)/$(OUT_NAME)
else
	$(BIN)/$(OUT_NAME)
endif
	
build: init  
	$(CC) -o $(BIN)/$(OUT_NAME) -Wall $(INCLUDES) $(FILES)

init:
	@mkdir -p $(BIN)
	
clean:
	rm -rf $(BIN)
