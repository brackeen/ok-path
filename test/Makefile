CC = clang -Wall -O3 -std=c99
BIN = ./bin
INCLUDES = -I .. 
FILES = *.c ../*.c
OUT_NAME = ok_path_test

test: build
	$(BIN)/$(OUT_NAME)
	
build: init  
	$(CC) -o $(BIN)/$(OUT_NAME) -Wall $(INCLUDES) $(FILES)

init:
	@mkdir -p $(BIN)
	
clean:
	rm -rf $(BIN)
