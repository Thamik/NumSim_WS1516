VERSION = 0.1
CC 	= g++
CFLAGS 	= -Wall -DVERSION=\"$(VERSION)\" -std=c++11 -I src -pthread
LDFLAGS = 

OBJ 	= obj/main.o obj/geom_gen.o
BIN 	= geomgen_release geomgen_debug

release: $(OBJ)
	$(CC) $(CFLAGS) -o geomgen_release -O3 $(OBJ) $(LDFLAGS) -flto

debug: $(OBJ)
	$(CC) $(CFLAGS) -o geomgen_debug -g -O0 $(OBJ) $(LDFLAGS)

.PHONY: clean
clean:
	rm -rf $(BIN) $(OBJ)

obj/%.o: src/%.cpp
	$(CC) $(CFLAGS) -o $@ -c $<
