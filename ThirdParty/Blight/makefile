CXX=g++
CC=gcc



#Here snippet should be your program
all: snippet bench_blight libblight.a
lib: libblight.a
#bench_blight is only here for benchmark purpose and can be removed



#Your compilation flags
CFLAGS_CUSTOM= -O9000  -std=c++11

#Other compilation flags
CFLAGS_LZ4+= -w -Wall -std=gnu99 -DUSE_THREADS  -fstrict-aliasing -Iext $(DEFS)
CFLAGS_BLIGHT+= -DNDEBUG -Ofast -flto -march=native -fPIC -mtune=native -g -std=c++11 -pipe -lz -fopenmp -msse4 -Ilz4

#Needed object files
LZ4O=lz4/lz4frame.o lz4/lz4.o lz4/xxhash.o lz4/lz4hc.o
BLO=blight.o utils.o



#Here  you should compile be your program using Blight instead of snippet
snippet: snippet.o $(BLO) $(LZ4O)
	$(CXX) -o $@ $^ $(CFLAGS_BLIGHT)

snippet.o: snippet.cpp
	$(CXX) -o $@ -c $< $(CFLAGS_CUSTOM)



#Benchmark compilation
bench_blight: bench_blight.o blight.o utils.o $(LZ4O)
	$(CXX) -o $@ $^ $(CFLAGS_BLIGHT)

bench_blight.o: bench_blight.cpp $(INC)
	$(CXX) -o $@ -c $< $(CFLAGS_BLIGHT)



#Blight object files (BLO)
utils.o: utils.cpp $(INC)
	$(CXX) -o $@ -c $< $(CFLAGS_BLIGHT) -fPIC

blight.o: blight.cpp $(INC)
	$(CXX) -o $@ -c $< $(CFLAGS_BLIGHT) -fPIC
	
	
	
#Lz4 object files (LZ4O)
lz4/lz4frame.o: lz4/lz4frame.c $(INC)
	$(CC) -o $@ -c $< $(CFLAGS_LZ4)

lz4/lz4.o: lz4/lz4.c $(INC)
	$(CC) -o $@ -c $< $(CFLAGS_LZ4)

lz4/lz4hc.o: lz4/lz4hc.c $(INC)
	$(CC) -o $@ -c $< $(CFLAGS_LZ4)
	
lz4/xxhash.o: lz4/xxhash.c $(INC)
	$(CC) -o $@ -c $< $(CFLAGS_LZ4) 


# libblight.a: blight.o utils.o lz4/lz4frame.o lz4/lz4.o lz4/xxhash.o lz4/lz4hc.o
# 	ar rcs libblight.a blight.o utils.o lz4/lz4frame.o lz4/lz4.o lz4/xxhash.o lz4/lz4hc.o
# 	ranlib libblight.a

# static library
libblight.a: $(BLO) $(LZ4O)
	gcc-ar rvs -o $@ $(BLO) $(LZ4O)

clean:
	rm -rf *.o
	rm -rf lib*a
	rm -rf snippet bench_blight
	rm -rf lz4/*o
	


rebuild: clean all
