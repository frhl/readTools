CC = g++
CFLAGS = -Wall -O2
INCLUDES = -I/usr/local/include
LIBS = -L/usr/local/lib -lhts -lz
TARGET1 = pairFinder
TARGET2 = vcfPairFinder
TARGET3 = readOverlap
SOURCE1 = pairFinder.cpp
SOURCE2 = vcfPairFinder.cpp
SOURCE3 = readOverlap.cpp

all: $(TARGET1) $(TARGET2) $(TARGET3)

$(TARGET1): $(SOURCE1)
		$(CC) $(CFLAGS) $(INCLUDES) $(SOURCE1) -o $(TARGET1) $(LIBS)

$(TARGET2): $(SOURCE2)
		$(CC) $(CFLAGS) $(INCLUDES) $(SOURCE2) -o $(TARGET2) $(LIBS)

$(TARGET3): $(SOURCE3)
		$(CC) $(CFLAGS) $(INCLUDES) $(SOURCE3) -o $(TARGET3) $(LIBS)


clean:
		rm -f $(TARGET1) $(TARGET2) $(TARGET3) 

