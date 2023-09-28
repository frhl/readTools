CC = g++
CFLAGS = -Wall -O2
INCLUDES = -I/usr/local/include
LIBS = -L/usr/local/lib -lhts -lz
TARGET1 = candidates
TARGET2 = phaseRead
SOURCE1 = candidates.cpp
SOURCE2 = phaseRead.cpp

all: $(TARGET1) $(TARGET2)

$(TARGET1): $(SOURCE1)
		$(CC) $(CFLAGS) $(INCLUDES) $(SOURCE1) -o $(TARGET1) $(LIBS)

$(TARGET2): $(SOURCE2)
		$(CC) $(CFLAGS) $(INCLUDES) $(SOURCE2) -o $(TARGET2) $(LIBS)

clean:
		rm -f $(TARGET1) $(TARGET2) 

