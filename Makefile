CXX = g++
CXXFLAGS = -std=c++17 -Iinclude -Wall -Wextra -O2

SRCS = src/main.cpp src/io.cpp src/probes.cpp
TARGET = probeDesign

all: $(TARGET)

$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) -o $@ $(SRCS)

clean:
	rm -f $(TARGET) *.o