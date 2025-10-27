CXX = g++
CXXFLAGS = -std=c++17 -Iinclude -I/opt/homebrew/include -Wall -Wextra -O2
LDFLAGS = -L/opt/homebrew/lib -lboost_program_options

SRCS = src/main.cpp src/io.cpp src/probes.cpp
OBJS = $(SRCS:.cpp=.o)
TARGET = probeDesign

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(TARGET) $(OBJS)