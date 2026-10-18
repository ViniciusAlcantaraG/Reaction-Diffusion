CXX = g++
CXXFLAGS = -O3 -ffast-math -march=native -flto -Wall -std=c++17
LDFLAGS = -lraylib -lgdi32 -lwinmm

TARGET = morphogenesis.exe
SRC = main.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET) $(LDFLAGS)

clean:
	rm -f $(TARGET)