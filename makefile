
CXX = g++-14
 
CXXFLAGS = -std=c++17 -Wall -fopenmp 


SRC_DIR = .
BIN_DIR = bin


TARGET = $(BIN_DIR)/1


SRCS = $(SRC_DIR)/conservationform.C \
	   $(SRC_DIR)/initstate.C \
	   $(SRC_DIR)/Solver.C \
	   $(SRC_DIR)/Solver_HLLC.C \
	   $(SRC_DIR)/Parallel.C \
       $(SRC_DIR)/main.C


OBJS = $(SRCS:.C=.o)


all: $(TARGET)


$(TARGET): $(OBJS)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)


%.o: %.C
	$(CXX) $(CXXFLAGS) -c $< -o $@


clean:
	rm -f $(OBJS)
	rm -f $(TARGET)


clean_exec:
	rm -f $(TARGET)

clean_objs:
	rm -f $(OBJS)

.PHONY: all clean clean_exec clean_objs
