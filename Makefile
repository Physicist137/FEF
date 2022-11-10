
CXX = g++
CXXFLAGS = --std=c++14
INCLUDE = -Iinclude -I/usr/include/eigen3 -Iexternal/SQLiteCpp/include
ERROR = -Wfatal-errors
WARNING = -Wall -Wextra 

SQLITE = -Lexternal/SQLiteCpp/build -lSQLiteCpp -l sqlite3
POSTGRESQL = -lpqxx -lpq
LIBS = -lquadmath -lboost_system $(SQLITE) $(POSTGRESQL)
NOTEBOOK = -fext-numeric-literals

SRC := $(wildcard src/*)
HPP := $(wildcard include/*)
OBJ := $(patsubst src/%.cpp, obj/%.o, $(SRC))
BIN := $(patsubst src/%.cpp, bin/%, $(SRC))
DEP := $(patsubst src/%.cpp, dep/%.d, $(SRC))

all: $(BIN)

bin/%: src/%.cpp
	@mkdir -p $(dir $@)
	@echo Creating $@...
	@$(CXX) $(CXXFLAGS) $< -o $@ $(INCLUDE) $(LIBS) $(WARNING) $(ERROR) $(NOTEBOOK) -O3 

#obj/%.o: src/%.cpp
#	@mkdir -p $(dir $@)
#	@echo Creating $@...
#	@$(CXX) $(CXXFLAGS) -c $< -o $@ $(INCLUDE) $(LIBS) $(WARNING)

#dep/%.d: src/%.cpp $(HPP)
#	@mkdir -p $(dir $@)
#	@$(CXX) $(CXXFLAGS) $< -MM -MG -MP -MF$@ -MT$@ -MT$(@:dep/%.d=obj/%.o) $(INCLUDE) $(LIBS) $(WARNING)

#-include $(DEP)


.PHONY: clear

clear:
	@rm -rf bin
	@rm -rf dep
	@rm -rf obj


